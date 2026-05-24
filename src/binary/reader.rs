//! Deterministic reader for `.kira-spatial.bin` (format v2).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crc32fast::Hasher as Crc32;
use serde_json::Value;

use crate::api::dataset::Dataset;
use crate::binary::bitmask::tissue_mask_from_u64_words;
use crate::binary::compress::zstd_decompress;
use crate::binary::format::{
    CSR_FLAG_INDPTR_U32, HEADER_SIZE, KIRA_SPATIAL_BIN_VERSION, MAGIC, MANDATORY_SECTION_IDS,
    MAX_SECTION_COUNT, MIN_SECTION_COUNT, SECTION_ENTRY_SIZE, SECTION_FLAG_ZSTD, SECTION_ID_CSR,
    SECTION_ID_FEATURE_TABLE, SECTION_ID_META_CORE, SECTION_ID_META_JSON,
    SECTION_ID_SPATIAL_DOMAIN, SectionEntry,
};
use crate::binary::hash::compute_dataset_hash;
use crate::config::LoadConfig;
use crate::determinism::float::{ensure_f32_finite_nonneg, total_cmp_f32};
use crate::determinism::json::{canonicalize_json, write_canonical_json};
use crate::error::SpatialIoError;
use crate::model::{
    coord::CoordSystem,
    csr::{BinsCsr, Indptr},
    features::{FeatureRow, FeatureTable},
    metadata::DatasetMetaCore,
    spatial_domain::SpatialDomain,
};

/// Reads and validates a deterministic `.kira-spatial.bin` dataset using the default config.
pub fn read_kira_bin<P: AsRef<Path>>(p: P) -> Result<Dataset, SpatialIoError> {
    read_kira_bin_with(p, &LoadConfig::default())
}

/// Reads a dataset with explicit [`LoadConfig`] knobs (hash validation toggle, memory cap).
pub fn read_kira_bin_with<P: AsRef<Path>>(
    p: P,
    cfg: &LoadConfig,
) -> Result<Dataset, SpatialIoError> {
    let path = p.as_ref();
    let mut file = File::open(path).map_err(|e| SpatialIoError::io_at(path, e))?;
    let file_meta = file.metadata().map_err(|e| SpatialIoError::io_at(path, e))?;
    let file_len = file_meta.len();

    let budget_bytes = (cfg.memory_budget_mb as u64)
        .checked_mul(1024 * 1024)
        .unwrap_or(u64::MAX);
    if file_len > budget_bytes {
        return Err(SpatialIoError::MemoryLimitExceeded(format!(
            "file size {file_len} exceeds memory budget {budget_bytes}"
        )));
    }

    if file_len < HEADER_SIZE {
        return Err(SpatialIoError::UnsupportedFormat(
            "file too small for header".to_string(),
        ));
    }

    let cap = usize::try_from(file_len)
        .map_err(|_| SpatialIoError::UnsupportedFormat("file too large for address space".to_string()))?;
    let mut bytes_vec = Vec::with_capacity(cap);
    file.read_to_end(&mut bytes_vec).map_err(|e| SpatialIoError::io_at(path, e))?;
    let bytes: &[u8] = &bytes_vec;

    if &bytes[0..8] != MAGIC.as_slice() {
        return Err(SpatialIoError::UnsupportedFormat(
            "invalid magic: expected KIRASPAT".to_string(),
        ));
    }

    let version = u16_at(bytes, 8)?;
    if version != KIRA_SPATIAL_BIN_VERSION {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "unsupported version: {version} (this build supports v{KIRA_SPATIAL_BIN_VERSION})"
        )));
    }

    let section_count = u16_at(bytes, 10)?;
    if section_count < MIN_SECTION_COUNT {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "section_count too small: {section_count}"
        )));
    }
    if section_count > MAX_SECTION_COUNT {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "section_count exceeds cap: {section_count} > {MAX_SECTION_COUNT}"
        )));
    }

    let mut header_hash = [0_u8; 16];
    header_hash.copy_from_slice(
        bytes
            .get(16..32)
            .ok_or_else(|| SpatialIoError::UnsupportedFormat("missing header hash".to_string()))?,
    );

    let table_start = HEADER_SIZE as usize;
    let table_len = (section_count as usize)
        .checked_mul(SECTION_ENTRY_SIZE as usize)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat("section table size overflow".to_string())
        })?;
    let table_end = table_start.checked_add(table_len).ok_or_else(|| {
        SpatialIoError::UnsupportedFormat("section table end overflow".to_string())
    })?;
    if table_end as u64 > file_len {
        return Err(SpatialIoError::UnsupportedFormat(
            "section table out of file bounds".to_string(),
        ));
    }

    let mut sections = Vec::with_capacity(section_count as usize);
    for i in 0..section_count as usize {
        let base = table_start + i * SECTION_ENTRY_SIZE as usize;
        let id = u16_at(bytes, base)?;
        let flags = u16_at(bytes, base + 2)?;
        let offset = u64_at(bytes, base + 4)?;
        let length = u64_at(bytes, base + 12)?;
        let crc32c = u32_at(bytes, base + 20)?;

        if length > 0 && offset % 64 != 0 {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "section {id} offset is not 64-byte aligned: {offset}"
            )));
        }

        if length > 0 {
            let end = offset.checked_add(length).ok_or_else(|| {
                SpatialIoError::UnsupportedFormat(format!("section {id} end overflow"))
            })?;
            if end > file_len {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "section {id} out of file bounds: {offset}..{end} > {file_len}"
                )));
            }
        }

        sections.push(SectionEntry::new(id, flags, offset, length, crc32c));
    }

    let mut ranges: Vec<(u64, u64, u16)> = sections
        .iter()
        .filter(|s| s.length > 0)
        .map(|s| (s.offset, s.offset + s.length, s.id))
        .collect();
    ranges.sort_unstable_by_key(|r| r.0);
    for win in ranges.windows(2) {
        let (a_start, a_end, a_id) = win[0];
        let (b_start, _b_end, b_id) = win[1];
        if a_end > b_start {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "section overlap: {a_id} [{a_start}..{a_end}) overlaps {b_id} [{b_start}..)"
            )));
        }
    }

    let mut mandatory: HashMap<u16, SectionEntry> = HashMap::new();
    for section in &sections {
        if MANDATORY_SECTION_IDS.contains(&section.id) {
            if section.offset == 0 || section.length == 0 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "mandatory section {} has zero offset/length",
                    section.id
                )));
            }
            if mandatory.insert(section.id, *section).is_some() {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "mandatory section {} appears more than once",
                    section.id
                )));
            }
        }
    }

    for required in MANDATORY_SECTION_IDS {
        if !mandatory.contains_key(&required) {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "missing mandatory section {required}"
            )));
        }
    }

    for entry in mandatory.values() {
        let raw = section_bytes(bytes, *entry)?;
        let actual = compute_crc32(raw);
        if actual != entry.crc32c {
            return Err(SpatialIoError::CrcMismatch {
                section_id: entry.id,
                expected: entry.crc32c,
                actual,
            });
        }
    }

    let spatial_raw = decode_section_payload(bytes, mandatory[&SECTION_ID_SPATIAL_DOMAIN])?;
    let csr_raw = decode_section_payload(bytes, mandatory[&SECTION_ID_CSR])?;
    let features_raw = decode_section_payload(bytes, mandatory[&SECTION_ID_FEATURE_TABLE])?;
    let meta_raw = decode_section_payload(bytes, mandatory[&SECTION_ID_META_CORE])?;
    let json_raw = decode_section_payload(bytes, mandatory[&SECTION_ID_META_JSON])?;

    let spatial = decode_spatial_domain(&spatial_raw, budget_bytes)?;
    let csr = decode_csr(&csr_raw, budget_bytes)?;
    let features = decode_feature_table(&features_raw, budget_bytes)?;
    let mut meta = decode_metadata_core(&meta_raw)?;
    let (metadata_json, canonical_json_bytes) = decode_metadata_json(&json_raw, budget_bytes)?;

    validate_cross_section_invariants(&spatial, &csr, &features, &meta)?;

    if cfg.validate_hash {
        let computed_hash =
            compute_dataset_hash(&spatial, &csr, &features, &meta, &canonical_json_bytes)?;
        if computed_hash != header_hash {
            return Err(SpatialIoError::HashMismatch);
        }
        meta.dataset_hash = computed_hash;
    } else {
        meta.dataset_hash = header_hash;
    }

    Ok(Dataset::from_parts(
        spatial,
        csr,
        features,
        meta,
        metadata_json,
    ))
}

fn decode_section_payload<'a>(
    bytes: &'a [u8],
    entry: SectionEntry,
) -> Result<std::borrow::Cow<'a, [u8]>, SpatialIoError> {
    let raw = section_bytes(bytes, entry)?;
    if entry.flags & SECTION_FLAG_ZSTD != 0 {
        let decoded = zstd_decompress(raw)?;
        Ok(std::borrow::Cow::Owned(decoded))
    } else {
        Ok(std::borrow::Cow::Borrowed(raw))
    }
}

fn compute_crc32(bytes: &[u8]) -> u32 {
    let mut h = Crc32::new();
    h.update(bytes);
    h.finalize()
}

fn validate_cross_section_invariants(
    spatial: &SpatialDomain,
    csr: &BinsCsr,
    features: &FeatureTable,
    meta: &DatasetMetaCore,
) -> Result<(), SpatialIoError> {
    if spatial.x.len() as u32 != csr.n_bins || spatial.x.len() as u32 != meta.n_bins {
        return Err(SpatialIoError::DimensionMismatch(
            "n_bins mismatch across spatial/csr/metadata".to_string(),
        ));
    }
    if features.rows.len() as u32 != csr.n_genes || features.rows.len() as u32 != meta.n_genes {
        return Err(SpatialIoError::DimensionMismatch(
            "n_genes mismatch across feature/csr/metadata".to_string(),
        ));
    }
    if csr.nnz != meta.nnz {
        return Err(SpatialIoError::DimensionMismatch(
            "nnz mismatch across csr/metadata".to_string(),
        ));
    }
    if spatial.coord_system != meta.coord_system {
        return Err(SpatialIoError::UnsupportedFormat(
            "coord_system mismatch between spatial and metadata".to_string(),
        ));
    }
    if spatial.bin_level != meta.bin_level {
        return Err(SpatialIoError::UnsupportedFormat(
            "bin_level mismatch between spatial and metadata".to_string(),
        ));
    }
    Ok(())
}

fn section_bytes(all: &[u8], entry: SectionEntry) -> Result<&[u8], SpatialIoError> {
    let start = usize::try_from(entry.offset)
        .map_err(|_| SpatialIoError::UnsupportedFormat("section offset too large".to_string()))?;
    let end_u64 = entry
        .offset
        .checked_add(entry.length)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("section end overflow".to_string()))?;
    let end = usize::try_from(end_u64)
        .map_err(|_| SpatialIoError::UnsupportedFormat("section end too large".to_string()))?;

    all.get(start..end)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("section slice out of bounds".to_string()))
}

fn ensure_alloc_within_budget(
    bytes_needed: u64,
    budget: u64,
    ctx: &str,
) -> Result<(), SpatialIoError> {
    if bytes_needed > budget {
        return Err(SpatialIoError::MemoryLimitExceeded(format!(
            "{ctx}: would allocate {bytes_needed} bytes > budget {budget}"
        )));
    }
    Ok(())
}

fn decode_spatial_domain(data: &[u8], budget: u64) -> Result<SpatialDomain, SpatialIoError> {
    let mut r = SliceCursor::new(data, "spatial");

    let n_bins = r.read_u32("n_bins")? as usize;
    let coord_system = u8_to_coord_system(r.read_u8("coord_system")?)?;
    let bin_level = r.read_u8("bin_level")?;
    let flags = r.read_u16("flags")?;
    let has_grid = (flags & 1) != 0;

    let needed = (n_bins as u64)
        .checked_mul((std::mem::size_of::<f32>() * 2 + std::mem::size_of::<u32>()) as u64)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat("spatial n_bins size overflow".to_string())
        })?;
    ensure_alloc_within_budget(needed, budget, "spatial domain declared size")?;

    let x = r.read_f32_vec(n_bins, "x")?;
    let y = r.read_f32_vec(n_bins, "y")?;

    let grid_row = if has_grid {
        Some(r.read_u32_vec(n_bins, "grid_row")?)
    } else {
        None
    };
    let grid_col = if has_grid {
        Some(r.read_u32_vec(n_bins, "grid_col")?)
    } else {
        None
    };

    let bin_id = r.read_u32_vec(n_bins, "bin_id")?;

    let n_bits = r.read_u64("tissue_mask.n_bits")?;
    let raw_bytes_len = r.read_u64("tissue_mask.raw_bytes_len")? as usize;
    if !raw_bytes_len.is_multiple_of(8) {
        return Err(SpatialIoError::UnsupportedFormat(
            "spatial tissue raw length is not 8-byte aligned".to_string(),
        ));
    }
    ensure_alloc_within_budget(raw_bytes_len as u64, budget, "tissue mask raw bytes")?;
    let raw_bytes = r.read_bytes(raw_bytes_len, "tissue_mask.raw")?;

    if n_bits as usize != n_bins {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "spatial tissue bits {n_bits} != n_bins {n_bins}"
        )));
    }

    let mut words = Vec::with_capacity(raw_bytes_len / 8);
    for chunk in raw_bytes.chunks_exact(8) {
        let mut arr = [0_u8; 8];
        arr.copy_from_slice(chunk);
        words.push(u64::from_le_bytes(arr));
    }

    let tissue_mask = tissue_mask_from_u64_words(&words, n_bins);

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "spatial section has trailing bytes".to_string(),
        ));
    }

    for (i, &id) in bin_id.iter().enumerate() {
        if id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "spatial domain invariant: bin_id[{i}] = {id}, expected {i}"
            )));
        }
    }

    if let (Some(rows), Some(cols)) = (&grid_row, &grid_col) {
        for i in 1..n_bins {
            let prev = (rows[i - 1], cols[i - 1]);
            let curr = (rows[i], cols[i]);
            if curr < prev {
                return Err(SpatialIoError::InvalidCsr(
                    "spatial domain invariant: grid coordinates are not sorted".to_string(),
                ));
            }
        }
    } else {
        for i in 1..n_bins {
            let by_y = total_cmp_f32(y[i], y[i - 1]);
            if by_y.is_lt() {
                return Err(SpatialIoError::InvalidCsr(
                    "spatial domain invariant: y is not sorted".to_string(),
                ));
            }
            if by_y.is_eq() && total_cmp_f32(x[i], x[i - 1]).is_lt() {
                return Err(SpatialIoError::InvalidCsr(
                    "spatial domain invariant: x is not sorted within equal y".to_string(),
                ));
            }
        }
    }

    SpatialDomain::new(
        x,
        y,
        grid_row,
        grid_col,
        bin_id,
        tissue_mask,
        coord_system,
        bin_level,
    )
}

fn decode_csr(data: &[u8], budget: u64) -> Result<BinsCsr, SpatialIoError> {
    let mut r = SliceCursor::new(data, "csr");

    let n_bins = r.read_u32("n_bins")?;
    let n_genes = r.read_u32("n_genes")?;
    let nnz = r.read_u64("nnz")?;
    let normalized = r.read_u8("normalized")? != 0;
    let flags = r.read_u8("flags")?;
    r.read_bytes(6, "reserved")?;

    let indptr_u32 = flags & CSR_FLAG_INDPTR_U32 != 0;

    let elem_bytes = 4_u64
        + 4
        + if indptr_u32 { 4 } else { 8 };
    let total = (nnz)
        .checked_mul(elem_bytes)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("csr nnz size overflow".to_string()))?;
    ensure_alloc_within_budget(total, budget, "csr declared size")?;

    let indptr: Indptr = if indptr_u32 {
        Indptr::U32(r.read_u32_vec(n_bins as usize + 1, "indptr")?)
    } else {
        Indptr::U64(r.read_u64_vec(n_bins as usize + 1, "indptr")?)
    };
    let indices = r.read_u32_vec(nnz as usize, "indices")?;
    let data_values = r.read_f32_vec(nnz as usize, "data")?;

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "csr section has trailing bytes".to_string(),
        ));
    }

    if indptr.first().unwrap_or(1) != 0 {
        return Err(SpatialIoError::InvalidCsr(
            "csr invariant: indptr[0] must be 0".to_string(),
        ));
    }
    if indptr.last().unwrap_or(0) != nnz {
        return Err(SpatialIoError::InvalidCsr(format!(
            "csr invariant: indptr[last] {} != nnz {nnz}",
            indptr.last().unwrap_or(0)
        )));
    }

    for i in 1..indptr.len() {
        if indptr.get(i) < indptr.get(i - 1) {
            return Err(SpatialIoError::InvalidCsr(format!(
                "csr invariant: indptr not monotonic at {i}"
            )));
        }
    }

    for row in 0..n_bins as usize {
        let start = indptr.get(row) as usize;
        let end = indptr.get(row + 1) as usize;

        let mut prev: Option<u32> = None;
        for idx in start..end {
            let gene = indices[idx];
            if gene >= n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "csr invariant: gene index out of bounds in row {row}: {gene} >= {n_genes}"
                )));
            }
            if let Some(p) = prev
                && gene <= p
            {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "csr invariant: row {row} indices are not strictly increasing"
                )));
            }
            prev = Some(gene);
            if !normalized {
                ensure_f32_finite_nonneg(data_values[idx])?;
            } else if !data_values[idx].is_finite() {
                return Err(SpatialIoError::InvalidFloat(
                    "csr data contains non-finite value".to_string(),
                ));
            }
        }
    }

    Ok(BinsCsr {
        indptr,
        indices,
        data: data_values,
        n_bins,
        n_genes,
        nnz,
        normalized,
    })
}

fn decode_feature_table(data: &[u8], budget: u64) -> Result<FeatureTable, SpatialIoError> {
    let mut r = SliceCursor::new(data, "feature_table");

    let n_genes = r.read_u32("n_genes")? as usize;
    ensure_alloc_within_budget(
        (n_genes as u64).saturating_mul(64),
        budget,
        "feature table declared size",
    )?;
    let mut rows = Vec::with_capacity(n_genes);
    let mut seen_names = HashSet::with_capacity(n_genes);
    let mut prev_name: Option<String> = None;

    for i in 0..n_genes {
        let gene_id = r.read_u32("gene_id")?;
        if gene_id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table invariant: gene_id {gene_id} != {i}"
            )));
        }

        let feature_id = r.read_len_prefixed_string("feature_id")?;
        let gene_name = r.read_len_prefixed_string("gene_name")?;
        let feature_type = r.read_len_prefixed_string("feature_type")?;

        if !seen_names.insert(gene_name.clone()) {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table invariant: duplicate gene_name {gene_name}"
            )));
        }
        if let Some(prev) = &prev_name
            && gene_name <= *prev
        {
            return Err(SpatialIoError::InvalidCsr(
                "feature table invariant: gene_name not strictly sorted".to_string(),
            ));
        }
        prev_name = Some(gene_name.clone());

        rows.push(FeatureRow {
            gene_id,
            feature_id,
            gene_name,
            feature_type,
        });
    }

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "feature table section has trailing bytes".to_string(),
        ));
    }

    Ok(FeatureTable { rows })
}

fn decode_metadata_core(data: &[u8]) -> Result<DatasetMetaCore, SpatialIoError> {
    let mut r = SliceCursor::new(data, "metadata_core");

    let dataset_name = r.read_len_prefixed_string("dataset_name")?;
    let source_format = r.read_len_prefixed_string("source_format")?;
    let bin_level = r.read_u8("bin_level")?;
    let n_bins = r.read_u32("n_bins")?;
    let n_genes = r.read_u32("n_genes")?;
    let nnz = r.read_u64("nnz")?;
    let coord_system = u8_to_coord_system(r.read_u8("coord_system")?)?;
    let normalized = r.read_u8("normalized")? != 0;

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "metadata_core section has trailing bytes".to_string(),
        ));
    }

    Ok(DatasetMetaCore {
        dataset_name,
        source_format,
        bin_level,
        n_bins,
        n_genes,
        nnz,
        coord_system,
        normalized,
        dataset_hash: [0_u8; 16],
    })
}

fn decode_metadata_json(data: &[u8], budget: u64) -> Result<(Value, Vec<u8>), SpatialIoError> {
    let mut r = SliceCursor::new(data, "metadata_json");

    let json_len = r.read_u64("json_len")?;
    ensure_alloc_within_budget(json_len, budget, "metadata json declared length")?;
    let json_bytes = r.read_bytes(json_len as usize, "json")?.to_vec();
    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "metadata_json section has trailing bytes".to_string(),
        ));
    }

    let parsed: Value = serde_json::from_slice(&json_bytes).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!("metadata json parse error: {e}"))
    })?;

    let canonical = canonicalize_json(&parsed);
    let mut canonical_bytes = Vec::new();
    write_canonical_json(&mut canonical_bytes, &canonical)?;

    if canonical_bytes != json_bytes {
        return Err(SpatialIoError::UnsupportedFormat(
            "metadata json not canonical".to_string(),
        ));
    }

    Ok((canonical, canonical_bytes))
}

fn u16_at(data: &[u8], offset: usize) -> Result<u16, SpatialIoError> {
    let slice = data
        .get(offset..offset + 2)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("read out of bounds".to_string()))?;
    let mut arr = [0_u8; 2];
    arr.copy_from_slice(slice);
    Ok(u16::from_le_bytes(arr))
}

fn u32_at(data: &[u8], offset: usize) -> Result<u32, SpatialIoError> {
    let slice = data
        .get(offset..offset + 4)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("read out of bounds".to_string()))?;
    let mut arr = [0_u8; 4];
    arr.copy_from_slice(slice);
    Ok(u32::from_le_bytes(arr))
}

fn u64_at(data: &[u8], offset: usize) -> Result<u64, SpatialIoError> {
    let slice = data
        .get(offset..offset + 8)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("read out of bounds".to_string()))?;
    let mut arr = [0_u8; 8];
    arr.copy_from_slice(slice);
    Ok(u64::from_le_bytes(arr))
}

fn u8_to_coord_system(v: u8) -> Result<CoordSystem, SpatialIoError> {
    match v {
        0 => Ok(CoordSystem::Grid),
        1 => Ok(CoordSystem::Pixel),
        2 => Ok(CoordSystem::Micron),
        _ => Err(SpatialIoError::UnsupportedFormat(format!(
            "invalid coord_system value: {v}"
        ))),
    }
}

struct SliceCursor<'a> {
    data: &'a [u8],
    pos: usize,
    section: &'static str,
}

impl<'a> SliceCursor<'a> {
    fn new(data: &'a [u8], section: &'static str) -> Self {
        Self {
            data,
            pos: 0,
            section,
        }
    }

    fn is_eof(&self) -> bool {
        self.pos == self.data.len()
    }

    fn read_bytes(&mut self, len: usize, field: &str) -> Result<&'a [u8], SpatialIoError> {
        let end = self.pos.checked_add(len).ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!(
                "{}.{} read overflow",
                self.section, field
            ))
        })?;
        let slice = self.data.get(self.pos..end).ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!(
                "{}.{} out of bounds at {}..{}",
                self.section, field, self.pos, end
            ))
        })?;
        self.pos = end;
        Ok(slice)
    }

    fn read_u8(&mut self, field: &str) -> Result<u8, SpatialIoError> {
        Ok(self.read_bytes(1, field)?[0])
    }

    fn read_u16(&mut self, field: &str) -> Result<u16, SpatialIoError> {
        let mut arr = [0_u8; 2];
        arr.copy_from_slice(self.read_bytes(2, field)?);
        Ok(u16::from_le_bytes(arr))
    }

    fn read_u32(&mut self, field: &str) -> Result<u32, SpatialIoError> {
        let mut arr = [0_u8; 4];
        arr.copy_from_slice(self.read_bytes(4, field)?);
        Ok(u32::from_le_bytes(arr))
    }

    fn read_u64(&mut self, field: &str) -> Result<u64, SpatialIoError> {
        let mut arr = [0_u8; 8];
        arr.copy_from_slice(self.read_bytes(8, field)?);
        Ok(u64::from_le_bytes(arr))
    }

    fn read_u32_vec(&mut self, n: usize, field: &str) -> Result<Vec<u32>, SpatialIoError> {
        let bytes = self.read_bytes(n * 4, field)?;
        let mut out: Vec<u32> = vec![0; n];
        bytemuck::cast_slice_mut::<u32, u8>(&mut out).copy_from_slice(bytes);
        Ok(out)
    }

    fn read_u64_vec(&mut self, n: usize, field: &str) -> Result<Vec<u64>, SpatialIoError> {
        let bytes = self.read_bytes(n * 8, field)?;
        let mut out: Vec<u64> = vec![0; n];
        bytemuck::cast_slice_mut::<u64, u8>(&mut out).copy_from_slice(bytes);
        Ok(out)
    }

    fn read_f32_vec(&mut self, n: usize, field: &str) -> Result<Vec<f32>, SpatialIoError> {
        let bytes = self.read_bytes(n * 4, field)?;
        let mut out: Vec<f32> = vec![0.0; n];
        bytemuck::cast_slice_mut::<f32, u8>(&mut out).copy_from_slice(bytes);
        Ok(out)
    }

    fn read_len_prefixed_string(&mut self, field: &str) -> Result<String, SpatialIoError> {
        let len = self.read_u32(&format!("{field}.len"))? as usize;
        let bytes = self.read_bytes(len, field)?;
        String::from_utf8(bytes.to_vec()).map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!(
                "{}.{} invalid utf8: {}",
                self.section, field, e
            ))
        })
    }
}
