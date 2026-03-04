use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;

use bitvec::vec::BitVec;
use serde_json::Value;

use crate::api::dataset::Dataset;
use crate::binary::format::{
    ENDIAN_LITTLE, HEADER_SIZE, KIRA_SPATIAL_BIN_VERSION, MAGIC, SECTION_ENTRY_SIZE,
    SECTION_ID_CSR, SECTION_ID_FEATURE_TABLE, SECTION_ID_META_CORE, SECTION_ID_META_JSON,
    SECTION_ID_SPATIAL_DOMAIN, SectionEntry,
};
use crate::binary::hash::compute_dataset_hash;
use crate::determinism::float::{ensure_f32_finite_nonneg, total_cmp_f32};
use crate::determinism::json::{canonicalize_json, write_canonical_json};
use crate::error::SpatialIoError;
use crate::model::{
    coord::CoordSystem,
    csr::BinsCsr,
    features::{FeatureRow, FeatureTable},
    metadata::DatasetMetaCore,
    spatial_domain::SpatialDomain,
};

/// Reads and validates a deterministic `.kira-spatial.bin` dataset.
pub fn read_kira_bin<P: AsRef<Path>>(p: P) -> Result<Dataset, SpatialIoError> {
    let mut file = File::open(p.as_ref())?;
    let mut bytes_vec = Vec::new();
    file.read_to_end(&mut bytes_vec)?;
    let bytes: &[u8] = &bytes_vec;
    if bytes.len() < HEADER_SIZE as usize {
        return Err(SpatialIoError::UnsupportedFormat(
            "file too small for header".to_string(),
        ));
    }

    if &bytes[0..8] != MAGIC.as_slice() {
        return Err(SpatialIoError::UnsupportedFormat(
            "invalid magic: expected KIRASPAT".to_string(),
        ));
    }

    let version = read_u16_at(bytes, 8, "header.version")?;
    if version != KIRA_SPATIAL_BIN_VERSION {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "unsupported version: {}",
            version
        )));
    }

    let endian = *bytes
        .get(10)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("missing endian byte".to_string()))?;
    if endian != ENDIAN_LITTLE {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "unsupported endian flag: {}",
            endian
        )));
    }

    let section_count = read_u16_at(bytes, 11, "header.section_count")? as usize;
    if section_count < 5 {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "section_count too small: {}",
            section_count
        )));
    }

    let mut header_hash = [0_u8; 16];
    header_hash.copy_from_slice(
        bytes
            .get(13..29)
            .ok_or_else(|| SpatialIoError::UnsupportedFormat("missing header hash".to_string()))?,
    );

    let table_start = HEADER_SIZE as usize;
    let table_len = section_count
        .checked_mul(SECTION_ENTRY_SIZE as usize)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat("section table size overflow".to_string())
        })?;
    let table_end = table_start.checked_add(table_len).ok_or_else(|| {
        SpatialIoError::UnsupportedFormat("section table end overflow".to_string())
    })?;
    if table_end > bytes.len() {
        return Err(SpatialIoError::UnsupportedFormat(
            "section table out of file bounds".to_string(),
        ));
    }

    let mut sections = Vec::with_capacity(section_count);
    for i in 0..section_count {
        let base = table_start + i * SECTION_ENTRY_SIZE as usize;
        let id = read_u16_at(bytes, base, "section.id")?;
        let offset = read_u64_at(bytes, base + 2, "section.offset")?;
        let length = read_u64_at(bytes, base + 10, "section.length")?;

        if offset % 64 != 0 {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "section {} offset is not 64-byte aligned: {}",
                id, offset
            )));
        }

        let end = offset.checked_add(length).ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!("section {} end overflow", id))
        })?;
        if end > bytes.len() as u64 {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "section {} out of file bounds: {}..{} > {}",
                id,
                offset,
                end,
                bytes.len()
            )));
        }

        sections.push(SectionEntry::new(id, offset, length));
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
                "section overlap: {} [{}..{}) overlaps {} [{}..)",
                a_id, a_start, a_end, b_id, b_start
            )));
        }
    }

    let mut mandatory: HashMap<u16, SectionEntry> = HashMap::new();
    for section in &sections {
        if matches!(
            section.id,
            SECTION_ID_SPATIAL_DOMAIN
                | SECTION_ID_CSR
                | SECTION_ID_FEATURE_TABLE
                | SECTION_ID_META_CORE
                | SECTION_ID_META_JSON
        ) {
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

    for required in [
        SECTION_ID_SPATIAL_DOMAIN,
        SECTION_ID_CSR,
        SECTION_ID_FEATURE_TABLE,
        SECTION_ID_META_CORE,
        SECTION_ID_META_JSON,
    ] {
        if !mandatory.contains_key(&required) {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "missing mandatory section {}",
                required
            )));
        }
    }

    let spatial =
        decode_spatial_domain(section_bytes(bytes, mandatory[&SECTION_ID_SPATIAL_DOMAIN])?)?;
    let csr = decode_csr(section_bytes(bytes, mandatory[&SECTION_ID_CSR])?)?;
    let features =
        decode_feature_table(section_bytes(bytes, mandatory[&SECTION_ID_FEATURE_TABLE])?)?;
    let meta = decode_metadata_core(section_bytes(bytes, mandatory[&SECTION_ID_META_CORE])?)?;
    let (metadata_json, canonical_json_bytes) =
        decode_metadata_json(section_bytes(bytes, mandatory[&SECTION_ID_META_JSON])?)?;

    validate_cross_section_invariants(&spatial, &csr, &features, &meta)?;

    let computed_hash =
        compute_dataset_hash(&spatial, &csr, &features, &meta, &canonical_json_bytes)?;
    if computed_hash != header_hash || computed_hash != meta.dataset_hash {
        return Err(SpatialIoError::UnsupportedFormat(
            "dataset hash mismatch".to_string(),
        ));
    }

    Ok(Dataset::from_parts(
        spatial,
        csr,
        features,
        meta,
        metadata_json,
    ))
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

fn decode_spatial_domain(data: &[u8]) -> Result<SpatialDomain, SpatialIoError> {
    let mut r = SliceCursor::new(data, "spatial");

    let n_bins = r.read_u32("n_bins")? as usize;
    let coord_system = u8_to_coord_system(r.read_u8("coord_system")?)?;
    let bin_level = r.read_u8("bin_level")?;
    let flags = r.read_u16("flags")?;
    let has_grid = (flags & 1) != 0;

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
    let raw_bytes = r.read_bytes(raw_bytes_len, "tissue_mask.raw")?;

    if n_bits as usize != n_bins {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "spatial tissue bits {} != n_bins {}",
            n_bits, n_bins
        )));
    }

    if !raw_bytes_len.is_multiple_of(std::mem::size_of::<usize>()) {
        return Err(SpatialIoError::UnsupportedFormat(
            "spatial tissue raw length is not word-aligned".to_string(),
        ));
    }

    let mut words = Vec::with_capacity(raw_bytes_len / std::mem::size_of::<usize>());
    for chunk in raw_bytes.chunks_exact(std::mem::size_of::<usize>()) {
        #[cfg(target_pointer_width = "64")]
        {
            let mut arr = [0_u8; 8];
            arr.copy_from_slice(chunk);
            words.push(usize::from_le_bytes(arr));
        }
        #[cfg(target_pointer_width = "32")]
        {
            let mut arr = [0_u8; 4];
            arr.copy_from_slice(chunk);
            words.push(usize::from_le_bytes(arr));
        }
    }

    let mut tissue_mask = BitVec::from_vec(words);
    if tissue_mask.len() < n_bits as usize {
        return Err(SpatialIoError::UnsupportedFormat(
            "spatial tissue mask bit length shorter than declared".to_string(),
        ));
    }
    tissue_mask.truncate(n_bits as usize);

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "spatial section has trailing bytes".to_string(),
        ));
    }

    for (i, &id) in bin_id.iter().enumerate() {
        if id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "spatial domain invariant: bin_id[{}] = {}, expected {}",
                i, id, i
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

fn decode_csr(data: &[u8]) -> Result<BinsCsr, SpatialIoError> {
    let mut r = SliceCursor::new(data, "csr");

    let n_bins = r.read_u32("n_bins")?;
    let n_genes = r.read_u32("n_genes")?;
    let nnz = r.read_u64("nnz")?;
    let normalized = r.read_u8("normalized")? != 0;
    r.read_bytes(7, "reserved")?;

    let indptr = r.read_u64_vec(n_bins as usize + 1, "indptr")?;
    let indices = r.read_u32_vec(nnz as usize, "indices")?;
    let data_values = r.read_f32_vec(nnz as usize, "data")?;

    if !r.is_eof() {
        return Err(SpatialIoError::UnsupportedFormat(
            "csr section has trailing bytes".to_string(),
        ));
    }

    if indptr.first().copied().unwrap_or(1) != 0 {
        return Err(SpatialIoError::InvalidCsr(
            "csr invariant: indptr[0] must be 0".to_string(),
        ));
    }
    if indptr.last().copied().unwrap_or(0) != nnz {
        return Err(SpatialIoError::InvalidCsr(format!(
            "csr invariant: indptr[last] {} != nnz {}",
            indptr.last().copied().unwrap_or(0),
            nnz
        )));
    }

    for i in 1..indptr.len() {
        if indptr[i] < indptr[i - 1] {
            return Err(SpatialIoError::InvalidCsr(format!(
                "csr invariant: indptr not monotonic at {}",
                i
            )));
        }
    }

    for row in 0..n_bins as usize {
        let start = indptr[row] as usize;
        let end = indptr[row + 1] as usize;

        let mut prev: Option<u32> = None;
        for idx in start..end {
            let gene = indices[idx];
            if gene >= n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "csr invariant: gene index out of bounds in row {}: {} >= {}",
                    row, gene, n_genes
                )));
            }
            if let Some(p) = prev
                && gene <= p
            {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "csr invariant: row {} indices are not strictly increasing",
                    row
                )));
            }
            prev = Some(gene);
            ensure_f32_finite_nonneg(data_values[idx])?;
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

fn decode_feature_table(data: &[u8]) -> Result<FeatureTable, SpatialIoError> {
    let mut r = SliceCursor::new(data, "feature_table");

    let n_genes = r.read_u32("n_genes")? as usize;
    let mut rows = Vec::with_capacity(n_genes);
    let mut seen_names = HashSet::with_capacity(n_genes);
    let mut prev_name: Option<String> = None;

    for i in 0..n_genes {
        let gene_id = r.read_u32("gene_id")?;
        if gene_id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table invariant: gene_id {} != {}",
                gene_id, i
            )));
        }

        let gene_name = r.read_len_prefixed_string("gene_name")?;
        let feature_type = r.read_len_prefixed_string("feature_type")?;

        if !seen_names.insert(gene_name.clone()) {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table invariant: duplicate gene_name {}",
                gene_name
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
    let hash = r.read_bytes(16, "dataset_hash")?;

    let mut dataset_hash = [0_u8; 16];
    dataset_hash.copy_from_slice(hash);

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
        dataset_hash,
    })
}

fn decode_metadata_json(data: &[u8]) -> Result<(Value, Vec<u8>), SpatialIoError> {
    let mut r = SliceCursor::new(data, "metadata_json");

    let json_len = r.read_u64("json_len")? as usize;
    let json_bytes = r.read_bytes(json_len, "json")?.to_vec();
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

fn read_u16_at(data: &[u8], offset: usize, ctx: &str) -> Result<u16, SpatialIoError> {
    let slice = data
        .get(offset..offset + 2)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat(format!("{} out of bounds", ctx)))?;
    let mut arr = [0_u8; 2];
    arr.copy_from_slice(slice);
    Ok(u16::from_le_bytes(arr))
}

fn read_u64_at(data: &[u8], offset: usize, ctx: &str) -> Result<u64, SpatialIoError> {
    let slice = data
        .get(offset..offset + 8)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat(format!("{} out of bounds", ctx)))?;
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
            "invalid coord_system value: {}",
            v
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
            SpatialIoError::UnsupportedFormat(format!("{}.{field} read overflow", self.section))
        })?;
        let slice = self.data.get(self.pos..end).ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!(
                "{}.{field} out of bounds at {}..{}",
                self.section, self.pos, end
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

    fn read_f32(&mut self, field: &str) -> Result<f32, SpatialIoError> {
        let mut arr = [0_u8; 4];
        arr.copy_from_slice(self.read_bytes(4, field)?);
        Ok(f32::from_le_bytes(arr))
    }

    fn read_u32_vec(&mut self, n: usize, field: &str) -> Result<Vec<u32>, SpatialIoError> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.read_u32(field)?);
        }
        Ok(out)
    }

    fn read_u64_vec(&mut self, n: usize, field: &str) -> Result<Vec<u64>, SpatialIoError> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.read_u64(field)?);
        }
        Ok(out)
    }

    fn read_f32_vec(&mut self, n: usize, field: &str) -> Result<Vec<f32>, SpatialIoError> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.read_f32(field)?);
        }
        Ok(out)
    }

    fn read_len_prefixed_string(&mut self, field: &str) -> Result<String, SpatialIoError> {
        let len = self.read_u32(&format!("{}.len", field))? as usize;
        let bytes = self.read_bytes(len, field)?;
        String::from_utf8(bytes.to_vec()).map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!(
                "{}.{} invalid utf8: {}",
                self.section, field, e
            ))
        })
    }
}
