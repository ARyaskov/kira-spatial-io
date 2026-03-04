use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::mem::size_of;
use std::path::Path;

use crate::api::dataset::Dataset;
use crate::binary::format::{
    ENDIAN_LITTLE, HEADER_SIZE, Header, MAGIC, SECTION_COUNT, SECTION_ENTRY_SIZE, SECTION_ID_CSR,
    SECTION_ID_FEATURE_TABLE, SECTION_ID_META_CORE, SECTION_ID_META_JSON,
    SECTION_ID_SPATIAL_DOMAIN, SectionEntry,
};
use crate::binary::hash::compute_dataset_hash;
use crate::determinism::json::{canonicalize_json, write_canonical_json};
use crate::error::SpatialIoError;
use crate::model::{
    coord::CoordSystem, csr::BinsCsr, metadata::DatasetMetaCore, spatial_domain::SpatialDomain,
};

/// Writes a dataset to deterministic `.kira-spatial.bin`.
pub fn write_kira_bin<P: AsRef<Path>>(p: P, ds: &Dataset) -> Result<(), SpatialIoError> {
    let canonical_json = canonicalize_json(ds.metadata_json());
    let mut canonical_json_bytes = Vec::new();
    write_canonical_json(&mut canonical_json_bytes, &canonical_json)?;

    let dataset_hash = compute_dataset_hash(
        ds.spatial_domain(),
        ds.expression_csr(),
        ds.features(),
        ds.metadata_core(),
        &canonical_json_bytes,
    )?;

    let file = File::create(p)?;
    let mut writer = BufWriter::new(file);

    let placeholder_entries = [
        SectionEntry::new(SECTION_ID_SPATIAL_DOMAIN, 0, 0),
        SectionEntry::new(SECTION_ID_CSR, 0, 0),
        SectionEntry::new(SECTION_ID_FEATURE_TABLE, 0, 0),
        SectionEntry::new(SECTION_ID_META_CORE, 0, 0),
        SectionEntry::new(SECTION_ID_META_JSON, 0, 0),
    ];

    write_header(&mut writer, &Header::new(dataset_hash))?;
    write_section_table(&mut writer, &placeholder_entries)?;

    let mut entries = Vec::with_capacity(SECTION_COUNT as usize);

    entries.push(write_section(
        &mut writer,
        SECTION_ID_SPATIAL_DOMAIN,
        |w| write_spatial_domain_section(w, ds.spatial_domain()),
    )?);

    entries.push(write_section(&mut writer, SECTION_ID_CSR, |w| {
        write_csr_section(w, ds.expression_csr())
    })?);

    entries.push(write_section(&mut writer, SECTION_ID_FEATURE_TABLE, |w| {
        write_feature_table_section(w, ds.features())
    })?);

    entries.push(write_section(&mut writer, SECTION_ID_META_CORE, |w| {
        write_metadata_core_section(w, ds.metadata_core(), dataset_hash)
    })?);

    entries.push(write_section(&mut writer, SECTION_ID_META_JSON, |w| {
        write_metadata_json_section(w, &canonical_json_bytes)
    })?);

    writer.flush()?;
    writer.seek(SeekFrom::Start(0))?;
    write_header(&mut writer, &Header::new(dataset_hash))?;
    write_section_table(&mut writer, &entries)?;
    writer.flush()?;
    writer.get_ref().sync_all()?;

    Ok(())
}

/// Pads writer position to next 64-byte boundary.
pub fn pad_to_64<W: Write + Seek>(w: &mut W) -> Result<(), SpatialIoError> {
    let pos = w.stream_position()?;
    let rem = (pos % 64) as usize;
    if rem == 0 {
        return Ok(());
    }
    let need = 64 - rem;
    let zeros = [0_u8; 64];
    w.write_all(&zeros[..need])?;
    Ok(())
}

fn write_section<W, F>(
    w: &mut W,
    id: u16,
    mut write_payload: F,
) -> Result<SectionEntry, SpatialIoError>
where
    W: Write + Seek,
    F: FnMut(&mut W) -> Result<(), SpatialIoError>,
{
    pad_to_64(w)?;
    let offset = w.stream_position()?;
    write_payload(w)?;
    let end = w.stream_position()?;
    let length = end
        .checked_sub(offset)
        .ok_or_else(|| SpatialIoError::InvalidCsr("section length underflow".to_string()))?;
    Ok(SectionEntry::new(id, offset, length))
}

fn write_header<W: Write>(w: &mut W, header: &Header) -> Result<(), SpatialIoError> {
    w.write_all(&MAGIC)?;
    write_u16_le(w, header.version)?;
    w.write_all(&[header.endian])?;
    write_u16_le(w, header.section_count)?;
    w.write_all(&header.dataset_hash)?;

    let written = 8_u64 + 2 + 1 + 2 + 16;
    let pad = HEADER_SIZE
        .checked_sub(written)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("invalid header size".to_string()))?;
    write_zeros(w, pad as usize)?;
    Ok(())
}

fn write_section_table<W: Write>(
    w: &mut W,
    entries: &[SectionEntry],
) -> Result<(), SpatialIoError> {
    if entries.len() != SECTION_COUNT as usize {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "invalid section table size: {}",
            entries.len()
        )));
    }

    for entry in entries {
        write_u16_le(w, entry.id)?;
        write_u64_le(w, entry.offset)?;
        write_u64_le(w, entry.length)?;
    }

    let expected = HEADER_SIZE
        .checked_add(SECTION_ENTRY_SIZE * SECTION_COUNT as u64)
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("table size overflow".to_string()))?;
    let current = HEADER_SIZE + entries.len() as u64 * SECTION_ENTRY_SIZE;
    if current != expected {
        return Err(SpatialIoError::UnsupportedFormat(
            "unexpected section table size".to_string(),
        ));
    }

    Ok(())
}

fn write_spatial_domain_section<W: Write>(
    w: &mut W,
    domain: &SpatialDomain,
) -> Result<(), SpatialIoError> {
    let n_bins = domain.x.len();
    if domain.y.len() != n_bins
        || domain.bin_id.len() != n_bins
        || domain.tissue_mask.len() != n_bins
    {
        return Err(SpatialIoError::DimensionMismatch(
            "spatial domain arrays have inconsistent lengths".to_string(),
        ));
    }

    let has_grid = match (&domain.grid_row, &domain.grid_col) {
        (Some(rows), Some(cols)) => {
            if rows.len() != n_bins || cols.len() != n_bins {
                return Err(SpatialIoError::DimensionMismatch(
                    "grid arrays length mismatch".to_string(),
                ));
            }
            true
        }
        (None, None) => false,
        _ => {
            return Err(SpatialIoError::DimensionMismatch(
                "grid_row and grid_col must both be present or absent".to_string(),
            ));
        }
    };

    write_u32_le(
        w,
        u32::try_from(n_bins).map_err(|_| {
            SpatialIoError::DimensionMismatch("n_bins does not fit u32".to_string())
        })?,
    )?;
    w.write_all(&[coord_system_to_u8(domain.coord_system)])?;
    w.write_all(&[domain.bin_level])?;

    let mut flags = 0_u16;
    if has_grid {
        flags |= 1;
    }
    write_u16_le(w, flags)?;

    write_f32_slice_le(w, &domain.x)?;
    write_f32_slice_le(w, &domain.y)?;

    if let (Some(rows), Some(cols)) = (&domain.grid_row, &domain.grid_col) {
        write_u32_slice_le(w, rows)?;
        write_u32_slice_le(w, cols)?;
    }

    write_u32_slice_le(w, &domain.bin_id)?;

    let n_bits = domain.tissue_mask.len() as u64;
    write_u64_le(w, n_bits)?;

    let raw = domain.tissue_mask.as_raw_slice();
    let raw_bytes_len = (raw.len() as u64)
        .checked_mul(size_of::<usize>() as u64)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat("tissue_mask raw size overflow".to_string())
        })?;
    write_u64_le(w, raw_bytes_len)?;
    for word in raw {
        w.write_all(&word.to_le_bytes())?;
    }

    Ok(())
}

fn write_csr_section<W: Write>(w: &mut W, csr: &BinsCsr) -> Result<(), SpatialIoError> {
    if csr.indptr.len() != csr.n_bins as usize + 1 {
        return Err(SpatialIoError::InvalidCsr(
            "indptr length does not match n_bins + 1".to_string(),
        ));
    }
    if csr.indptr.last().copied().unwrap_or(1) != csr.nnz {
        return Err(SpatialIoError::InvalidCsr(
            "indptr[last] does not match nnz".to_string(),
        ));
    }
    if csr.indices.len() != csr.nnz as usize || csr.data.len() != csr.nnz as usize {
        return Err(SpatialIoError::InvalidCsr(
            "indices/data length does not match nnz".to_string(),
        ));
    }

    write_u32_le(w, csr.n_bins)?;
    write_u32_le(w, csr.n_genes)?;
    write_u64_le(w, csr.nnz)?;
    w.write_all(&[u8::from(csr.normalized)])?;
    write_zeros(w, 7)?;

    write_u64_slice_le(w, &csr.indptr)?;
    write_u32_slice_le(w, &csr.indices)?;
    write_f32_slice_le(w, &csr.data)?;
    Ok(())
}

fn write_feature_table_section<W: Write>(
    w: &mut W,
    table: &crate::model::features::FeatureTable,
) -> Result<(), SpatialIoError> {
    write_u32_le(
        w,
        u32::try_from(table.rows.len()).map_err(|_| {
            SpatialIoError::DimensionMismatch("n_genes does not fit u32".to_string())
        })?,
    )?;

    for (i, row) in table.rows.iter().enumerate() {
        if row.gene_id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table gene_id is not canonical at row {i}"
            )));
        }
        write_u32_le(w, row.gene_id)?;
        write_len_prefixed_str(w, &row.gene_name)?;
        write_len_prefixed_str(w, &row.feature_type)?;
    }

    Ok(())
}

fn write_metadata_core_section<W: Write>(
    w: &mut W,
    meta: &DatasetMetaCore,
    dataset_hash: [u8; 16],
) -> Result<(), SpatialIoError> {
    write_len_prefixed_str(w, &meta.dataset_name)?;
    write_len_prefixed_str(w, &meta.source_format)?;
    w.write_all(&[meta.bin_level])?;
    write_u32_le(w, meta.n_bins)?;
    write_u32_le(w, meta.n_genes)?;
    write_u64_le(w, meta.nnz)?;
    w.write_all(&[coord_system_to_u8(meta.coord_system)])?;
    w.write_all(&[u8::from(meta.normalized)])?;
    w.write_all(&dataset_hash)?;
    Ok(())
}

fn write_metadata_json_section<W: Write>(
    w: &mut W,
    canonical_json_bytes: &[u8],
) -> Result<(), SpatialIoError> {
    write_u64_le(w, canonical_json_bytes.len() as u64)?;
    w.write_all(canonical_json_bytes)?;
    Ok(())
}

fn write_len_prefixed_str<W: Write>(w: &mut W, s: &str) -> Result<(), SpatialIoError> {
    let len = u32::try_from(s.len()).map_err(|_| {
        SpatialIoError::DimensionMismatch("string length does not fit u32".to_string())
    })?;
    write_u32_le(w, len)?;
    w.write_all(s.as_bytes())?;
    Ok(())
}

fn write_f32_slice_le<W: Write>(w: &mut W, values: &[f32]) -> Result<(), SpatialIoError> {
    for &v in values {
        w.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

fn write_u32_slice_le<W: Write>(w: &mut W, values: &[u32]) -> Result<(), SpatialIoError> {
    for &v in values {
        write_u32_le(w, v)?;
    }
    Ok(())
}

fn write_u64_slice_le<W: Write>(w: &mut W, values: &[u64]) -> Result<(), SpatialIoError> {
    for &v in values {
        write_u64_le(w, v)?;
    }
    Ok(())
}

fn write_u16_le<W: Write>(w: &mut W, v: u16) -> Result<(), SpatialIoError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_u32_le<W: Write>(w: &mut W, v: u32) -> Result<(), SpatialIoError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_u64_le<W: Write>(w: &mut W, v: u64) -> Result<(), SpatialIoError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_zeros<W: Write>(w: &mut W, n: usize) -> Result<(), SpatialIoError> {
    const ZERO_BLOCK: [u8; 64] = [0_u8; 64];
    let mut remaining = n;
    while remaining > 0 {
        let chunk = remaining.min(ZERO_BLOCK.len());
        w.write_all(&ZERO_BLOCK[..chunk])?;
        remaining -= chunk;
    }
    Ok(())
}

fn coord_system_to_u8(coord: CoordSystem) -> u8 {
    match coord {
        CoordSystem::Grid => 0,
        CoordSystem::Pixel => 1,
        CoordSystem::Micron => 2,
    }
}

#[allow(dead_code)]
fn _endianness_guard() -> Result<(), SpatialIoError> {
    if ENDIAN_LITTLE != 1 {
        return Err(SpatialIoError::UnsupportedFormat(
            "only little-endian format is supported".to_string(),
        ));
    }
    Ok(())
}
