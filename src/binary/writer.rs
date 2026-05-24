//! Deterministic writer for `.kira-spatial.bin` (format v2).

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use crc32fast::Hasher as Crc32;

use crate::api::dataset::Dataset;
use crate::binary::bitmask::tissue_mask_to_u64_words;
use crate::binary::compress::zstd_compress;
use crate::binary::format::{
    CSR_FLAG_INDPTR_U32, HEADER_SIZE, Header, MAGIC, MANDATORY_SECTION_IDS, SECTION_ENTRY_SIZE,
    SECTION_FLAG_ZSTD, SECTION_ID_CSR, SECTION_ID_FEATURE_TABLE, SECTION_ID_META_CORE,
    SECTION_ID_META_JSON, SECTION_ID_SPATIAL_DOMAIN, SectionEntry,
};
use crate::binary::hash::compute_dataset_hash;
use crate::config::CompressionPolicy;
use crate::determinism::json::{canonicalize_json, write_canonical_json};
use crate::error::{IoPathExt, SpatialIoError};
use crate::model::{
    coord::CoordSystem,
    csr::{BinsCsr, Indptr},
    metadata::DatasetMetaCore,
    spatial_domain::SpatialDomain,
};

/// Writes a dataset to a deterministic `.kira-spatial.bin` file.
pub fn write_kira_bin<P: AsRef<Path>>(p: P, ds: &Dataset) -> Result<(), SpatialIoError> {
    write_kira_bin_with_compression(p, ds, CompressionPolicy::None)
}

/// Writes a dataset to a deterministic `.kira-spatial.bin` file with the given compression policy.
pub fn write_kira_bin_with_compression<P: AsRef<Path>>(
    p: P,
    ds: &Dataset,
    compression: CompressionPolicy,
) -> Result<(), SpatialIoError> {
    let path = p.as_ref();

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

    let mut spatial_buf = Vec::new();
    write_spatial_domain_section(&mut spatial_buf, ds.spatial_domain())?;

    let mut csr_buf = Vec::new();
    write_csr_section(&mut csr_buf, ds.expression_csr())?;

    let mut features_buf = Vec::new();
    write_feature_table_section(&mut features_buf, ds.features())?;

    let mut meta_buf = Vec::new();
    write_metadata_core_section(&mut meta_buf, ds.metadata_core())?;

    let mut json_buf = Vec::new();
    write_metadata_json_section(&mut json_buf, &canonical_json_bytes)?;

    let raw_sections: [(u16, Vec<u8>); 5] = [
        (SECTION_ID_SPATIAL_DOMAIN, spatial_buf),
        (SECTION_ID_CSR, csr_buf),
        (SECTION_ID_FEATURE_TABLE, features_buf),
        (SECTION_ID_META_CORE, meta_buf),
        (SECTION_ID_META_JSON, json_buf),
    ];

    let compress_level = match compression {
        CompressionPolicy::None => None,
        CompressionPolicy::Zstd(level) => Some(level),
    };

    let mut prepared: Vec<(u16, u16, Vec<u8>)> = Vec::with_capacity(raw_sections.len());
    for (id, raw) in raw_sections {
        match compress_level {
            Some(level) => {
                let compressed = zstd_compress(&raw, level)?;
                prepared.push((id, SECTION_FLAG_ZSTD, compressed));
            }
            None => prepared.push((id, 0, raw)),
        }
    }

    let section_count = prepared.len() as u16;
    let file = File::create(path).map_err(|e| SpatialIoError::io_at(path, e))?;
    let mut writer = BufWriter::new(file);

    write_header(
        &mut writer,
        &Header::new(section_count, dataset_hash),
    )?;
    let placeholder_table = vec![SectionEntry::new(0, 0, 0, 0, 0); prepared.len()];
    write_section_table(&mut writer, &placeholder_table).io_path(path)?;

    let mut entries: Vec<SectionEntry> = Vec::with_capacity(prepared.len());
    for (id, flags, payload) in &prepared {
        pad_to_64(&mut writer).io_path(path)?;
        let offset = writer.stream_position().map_err(|e| SpatialIoError::io_at(path, e))?;
        writer.write_all(payload).map_err(|e| SpatialIoError::io_at(path, e))?;
        let crc = compute_crc32c(payload);
        entries.push(SectionEntry::new(
            *id,
            *flags,
            offset,
            payload.len() as u64,
            crc,
        ));
    }

    writer.flush().map_err(|e| SpatialIoError::io_at(path, e))?;
    writer
        .seek(SeekFrom::Start(0))
        .map_err(|e| SpatialIoError::io_at(path, e))?;
    write_header(
        &mut writer,
        &Header::new(section_count, dataset_hash),
    )?;
    write_section_table(&mut writer, &entries).io_path(path)?;
    writer.flush().map_err(|e| SpatialIoError::io_at(path, e))?;
    writer
        .get_ref()
        .sync_all()
        .map_err(|e| SpatialIoError::io_at(path, e))?;

    Ok(())
}

fn compute_crc32c(bytes: &[u8]) -> u32 {
    let mut h = Crc32::new();
    h.update(bytes);
    h.finalize()
}

/// Pads writer position to next 64-byte boundary.
fn pad_to_64<W: Write + Seek>(w: &mut W) -> std::io::Result<()> {
    let pos = w.stream_position()?;
    let rem = (pos % 64) as usize;
    if rem == 0 {
        return Ok(());
    }
    let need = 64 - rem;
    let zeros = [0_u8; 64];
    w.write_all(&zeros[..need])
}

fn write_header<W: Write>(w: &mut W, header: &Header) -> Result<(), SpatialIoError> {
    let mut buf = [0_u8; HEADER_SIZE as usize];
    buf[0..8].copy_from_slice(&MAGIC);
    buf[8..10].copy_from_slice(&header.version.to_le_bytes());
    buf[10..12].copy_from_slice(&header.section_count.to_le_bytes());
    buf[16..32].copy_from_slice(&header.dataset_hash);
    w.write_all(&buf)?;
    Ok(())
}

fn write_section_table<W: Write>(w: &mut W, entries: &[SectionEntry]) -> std::io::Result<()> {
    debug_assert_eq!(SECTION_ENTRY_SIZE, 24);
    let mut buf = vec![0_u8; entries.len() * SECTION_ENTRY_SIZE as usize];
    for (i, entry) in entries.iter().enumerate() {
        let base = i * SECTION_ENTRY_SIZE as usize;
        buf[base..base + 2].copy_from_slice(&entry.id.to_le_bytes());
        buf[base + 2..base + 4].copy_from_slice(&entry.flags.to_le_bytes());
        buf[base + 4..base + 12].copy_from_slice(&entry.offset.to_le_bytes());
        buf[base + 12..base + 20].copy_from_slice(&entry.length.to_le_bytes());
        buf[base + 20..base + 24].copy_from_slice(&entry.crc32c.to_le_bytes());
    }
    w.write_all(&buf)
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

    let n_bins_u32 = u32::try_from(n_bins).map_err(|_| {
        SpatialIoError::DimensionMismatch("n_bins does not fit u32".to_string())
    })?;
    w.write_all(&n_bins_u32.to_le_bytes())?;
    w.write_all(&[coord_system_to_u8(domain.coord_system)])?;
    w.write_all(&[domain.bin_level])?;

    let mut flags = 0_u16;
    if has_grid {
        flags |= 1;
    }
    w.write_all(&flags.to_le_bytes())?;

    w.write_all(bytemuck::cast_slice(&domain.x))?;
    w.write_all(bytemuck::cast_slice(&domain.y))?;

    if let (Some(rows), Some(cols)) = (&domain.grid_row, &domain.grid_col) {
        w.write_all(bytemuck::cast_slice(rows))?;
        w.write_all(bytemuck::cast_slice(cols))?;
    }

    w.write_all(bytemuck::cast_slice(&domain.bin_id))?;

    let words = tissue_mask_to_u64_words(&domain.tissue_mask);
    w.write_all(&(domain.tissue_mask.len() as u64).to_le_bytes())?;
    w.write_all(&((words.len() * 8) as u64).to_le_bytes())?;
    w.write_all(bytemuck::cast_slice(&words))?;

    Ok(())
}

fn write_csr_section<W: Write>(w: &mut W, csr: &BinsCsr) -> Result<(), SpatialIoError> {
    if csr.indptr.len() != csr.n_bins as usize + 1 {
        return Err(SpatialIoError::InvalidCsr(
            "indptr length does not match n_bins + 1".to_string(),
        ));
    }
    if csr.indptr.last().unwrap_or(1) != csr.nnz {
        return Err(SpatialIoError::InvalidCsr(
            "indptr[last] does not match nnz".to_string(),
        ));
    }
    if csr.indices.len() != csr.nnz as usize || csr.data.len() != csr.nnz as usize {
        return Err(SpatialIoError::InvalidCsr(
            "indices/data length does not match nnz".to_string(),
        ));
    }

    w.write_all(&csr.n_bins.to_le_bytes())?;
    w.write_all(&csr.n_genes.to_le_bytes())?;
    w.write_all(&csr.nnz.to_le_bytes())?;

    let mut flags = 0_u8;
    if csr.indptr.is_u32() {
        flags |= CSR_FLAG_INDPTR_U32;
    }
    w.write_all(&[u8::from(csr.normalized)])?;
    w.write_all(&[flags])?;
    w.write_all(&[0_u8; 6])?;

    match &csr.indptr {
        Indptr::U32(v) => w.write_all(bytemuck::cast_slice(v))?,
        Indptr::U64(v) => w.write_all(bytemuck::cast_slice(v))?,
    }
    w.write_all(bytemuck::cast_slice(&csr.indices))?;
    w.write_all(bytemuck::cast_slice(&csr.data))?;
    Ok(())
}

fn write_feature_table_section<W: Write>(
    w: &mut W,
    table: &crate::model::features::FeatureTable,
) -> Result<(), SpatialIoError> {
    let n_genes = u32::try_from(table.rows.len()).map_err(|_| {
        SpatialIoError::DimensionMismatch("n_genes does not fit u32".to_string())
    })?;
    w.write_all(&n_genes.to_le_bytes())?;

    for (i, row) in table.rows.iter().enumerate() {
        if row.gene_id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table gene_id is not canonical at row {i}"
            )));
        }
        w.write_all(&row.gene_id.to_le_bytes())?;
        write_len_prefixed_str(w, &row.feature_id)?;
        write_len_prefixed_str(w, &row.gene_name)?;
        write_len_prefixed_str(w, &row.feature_type)?;
    }

    Ok(())
}

fn write_metadata_core_section<W: Write>(
    w: &mut W,
    meta: &DatasetMetaCore,
) -> Result<(), SpatialIoError> {
    write_len_prefixed_str(w, &meta.dataset_name)?;
    write_len_prefixed_str(w, &meta.source_format)?;
    w.write_all(&[meta.bin_level])?;
    w.write_all(&meta.n_bins.to_le_bytes())?;
    w.write_all(&meta.n_genes.to_le_bytes())?;
    w.write_all(&meta.nnz.to_le_bytes())?;
    w.write_all(&[coord_system_to_u8(meta.coord_system)])?;
    w.write_all(&[u8::from(meta.normalized)])?;
    Ok(())
}

fn write_metadata_json_section<W: Write>(
    w: &mut W,
    canonical_json_bytes: &[u8],
) -> Result<(), SpatialIoError> {
    w.write_all(&(canonical_json_bytes.len() as u64).to_le_bytes())?;
    w.write_all(canonical_json_bytes)?;
    Ok(())
}

fn write_len_prefixed_str<W: Write>(w: &mut W, s: &str) -> Result<(), SpatialIoError> {
    let len = u32::try_from(s.len()).map_err(|_| {
        SpatialIoError::DimensionMismatch("string length does not fit u32".to_string())
    })?;
    w.write_all(&len.to_le_bytes())?;
    w.write_all(s.as_bytes())?;
    Ok(())
}

fn coord_system_to_u8(coord: CoordSystem) -> u8 {
    match coord {
        CoordSystem::Grid => 0,
        CoordSystem::Pixel => 1,
        CoordSystem::Micron => 2,
    }
}

#[doc(hidden)]
#[allow(dead_code)]
const _: () = {
    if MANDATORY_SECTION_IDS.len() != 5 {
        panic!("mandatory section id count must remain 5");
    }
};
