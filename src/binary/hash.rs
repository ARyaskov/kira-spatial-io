//! Canonical dataset hash (leading 16 bytes of BLAKE3 over canonical section streams).

use blake3::Hasher;

use crate::binary::bitmask::tissue_mask_to_u64_words;
use crate::binary::format::{
    SECTION_ID_CSR, SECTION_ID_FEATURE_TABLE, SECTION_ID_META_CORE, SECTION_ID_META_JSON,
    SECTION_ID_SPATIAL_DOMAIN,
};
use crate::error::SpatialIoError;
use crate::model::{
    coord::CoordSystem, csr::BinsCsr, features::FeatureTable, metadata::DatasetMetaCore,
    spatial_domain::SpatialDomain,
};

/// Length of the canonical dataset hash in bytes.
pub const DATASET_HASH_BYTES: usize = 16;

/// Computes the canonical dataset hash from normalized payloads.
pub fn compute_dataset_hash(
    spatial: &SpatialDomain,
    csr: &BinsCsr,
    features: &FeatureTable,
    meta: &DatasetMetaCore,
    canonical_json_bytes: &[u8],
) -> Result<[u8; DATASET_HASH_BYTES], SpatialIoError> {
    let mut hasher = Hasher::new();

    hash_u16(&mut hasher, SECTION_ID_SPATIAL_DOMAIN);
    hash_spatial_domain(&mut hasher, spatial)?;

    hash_u16(&mut hasher, SECTION_ID_CSR);
    hash_csr(&mut hasher, csr)?;

    hash_u16(&mut hasher, SECTION_ID_FEATURE_TABLE);
    hash_feature_table(&mut hasher, features)?;

    hash_u16(&mut hasher, SECTION_ID_META_CORE);
    hash_metadata_core(&mut hasher, meta)?;

    hash_u16(&mut hasher, SECTION_ID_META_JSON);
    hash_u64(&mut hasher, canonical_json_bytes.len() as u64);
    hasher.update(canonical_json_bytes);

    let full = hasher.finalize();
    let mut out = [0_u8; DATASET_HASH_BYTES];
    out.copy_from_slice(&full.as_bytes()[..DATASET_HASH_BYTES]);
    Ok(out)
}

fn hash_spatial_domain(hasher: &mut Hasher, domain: &SpatialDomain) -> Result<(), SpatialIoError> {
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

    hash_u32(hasher, n_bins as u32);
    hash_u8(hasher, coord_system_to_u8(domain.coord_system));
    hash_u8(hasher, domain.bin_level);
    hash_u16(hasher, if has_grid { 1 } else { 0 });

    hasher.update(bytemuck::cast_slice(&domain.x));
    hasher.update(bytemuck::cast_slice(&domain.y));

    if let (Some(rows), Some(cols)) = (&domain.grid_row, &domain.grid_col) {
        hasher.update(bytemuck::cast_slice(rows));
        hasher.update(bytemuck::cast_slice(cols));
    }

    hasher.update(bytemuck::cast_slice(&domain.bin_id));

    let words = tissue_mask_to_u64_words(&domain.tissue_mask);
    hash_u64(hasher, domain.tissue_mask.len() as u64);
    hash_u64(hasher, (words.len() * 8) as u64);
    hasher.update(bytemuck::cast_slice(&words));

    Ok(())
}

fn hash_csr(hasher: &mut Hasher, csr: &BinsCsr) -> Result<(), SpatialIoError> {
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

    hash_u32(hasher, csr.n_bins);
    hash_u32(hasher, csr.n_genes);
    hash_u64(hasher, csr.nnz);
    hash_u8(hasher, u8::from(csr.normalized));

    let indptr_u64 = csr.indptr.to_u64_vec();
    hasher.update(bytemuck::cast_slice(&indptr_u64));
    hasher.update(bytemuck::cast_slice(&csr.indices));
    hasher.update(bytemuck::cast_slice(&csr.data));

    Ok(())
}

fn hash_feature_table(hasher: &mut Hasher, table: &FeatureTable) -> Result<(), SpatialIoError> {
    hash_u32(hasher, table.rows.len() as u32);
    for (i, row) in table.rows.iter().enumerate() {
        if row.gene_id != i as u32 {
            return Err(SpatialIoError::InvalidCsr(format!(
                "feature table gene_id is not canonical at row {i}"
            )));
        }
        hash_u32(hasher, row.gene_id);
        hash_len_prefixed_str(hasher, &row.feature_id)?;
        hash_len_prefixed_str(hasher, &row.gene_name)?;
        hash_len_prefixed_str(hasher, &row.feature_type)?;
    }
    Ok(())
}

fn hash_metadata_core(hasher: &mut Hasher, meta: &DatasetMetaCore) -> Result<(), SpatialIoError> {
    hash_len_prefixed_str(hasher, &meta.dataset_name)?;
    hash_len_prefixed_str(hasher, &meta.source_format)?;
    hash_u8(hasher, meta.bin_level);
    hash_u32(hasher, meta.n_bins);
    hash_u32(hasher, meta.n_genes);
    hash_u64(hasher, meta.nnz);
    hash_u8(hasher, coord_system_to_u8(meta.coord_system));
    hash_u8(hasher, u8::from(meta.normalized));
    Ok(())
}

fn hash_len_prefixed_str(hasher: &mut Hasher, s: &str) -> Result<(), SpatialIoError> {
    let len = u32::try_from(s.len()).map_err(|_| {
        SpatialIoError::DimensionMismatch("string length does not fit u32".to_string())
    })?;
    hash_u32(hasher, len);
    hasher.update(s.as_bytes());
    Ok(())
}

fn coord_system_to_u8(coord: CoordSystem) -> u8 {
    match coord {
        CoordSystem::Grid => 0,
        CoordSystem::Pixel => 1,
        CoordSystem::Micron => 2,
    }
}

fn hash_u8(hasher: &mut Hasher, v: u8) {
    hasher.update(&[v]);
}

fn hash_u16(hasher: &mut Hasher, v: u16) {
    hasher.update(&v.to_le_bytes());
}

fn hash_u32(hasher: &mut Hasher, v: u32) {
    hasher.update(&v.to_le_bytes());
}

fn hash_u64(hasher: &mut Hasher, v: u64) {
    hasher.update(&v.to_le_bytes());
}
