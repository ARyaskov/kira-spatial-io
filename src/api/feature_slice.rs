//! Targeted single-gene loader from Visium HD `feature_slice.h5`.

use std::collections::HashMap;
use std::path::Path;

use hdf5::File;
use ndarray::s;

use crate::config::LoadConfig;
use crate::determinism::sort::sort_bins;
use crate::error::SpatialIoError;
use crate::input::h5::discover;
use crate::input::h5::strings::read_string_dataset_any;
use crate::input::parquet::barcode_mapping::load_barcode_mapping_parquet;
use crate::input::parquet::spatial_mapping::load_spatial_domain_from_mapping_table;
use crate::model::spatial_domain::SpatialDomain;

/// Targeted feature-slice load for a single gene from Visium HD `feature_slice.h5`.
#[derive(Debug)]
pub struct FeatureSliceGene {
    /// Canonical spatial domain sorted by grid row/col.
    pub spatial_domain: SpatialDomain,
    /// Sparse gene values scattered into `spatial_domain` order.
    pub values: Vec<f32>,
}

/// Loads one gene from a Visium HD `feature_slice.h5` using parquet barcode mappings.
pub fn load_feature_slice_gene<P: AsRef<Path>>(
    path: P,
    gene_name: &str,
    cfg: LoadConfig,
) -> Result<FeatureSliceGene, SpatialIoError> {
    let paths = discover::discover_h5_paths(path.as_ref())?;
    let parquet_path = match &paths.spatial_input {
        discover::SpatialInput::Csv(_) => {
            return Err(SpatialIoError::UnsupportedFormat(
                "feature_slice fallback requires barcode_mappings.parquet".to_string(),
            ));
        }
        discover::SpatialInput::Parquet(path) => path.clone(),
    };

    let h5 = File::open(&paths.h5_path).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!(
            "failed to open h5 file {}: {e}",
            paths.h5_path.display()
        ))
    })?;

    let feature_index = resolve_feature_index(&h5, gene_name)?;

    let slice_base = format!("/feature_slices/{feature_index}");
    let row_ds = h5.dataset(&(slice_base.clone() + "/row")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {slice_base}/row dataset for gene {gene_name}"
        ))
    })?;
    let col_ds = h5.dataset(&(slice_base.clone() + "/col")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {slice_base}/col dataset for gene {gene_name}"
        ))
    })?;
    let data_ds = h5.dataset(&(slice_base.clone() + "/data")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {slice_base}/data dataset for gene {gene_name}"
        ))
    })?;

    let rows = read_u32_dataset(&row_ds, &(slice_base.clone() + "/row"))?;
    let cols = read_u32_dataset(&col_ds, &(slice_base.clone() + "/col"))?;
    let counts = read_u32_dataset(&data_ds, &(slice_base.clone() + "/data"))?;

    if rows.len() != cols.len() || rows.len() != counts.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "feature slice arrays length mismatch for {gene_name}: row={}, col={}, data={}",
            rows.len(),
            cols.len(),
            counts.len()
        )));
    }

    let (table, summary) = load_barcode_mapping_parquet(&parquet_path, &cfg)?;
    if !summary.has_grid && !summary.has_xy {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "parquet mapping {} does not contain usable grid or xy coordinates",
            parquet_path.display()
        )));
    }

    let mut spatial_domain = load_spatial_domain_from_mapping_table(&table, &cfg)?;
    sort_bins(&mut spatial_domain)?;

    let grid_row = spatial_domain.grid_row.as_ref().ok_or_else(|| {
        SpatialIoError::DimensionMismatch("grid_row missing after mapping load".to_string())
    })?;
    let grid_col = spatial_domain.grid_col.as_ref().ok_or_else(|| {
        SpatialIoError::DimensionMismatch("grid_col missing after mapping load".to_string())
    })?;

    let mut coord_to_idx = HashMap::<(u32, u32), usize>::with_capacity(spatial_domain.len());
    for (idx, (&r, &c)) in grid_row.iter().zip(grid_col.iter()).enumerate() {
        coord_to_idx.entry((r, c)).or_insert(idx);
    }

    let mut values = vec![0.0_f32; spatial_domain.len()];
    let mut missing_coords = 0_u64;
    for ((&r, &c), &v) in rows.iter().zip(cols.iter()).zip(counts.iter()) {
        if let Some(&idx) = coord_to_idx.get(&(r, c)) {
            values[idx] += v as f32;
        } else {
            missing_coords += 1;
        }
    }

    if cfg.validate_strict && missing_coords > 0 {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "feature slice entries without coordinate match in parquet mapping: {missing_coords}"
        )));
    }

    Ok(FeatureSliceGene {
        spatial_domain,
        values,
    })
}

fn resolve_feature_index(file: &hdf5::File, gene_name: &str) -> Result<usize, SpatialIoError> {
    let names_ds = file.dataset("/features/name").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /features/name dataset".to_string())
    })?;
    let feature_names = read_string_dataset_any(&names_ds, "/features/name")?;
    feature_names
        .iter()
        .position(|name| name == gene_name)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!(
                "gene not found in /features/name: {gene_name}"
            ))
        })
}

fn read_u32_dataset(ds: &hdf5::Dataset, path: &str) -> Result<Vec<u32>, SpatialIoError> {
    if let Ok(v) = ds.read_raw::<u32>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_raw::<u64>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x > u32::MAX as u64 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "value too large for u32 in {path}"
                )));
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    if let Ok(v) = ds.read_raw::<i32>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x < 0 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "negative value in {path}"
                )));
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    if let Ok(v) = ds.read_raw::<i64>() {
        let mut out = Vec::with_capacity(v.len());
        for x in v {
            if x < 0 || x > u32::MAX as i64 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "value out of u32 range in {path}"
                )));
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    let len_guess = ds.shape().into_iter().product::<usize>();
    if let Ok(v) = ds.read_slice_1d::<u32, _>(s![..len_guess]) {
        return Ok(v.to_vec());
    }
    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported dtype in {path}"
    )))
}
