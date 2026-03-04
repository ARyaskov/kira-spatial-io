use std::collections::HashMap;
use std::path::Path;

use bitvec::vec::BitVec;

use crate::config::LoadConfig;
use crate::error::SpatialIoError;
use crate::input::parquet::barcode_mapping::load_barcode_mapping_parquet;
use crate::model::coord::CoordSystem;
use crate::model::mapping::{BarcodeMappingRow, BarcodeMappingTable};
use crate::model::spatial_domain::SpatialDomain;

pub(crate) fn load_spatial_domain_from_mapping_parquet(
    path: &Path,
    barcodes: &[String],
    cfg: &LoadConfig,
) -> Result<SpatialDomain, SpatialIoError> {
    let (table, summary) = load_barcode_mapping_parquet(path, cfg)?;
    if !summary.has_grid && !summary.has_xy {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "parquet mapping {} does not contain usable grid or xy coordinates",
            path.display()
        )));
    }

    let mut barcode_to_idx = HashMap::with_capacity(barcodes.len());
    for (idx, barcode) in barcodes.iter().enumerate() {
        barcode_to_idx.insert(barcode.as_str(), idx);
    }

    let n = barcodes.len();
    let mut x = vec![0.0_f32; n];
    let mut y = vec![0.0_f32; n];
    let mut grid_row = vec![0_u32; n];
    let mut grid_col = vec![0_u32; n];
    let mut tissue_mask = BitVec::repeat(false, n);
    let mut seen = vec![false; n];

    for row in &table.rows {
        let Some(&idx) = barcode_to_idx.get(row.barcode.as_str()) else {
            continue;
        };
        if seen[idx] {
            continue;
        }

        let Some((xv, yv, r, c)) = extract_coords(row)? else {
            continue;
        };

        x[idx] = xv;
        y[idx] = yv;
        grid_row[idx] = r;
        grid_col[idx] = c;
        tissue_mask.set(idx, true);
        seen[idx] = true;
    }

    if cfg.validate_strict
        && let Some((missing_idx, _)) = seen.iter().enumerate().find(|(_, is_seen)| !**is_seen)
    {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "missing coordinates for barcode: {}",
            barcodes[missing_idx]
        )));
    }

    let bin_id: Vec<u32> = (0..(n as u32)).collect();
    SpatialDomain::new(
        x,
        y,
        Some(grid_row),
        Some(grid_col),
        bin_id,
        tissue_mask,
        CoordSystem::Pixel,
        cfg.bin_level.unwrap_or(0),
    )
}

fn extract_coords(row: &BarcodeMappingRow) -> Result<Option<(f32, f32, u32, u32)>, SpatialIoError> {
    match (row.x, row.y, row.grid_row, row.grid_col) {
        (Some(x), Some(y), Some(r), Some(c)) => Ok(Some((x, y, r, c))),
        (None, None, Some(r), Some(c)) => Ok(Some((c as f32, r as f32, r, c))),
        (Some(x), Some(y), None, None) => {
            let r = f32_to_u32_round_nonneg(y, "y")?;
            let c = f32_to_u32_round_nonneg(x, "x")?;
            Ok(Some((x, y, r, c)))
        }
        _ => Ok(None),
    }
}

fn f32_to_u32_round_nonneg(value: f32, field: &str) -> Result<u32, SpatialIoError> {
    if !value.is_finite() || value < 0.0 {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "invalid {} value in parquet mapping: {}",
            field, value
        )));
    }
    let rounded = value.round();
    if rounded > u32::MAX as f32 {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "{} value in parquet mapping out of u32 range: {}",
            field, value
        )));
    }
    Ok(rounded as u32)
}

pub(crate) fn load_spatial_domain_from_mapping_table(
    table: &BarcodeMappingTable,
    cfg: &LoadConfig,
) -> Result<SpatialDomain, SpatialIoError> {
    let mut x = Vec::<f32>::new();
    let mut y = Vec::<f32>::new();
    let mut grid_row = Vec::<u32>::new();
    let mut grid_col = Vec::<u32>::new();
    let mut tissue_mask = BitVec::new();

    for row in &table.rows {
        let Some((xv, yv, r, c)) = extract_coords(row)? else {
            continue;
        };
        x.push(xv);
        y.push(yv);
        grid_row.push(r);
        grid_col.push(c);
        tissue_mask.push(true);
    }

    if x.is_empty() {
        return Err(SpatialIoError::UnsupportedFormat(
            "parquet mapping does not contain usable coordinates".to_string(),
        ));
    }

    if cfg.validate_strict && x.len() != table.rows.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "missing coordinates in parquet mapping: {} rows without usable coordinates",
            table.rows.len() - x.len()
        )));
    }

    let bin_id: Vec<u32> = (0..(x.len() as u32)).collect();
    SpatialDomain::new(
        x,
        y,
        Some(grid_row),
        Some(grid_col),
        bin_id,
        tissue_mask,
        CoordSystem::Pixel,
        cfg.bin_level.unwrap_or(0),
    )
}
