//! Reader for Visium `tissue_positions[_list].csv`.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use bitvec::vec::BitVec;

use crate::config::LoadConfig;
use crate::error::{IoPathExt, SpatialIoError};
use crate::input::util::open_text_maybe_gz;
use crate::model::coord::CoordSystem;
use crate::model::spatial_domain::SpatialDomain;

pub fn load_spatial_domain(
    csv_path: &Path,
    barcodes: &[String],
    cfg: &LoadConfig,
) -> Result<SpatialDomain, SpatialIoError> {
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

    let reader = open_text_maybe_gz(csv_path)?;
    let mut first_record = true;
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.io_path(csv_path)?;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() {
            continue;
        }

        let mut cols = trimmed.split(',');
        let barcode = cols.next().unwrap_or("");
        let c1 = cols.next().unwrap_or("");
        let c2 = cols.next().unwrap_or("");
        let c3 = cols.next().unwrap_or("");
        let c4 = cols.next().unwrap_or("");
        let c5 = cols.next().unwrap_or("");
        if c5.is_empty() {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "{}: row {} has fewer than 6 columns",
                csv_path.display(),
                lineno + 1
            )));
        }

        if first_record {
            first_record = false;
            if barcode == "barcode" {
                continue;
            }
        }

        let in_tissue = parse_u8(c1, "in_tissue", csv_path, lineno)?;
        let array_row = parse_u32(c2, "array_row", csv_path, lineno)?;
        let array_col = parse_u32(c3, "array_col", csv_path, lineno)?;
        let pxl_row = parse_f32(c4, "pxl_row_in_fullres", csv_path, lineno)?;
        let pxl_col = parse_f32(c5, "pxl_col_in_fullres", csv_path, lineno)?;

        let Some(&idx) = barcode_to_idx.get(barcode) else {
            if cfg.validate_strict {
                return Err(SpatialIoError::DimensionMismatch(format!(
                    "{}: barcode in spatial not in barcodes: {barcode}",
                    csv_path.display()
                )));
            }
            continue;
        };

        x[idx] = pxl_col;
        y[idx] = pxl_row;
        grid_row[idx] = array_row;
        grid_col[idx] = array_col;
        tissue_mask.set(idx, in_tissue == 1);
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

fn parse_u8(
    s: &str,
    field: &str,
    path: &Path,
    lineno: usize,
) -> Result<u8, SpatialIoError> {
    s.parse::<u8>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "{}: invalid {field} at line {}: {s}",
            path.display(),
            lineno + 1
        ))
    })
}

fn parse_u32(
    s: &str,
    field: &str,
    path: &Path,
    lineno: usize,
) -> Result<u32, SpatialIoError> {
    s.parse::<u32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "{}: invalid {field} at line {}: {s}",
            path.display(),
            lineno + 1
        ))
    })
}

fn parse_f32(
    s: &str,
    field: &str,
    path: &Path,
    lineno: usize,
) -> Result<f32, SpatialIoError> {
    s.parse::<f32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "{}: invalid {field} at line {}: {s}",
            path.display(),
            lineno + 1
        ))
    })
}
