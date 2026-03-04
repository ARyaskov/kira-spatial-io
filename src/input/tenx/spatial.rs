use std::collections::HashMap;
use std::path::Path;

use bitvec::vec::BitVec;
use csv::ReaderBuilder;

use crate::config::LoadConfig;
use crate::error::SpatialIoError;
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
    let mut csv = ReaderBuilder::new().has_headers(false).from_reader(reader);

    let mut first_record = true;
    for result in csv.records() {
        let record = result.map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!("invalid spatial csv record: {e}"))
        })?;

        if record.len() < 6 {
            return Err(SpatialIoError::UnsupportedFormat(
                "spatial csv row has fewer than 6 columns".to_string(),
            ));
        }

        let barcode = record.get(0).unwrap_or_default();
        if first_record {
            first_record = false;
            if barcode == "barcode" {
                continue;
            }
        }

        let in_tissue = parse_u8(record.get(1).unwrap_or_default(), "in_tissue")?;
        let array_row = parse_u32(record.get(2).unwrap_or_default(), "array_row")?;
        let array_col = parse_u32(record.get(3).unwrap_or_default(), "array_col")?;
        let pxl_row = parse_f32(record.get(4).unwrap_or_default(), "pxl_row_in_fullres")?;
        let pxl_col = parse_f32(record.get(5).unwrap_or_default(), "pxl_col_in_fullres")?;

        let Some(&idx) = barcode_to_idx.get(barcode) else {
            if cfg.validate_strict {
                return Err(SpatialIoError::DimensionMismatch(format!(
                    "barcode in spatial not in barcodes: {barcode}"
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

fn parse_u8(s: &str, field: &str) -> Result<u8, SpatialIoError> {
    s.parse::<u8>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} value in spatial csv: {s}"))
    })
}

fn parse_u32(s: &str, field: &str) -> Result<u32, SpatialIoError> {
    s.parse::<u32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} value in spatial csv: {s}"))
    })
}

fn parse_f32(s: &str, field: &str) -> Result<f32, SpatialIoError> {
    s.parse::<f32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} value in spatial csv: {s}"))
    })
}
