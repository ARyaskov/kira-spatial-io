use std::cmp::Ordering;
use std::fs::File;
use std::path::Path;

use arrow2::array::{Array, BinaryArray, PrimitiveArray, Utf8Array};
use arrow2::datatypes::DataType;
use arrow2::io::parquet::read::{FileReader, infer_schema, read_metadata};

use crate::config::LoadConfig;
use crate::determinism::float::ensure_f32_finite_nonneg;
use crate::error::SpatialIoError;
use crate::input::tenx::mtx::ensure_budget;
use crate::model::mapping::{BarcodeMappingRow, BarcodeMappingTable};

const MAPPING_ROW_ESTIMATE_BYTES: u64 = 160;

#[derive(Clone, Debug)]
pub(crate) struct MappingSummary {
    pub row_count: u64,
    pub has_cell_id: bool,
    pub has_grid: bool,
    pub has_xy: bool,
    pub mapping_has_duplicates: bool,
}

pub(crate) fn load_barcode_mapping_parquet(
    path: &Path,
    cfg: &LoadConfig,
) -> Result<(BarcodeMappingTable, MappingSummary), SpatialIoError> {
    let mut file = File::open(path)?;
    let metadata = read_metadata(&mut file).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!(
            "failed reading parquet metadata {}: {e}",
            path.display()
        ))
    })?;
    let schema = infer_schema(&metadata).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!(
            "failed inferring parquet schema {}: {e}",
            path.display()
        ))
    })?;

    let columns = schema
        .fields
        .iter()
        .map(|f| f.name.as_str())
        .collect::<Vec<_>>();
    let selected = discover_columns(&columns)?;

    let mut rows = Vec::new();
    let row_groups = metadata.row_groups.clone();
    let reader = FileReader::new(file, row_groups, schema, None, None, None);

    for chunk_result in reader {
        let chunk = chunk_result.map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!(
                "failed reading parquet row group {}: {e}",
                path.display()
            ))
        })?;

        for row_idx in 0..chunk.len() {
            let barcode = read_required_string(chunk.arrays()[selected.barcode].as_ref(), row_idx)?;
            let cell_id = selected
                .cell_id
                .and_then(|i| {
                    read_optional_u64(chunk.arrays()[i].as_ref(), row_idx, "cell_id").ok()
                })
                .flatten();
            let mut x = selected
                .x
                .map(|i| read_optional_f32(chunk.arrays()[i].as_ref(), row_idx, "x"))
                .transpose()?
                .flatten();
            let mut y = selected
                .y
                .map(|i| read_optional_f32(chunk.arrays()[i].as_ref(), row_idx, "y"))
                .transpose()?
                .flatten();

            let mut grid_row = selected
                .grid_row
                .map(|i| read_optional_u32(chunk.arrays()[i].as_ref(), row_idx, "grid_row"))
                .transpose()?
                .flatten();
            let mut grid_col = selected
                .grid_col
                .map(|i| read_optional_u32(chunk.arrays()[i].as_ref(), row_idx, "grid_col"))
                .transpose()?
                .flatten();

            if (grid_row.is_none() || grid_col.is_none())
                && let Some((r, c)) = parse_square_barcode_coords(&barcode)
            {
                if grid_row.is_none() {
                    grid_row = Some(r);
                }
                if grid_col.is_none() {
                    grid_col = Some(c);
                }
                if x.is_none() {
                    x = Some(c as f32);
                }
                if y.is_none() {
                    y = Some(r as f32);
                }
            }

            rows.push(BarcodeMappingRow {
                barcode,
                cell_id,
                grid_row,
                grid_col,
                x,
                y,
            });
        }

        let estimated = (rows.len() as u64)
            .checked_mul(MAPPING_ROW_ESTIMATE_BYTES)
            .ok_or_else(|| {
                SpatialIoError::MemoryLimitExceeded("parquet mapping size overflow".to_string())
            })?;
        ensure_budget(estimated, cfg, "parquet mapping too large")?;
    }

    rows.sort_unstable_by(compare_mapping_rows);

    let mapping_has_duplicates = rows.windows(2).any(|w| w[0].barcode == w[1].barcode);
    let summary = MappingSummary {
        row_count: rows.len() as u64,
        has_cell_id: rows.iter().any(|r| r.cell_id.is_some()),
        has_grid: rows
            .iter()
            .any(|r| r.grid_row.is_some() && r.grid_col.is_some()),
        has_xy: rows.iter().any(|r| r.x.is_some() && r.y.is_some()),
        mapping_has_duplicates,
    };

    Ok((BarcodeMappingTable::new(rows), summary))
}

#[derive(Debug, Clone, Copy)]
struct SelectedColumns {
    barcode: usize,
    cell_id: Option<usize>,
    grid_row: Option<usize>,
    grid_col: Option<usize>,
    x: Option<usize>,
    y: Option<usize>,
}

fn discover_columns(column_names: &[&str]) -> Result<SelectedColumns, SpatialIoError> {
    let normalized = column_names
        .iter()
        .map(|s| s.to_ascii_lowercase())
        .collect::<Vec<_>>();

    let barcode = find_first(
        &normalized,
        &[
            "barcode",
            "barcodes",
            "square_002um",
            "square_008um",
            "square_016um",
        ],
    )
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(
            "parquet barcode mapping missing required columns".to_string(),
        )
    })?;

    let cell_id = find_first(
        &normalized,
        &[
            "cell_id",
            "cell",
            "segmentation_label",
            "segmentation_id",
            "label",
        ],
    );

    let grid_row = find_first(
        &normalized,
        &["grid_row", "array_row", "row", "bin_row", "bin_y"],
    );
    let grid_col = find_first(
        &normalized,
        &["grid_col", "array_col", "col", "bin_col", "bin_x"],
    );

    let x = find_first(
        &normalized,
        &["x", "pxl_col_in_fullres", "pixel_x", "coord_x"],
    );
    let y = find_first(
        &normalized,
        &["y", "pxl_row_in_fullres", "pixel_y", "coord_y"],
    );

    let has_required_mapping = cell_id.is_some()
        || (grid_row.is_some() && grid_col.is_some())
        || (x.is_some() && y.is_some());

    if !has_required_mapping {
        return Err(SpatialIoError::UnsupportedFormat(
            "parquet barcode mapping missing required columns".to_string(),
        ));
    }

    Ok(SelectedColumns {
        barcode,
        cell_id,
        grid_row,
        grid_col,
        x,
        y,
    })
}

fn parse_square_barcode_coords(value: &str) -> Option<(u32, u32)> {
    // Example: s_008um_00012_00345-1
    if !value.starts_with("s_") {
        return None;
    }
    let mut parts = value.split('_');
    let _prefix = parts.next()?;
    let _bin = parts.next()?;
    let row = parts.next()?;
    let col_with_suffix = parts.next()?;
    let col = col_with_suffix.split('-').next()?;

    let r = row.parse::<u32>().ok()?;
    let c = col.parse::<u32>().ok()?;
    Some((r, c))
}

fn find_first(normalized: &[String], candidates: &[&str]) -> Option<usize> {
    candidates
        .iter()
        .find_map(|candidate| normalized.iter().position(|name| name == candidate))
}

fn read_required_string(array: &dyn Array, row_idx: usize) -> Result<String, SpatialIoError> {
    let value = read_optional_string(array, row_idx)?.ok_or_else(|| {
        SpatialIoError::UnsupportedFormat("parquet barcode contains null value".to_string())
    })?;
    if value.is_empty() {
        return Err(SpatialIoError::UnsupportedFormat(
            "parquet barcode contains empty value".to_string(),
        ));
    }
    Ok(value)
}

fn read_optional_string(
    array: &dyn Array,
    row_idx: usize,
) -> Result<Option<String>, SpatialIoError> {
    if array.is_null(row_idx) {
        return Ok(None);
    }

    if let Some(a) = array.as_any().downcast_ref::<Utf8Array<i32>>() {
        return Ok(Some(a.value(row_idx).to_string()));
    }
    if let Some(a) = array.as_any().downcast_ref::<Utf8Array<i64>>() {
        return Ok(Some(a.value(row_idx).to_string()));
    }
    if let Some(a) = array.as_any().downcast_ref::<BinaryArray<i32>>() {
        return Ok(Some(String::from_utf8(a.value(row_idx).to_vec()).map_err(
            |_| {
                SpatialIoError::UnsupportedFormat(
                    "parquet barcode column is not valid UTF-8".to_string(),
                )
            },
        )?));
    }
    if let Some(a) = array.as_any().downcast_ref::<BinaryArray<i64>>() {
        return Ok(Some(String::from_utf8(a.value(row_idx).to_vec()).map_err(
            |_| {
                SpatialIoError::UnsupportedFormat(
                    "parquet barcode column is not valid UTF-8".to_string(),
                )
            },
        )?));
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported parquet string column type: {:?}",
        array.data_type()
    )))
}

fn read_optional_u64(
    array: &dyn Array,
    row_idx: usize,
    column: &str,
) -> Result<Option<u64>, SpatialIoError> {
    if array.is_null(row_idx) {
        return Ok(None);
    }

    let value = match array.data_type() {
        DataType::UInt8 => read_primitive::<u8>(array, row_idx).map(u64::from),
        DataType::UInt16 => read_primitive::<u16>(array, row_idx).map(u64::from),
        DataType::UInt32 => read_primitive::<u32>(array, row_idx).map(u64::from),
        DataType::UInt64 => read_primitive::<u64>(array, row_idx),
        DataType::Int8 => read_primitive::<i8>(array, row_idx).and_then(to_u64_i128),
        DataType::Int16 => read_primitive::<i16>(array, row_idx).and_then(to_u64_i128),
        DataType::Int32 => read_primitive::<i32>(array, row_idx).and_then(to_u64_i128),
        DataType::Int64 => read_primitive::<i64>(array, row_idx).and_then(to_u64_i128),
        _ => None,
    }
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(format!(
            "unsupported parquet numeric column for {}: {:?}",
            column,
            array.data_type()
        ))
    })?;

    Ok(Some(value))
}

fn read_optional_u32(
    array: &dyn Array,
    row_idx: usize,
    column: &str,
) -> Result<Option<u32>, SpatialIoError> {
    let value = read_optional_u64(array, row_idx, column)?;
    match value {
        Some(v) => {
            let value_u32 = u32::try_from(v).map_err(|_| {
                SpatialIoError::UnsupportedFormat(format!(
                    "{} value out of range for u32: {}",
                    column, v
                ))
            })?;
            Ok(Some(value_u32))
        }
        None => Ok(None),
    }
}

fn read_optional_f32(
    array: &dyn Array,
    row_idx: usize,
    column: &str,
) -> Result<Option<f32>, SpatialIoError> {
    if array.is_null(row_idx) {
        return Ok(None);
    }

    let value = match array.data_type() {
        DataType::Float32 => read_primitive::<f32>(array, row_idx),
        DataType::Float64 => read_primitive::<f64>(array, row_idx).map(|v| v as f32),
        DataType::UInt8 => read_primitive::<u8>(array, row_idx).map(|v| v as f32),
        DataType::UInt16 => read_primitive::<u16>(array, row_idx).map(|v| v as f32),
        DataType::UInt32 => read_primitive::<u32>(array, row_idx).map(|v| v as f32),
        DataType::UInt64 => read_primitive::<u64>(array, row_idx).map(|v| v as f32),
        DataType::Int8 => read_primitive::<i8>(array, row_idx).map(|v| v as f32),
        DataType::Int16 => read_primitive::<i16>(array, row_idx).map(|v| v as f32),
        DataType::Int32 => read_primitive::<i32>(array, row_idx).map(|v| v as f32),
        DataType::Int64 => read_primitive::<i64>(array, row_idx).map(|v| v as f32),
        _ => None,
    }
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(format!(
            "unsupported parquet numeric column for {}: {:?}",
            column,
            array.data_type()
        ))
    })?;

    ensure_f32_finite_nonneg(value)?;
    Ok(Some(value))
}

fn read_primitive<T: arrow2::types::NativeType>(array: &dyn Array, row_idx: usize) -> Option<T> {
    array
        .as_any()
        .downcast_ref::<PrimitiveArray<T>>()
        .map(|a| a.value(row_idx))
}

fn to_u64_i128<T: Into<i128>>(value: T) -> Option<u64> {
    let v: i128 = value.into();
    if v < 0 || v > u64::MAX as i128 {
        None
    } else {
        Some(v as u64)
    }
}

fn compare_mapping_rows(a: &BarcodeMappingRow, b: &BarcodeMappingRow) -> Ordering {
    a.barcode
        .cmp(&b.barcode)
        .then_with(|| cmp_option_ord(&a.grid_row, &b.grid_row))
        .then_with(|| cmp_option_ord(&a.grid_col, &b.grid_col))
        .then_with(|| cmp_option_ord(&a.cell_id, &b.cell_id))
        .then_with(|| cmp_option_f32(&a.x, &b.x))
        .then_with(|| cmp_option_f32(&a.y, &b.y))
}

fn cmp_option_ord<T: Ord>(a: &Option<T>, b: &Option<T>) -> Ordering {
    match (a, b) {
        (Some(x), Some(y)) => x.cmp(y),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    }
}

fn cmp_option_f32(a: &Option<f32>, b: &Option<f32>) -> Ordering {
    match (a, b) {
        (Some(x), Some(y)) => x.total_cmp(y),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    }
}
