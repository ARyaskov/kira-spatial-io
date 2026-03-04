use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use hdf5::File;

use crate::config::LoadConfig;
use crate::determinism::sort::sort_bins;
use crate::error::SpatialIoError;
use crate::input::h5::discover;
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
        #[cfg(feature = "parquet")]
        discover::SpatialInput::Parquet(path) => path.clone(),
    };

    let h5 = File::open(&paths.h5_path).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!(
            "failed to open h5 file {}: {e}",
            paths.h5_path.display()
        ))
    })?;

    let feature_index = match resolve_feature_index_from_h5(&h5, gene_name) {
        Ok(idx) => idx,
        Err(_) => resolve_feature_index_with_h5dump(&paths.h5_path, gene_name)?,
    };

    let slice_base = format!("/feature_slices/{feature_index}");
    let row_ds = h5.dataset(&(slice_base.clone() + "/row")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {}/row dataset for gene {}",
            slice_base, gene_name
        ))
    })?;
    let col_ds = h5.dataset(&(slice_base.clone() + "/col")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {}/col dataset for gene {}",
            slice_base, gene_name
        ))
    })?;
    let data_ds = h5.dataset(&(slice_base.clone() + "/data")).map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing {}/data dataset for gene {}",
            slice_base, gene_name
        ))
    })?;

    let rows =
        read_u32_dataset_with_fallback(&row_ds, &paths.h5_path, &(slice_base.clone() + "/row"))?;
    let cols =
        read_u32_dataset_with_fallback(&col_ds, &paths.h5_path, &(slice_base.clone() + "/col"))?;
    let counts =
        read_u32_dataset_with_fallback(&data_ds, &paths.h5_path, &(slice_base.clone() + "/data"))?;

    if rows.len() != cols.len() || rows.len() != counts.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "feature slice arrays length mismatch for {}: row={}, col={}, data={}",
            gene_name,
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
            "feature slice entries without coordinate match in parquet mapping: {}",
            missing_coords
        )));
    }

    Ok(FeatureSliceGene {
        spatial_domain,
        values,
    })
}

fn resolve_feature_index_from_h5(
    file: &hdf5::File,
    gene_name: &str,
) -> Result<usize, SpatialIoError> {
    let names_ds = file.dataset("/features/name").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /features/name dataset".to_string())
    })?;
    let feature_names = read_string_dataset_any(&names_ds, "/features/name")?;
    feature_names
        .iter()
        .position(|name| name == gene_name)
        .ok_or_else(|| {
            SpatialIoError::UnsupportedFormat(format!(
                "gene not found in /features/name: {}",
                gene_name
            ))
        })
}

fn resolve_feature_index_with_h5dump(
    h5_path: &Path,
    gene_name: &str,
) -> Result<usize, SpatialIoError> {
    let out = Command::new("h5dump")
        .arg("-d")
        .arg("/features/name")
        .arg(h5_path)
        .output()
        .map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!(
                "failed to run h5dump for /features/name: {e}"
            ))
        })?;

    if !out.status.success() {
        return Err(SpatialIoError::UnsupportedFormat(
            "h5dump failed reading /features/name".to_string(),
        ));
    }

    let stdout = String::from_utf8_lossy(&out.stdout);
    let marker_with_null = format!("\"{}\\000", gene_name);
    let marker_exact = format!("\"{}\"", gene_name);
    for line in stdout.lines() {
        if !(line.contains(&marker_with_null) || line.contains(&marker_exact)) {
            continue;
        }
        if let Some((idx, _name)) = parse_h5dump_feature_line(line) {
            return Ok(idx);
        }
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "gene not found in /features/name: {}",
        gene_name
    )))
}

fn parse_h5dump_feature_line(line: &str) -> Option<(usize, String)> {
    let start = line.find('(')?;
    let end = line[start + 1..].find(')')? + start + 1;
    let idx = line[start + 1..end].trim().parse::<usize>().ok()?;

    let quote_start = line[end..].find('"')? + end + 1;
    let rest = &line[quote_start..];
    let quote_end = rest.find('"')?;
    let mut value = rest[..quote_end].to_string();

    if let Some(null_pos) = value.find("\\000") {
        value.truncate(null_pos);
    }

    Some((idx, value))
}

fn read_u32_dataset_with_fallback(
    ds: &hdf5::Dataset,
    h5_path: &Path,
    dataset_path: &str,
) -> Result<Vec<u32>, SpatialIoError> {
    if let Ok(vals) = ds.read_raw::<u32>() {
        return Ok(vals);
    }

    let out = Command::new("h5dump")
        .arg("-d")
        .arg(dataset_path)
        .arg(h5_path)
        .output()
        .map_err(|e| {
            SpatialIoError::UnsupportedFormat(format!(
                "failed to run h5dump for {}: {e}",
                dataset_path
            ))
        })?;

    if !out.status.success() {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "h5dump failed reading {}",
            dataset_path
        )));
    }

    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut values = Vec::<u32>::new();
    for line in stdout.lines() {
        let Some(colon) = line.find(':') else {
            continue;
        };
        let payload = &line[colon + 1..];
        for token in payload.split(|ch: char| !ch.is_ascii_digit()) {
            if token.is_empty() {
                continue;
            }
            values.push(token.parse::<u32>().map_err(|_| {
                SpatialIoError::UnsupportedFormat(format!(
                    "invalid numeric token in h5dump output for {}",
                    dataset_path
                ))
            })?);
        }
    }

    if values.is_empty() {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "no values parsed from h5dump output for {}",
            dataset_path
        )));
    }

    Ok(values)
}

fn read_string_dataset_any(ds: &hdf5::Dataset, path: &str) -> Result<Vec<String>, SpatialIoError> {
    if let Ok(vals) = ds.read_raw::<hdf5::types::VarLenUnicode>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }
    if let Ok(vals) = ds.read_raw::<hdf5::types::VarLenAscii>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }

    macro_rules! try_fixed_ascii {
        ($n:expr) => {
            if let Ok(vals) = ds.read_raw::<hdf5::types::FixedAscii<$n>>() {
                return Ok(vals
                    .into_iter()
                    .map(|v| v.as_str().trim_matches(char::from(0)).to_string())
                    .collect());
            }
            if let Ok(vals) = ds.read_raw::<hdf5::types::FixedUnicode<$n>>() {
                return Ok(vals
                    .into_iter()
                    .map(|v| v.as_str().trim_matches(char::from(0)).to_string())
                    .collect());
            }
            if let Ok(vals) = ds.read_raw::<[u8; $n]>() {
                return Ok(vals
                    .into_iter()
                    .map(|buf| {
                        let end = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
                        String::from_utf8_lossy(&buf[..end]).trim().to_string()
                    })
                    .collect());
            }
        };
    }
    try_fixed_ascii!(32);
    try_fixed_ascii!(64);
    try_fixed_ascii!(96);
    try_fixed_ascii!(128);
    try_fixed_ascii!(192);
    try_fixed_ascii!(256);
    try_fixed_ascii!(384);
    try_fixed_ascii!(512);
    try_fixed_ascii!(1024);

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported string dataset type at {}",
        path
    )))
}
