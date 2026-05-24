use std::path::Path;

use hdf5::File;
use serde_json::json;

use crate::api::dataset::Dataset;
use crate::config::LoadConfig;
use crate::determinism::{json::canonicalize_json, sort::sort_bins};
use crate::error::SpatialIoError;
use crate::input::h5::{barcodes, discover, features, matrix};
use crate::input::tenx::spatial;
use crate::model::metadata::DatasetMetaCore;

pub fn load_10x_h5<P: AsRef<Path>>(path: P, cfg: LoadConfig) -> Result<Dataset, SpatialIoError> {
    let paths = discover::discover_h5_paths(path.as_ref())?;
    let h5 = File::open(&paths.h5_path).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!(
            "failed to open h5 file {}: {e}",
            paths.h5_path.display()
        ))
    })?;

    let barcodes = barcodes::load_barcodes(&h5)?;
    let feature_build = features::load_features(&h5)?;

    let mut spatial_domain = match &paths.spatial_input {
        discover::SpatialInput::Csv(csv_path) => {
            spatial::load_spatial_domain(csv_path, &barcodes, &cfg)?
        }
        #[cfg(feature = "parquet")]
        discover::SpatialInput::Parquet(parquet_path) => {
            crate::input::parquet::spatial_mapping::load_spatial_domain_from_mapping_parquet(
                parquet_path,
                &barcodes,
                &cfg,
            )?
        }
    };
    let sort_result = sort_bins(&mut spatial_domain)?;

    let csr = matrix::build_csr_from_h5(
        &h5,
        &cfg,
        &sort_result.old_to_new,
        &feature_build.old_to_new,
    )?;

    let n_bins = barcodes.len() as u32;
    let n_genes = feature_build.table.rows.len() as u32;

    let metadata_core = DatasetMetaCore {
        dataset_name: paths
            .root_dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("dataset")
            .to_string(),
        source_format: "10x-h5".to_string(),
        bin_level: spatial_domain.bin_level,
        n_bins,
        n_genes,
        nnz: csr.nnz,
        coord_system: spatial_domain.coord_system,
        normalized: false,
        dataset_hash: [0_u8; 16],
    };

    let spatial_source = match &paths.spatial_input {
        discover::SpatialInput::Csv(path) => path.display().to_string(),
        #[cfg(feature = "parquet")]
        discover::SpatialInput::Parquet(path) => path.display().to_string(),
    };

    let metadata_json = canonicalize_json(&json!({
        "source": {
            "path": paths.root_dir.display().to_string(),
            "h5_file": paths.h5_path.display().to_string(),
            "spatial_source": spatial_source,
        },
        "tenx": {
            "has_spatial": true,
            "duplicate_policy": "sum-per-bin-gene"
        }
    }));

    Ok(Dataset::from_parts(
        spatial_domain,
        csr,
        feature_build.table,
        metadata_core,
        metadata_json,
    ))
}
