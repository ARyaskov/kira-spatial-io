use std::path::Path;

use serde_json::json;

use crate::api::dataset::Dataset;
use crate::config::LoadConfig;
use crate::determinism::{json::canonicalize_json, sort::sort_bins};
use crate::error::SpatialIoError;
use crate::model::metadata::DatasetMetaCore;
#[cfg(feature = "parquet")]
use crate::{
    input::parquet::barcode_mapping::load_barcode_mapping_parquet,
    model::mapping::BarcodeMappingTable,
};

pub(crate) mod barcodes;
pub(crate) mod discover;
pub(crate) mod features;
pub(crate) mod hd_binned;
pub(crate) mod mtx;
pub(crate) mod spatial;

pub(crate) fn load_10x_mtx<P: AsRef<Path>>(
    path: P,
    cfg: LoadConfig,
) -> Result<Dataset, SpatialIoError> {
    let root = path.as_ref();

    let (dataset_root, effective_bin_level, source_layout, bin_size_um) =
        if let Some(layout) = hd_binned::discover_hd_binned(root)? {
            let selected = hd_binned::select_bin_folder(&layout, &cfg)?;
            (
                selected.path,
                Some(selected.level_code),
                "visium-hd-binned".to_string(),
                Some(selected.um),
            )
        } else {
            (
                root.to_path_buf(),
                cfg.bin_level,
                "tenx-mtx".to_string(),
                None,
            )
        };

    let paths = discover::discover_tenx_paths(&dataset_root)?;

    let barcodes = barcodes::load_barcodes(&paths.barcodes_tsv)?;
    let feature_rows_raw = features::load_feature_rows(&paths.features_tsv)?;
    let feature_build = features::build_feature_table(feature_rows_raw)?;

    let effective_cfg = LoadConfig {
        memory_budget_mb: cfg.memory_budget_mb,
        bin_level: effective_bin_level,
        validate_strict: cfg.validate_strict,
    };

    let mut spatial_domain =
        spatial::load_spatial_domain(&paths.spatial_csv, &barcodes, &effective_cfg)?;
    let sort_result = sort_bins(&mut spatial_domain)?;

    let csr = mtx::build_csr_from_mtx_path(
        &paths.matrix_mtx,
        &effective_cfg,
        &sort_result.old_to_new,
        &feature_build.old_to_new,
    )?;

    let n_bins = barcodes.len() as u32;
    let n_genes = feature_build.table.rows.len() as u32;
    let metadata_core = DatasetMetaCore {
        dataset_name: root
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("dataset")
            .to_string(),
        source_format: if bin_size_um.is_some() {
            "10x-hd-binned-mtx".to_string()
        } else {
            "10x-mtx".to_string()
        },
        bin_level: spatial_domain.bin_level,
        n_bins,
        n_genes,
        nnz: csr.nnz,
        coord_system: spatial_domain.coord_system,
        normalized: false,
        dataset_hash: [0_u8; 16],
    };

    let metadata_json = json!({
        "source": {
            "path": root.display().to_string(),
            "dataset_root": dataset_root.display().to_string(),
            "matrix_dir": paths.matrix_dir.display().to_string(),
            "matrix_mtx": paths.matrix_mtx.display().to_string(),
            "layout": source_layout
        },
        "tenx": {
            "has_spatial": true,
            "duplicate_policy": "sum-per-bin-gene",
            "hd": bin_size_um.is_some(),
            "bin_size_um": bin_size_um.unwrap_or(0),
            "bin_level_code": spatial_domain.bin_level
        }
    });

    #[cfg(feature = "parquet")]
    let mut metadata_json = metadata_json;

    #[cfg(feature = "parquet")]
    let barcode_mapping: Option<BarcodeMappingTable> = {
        if bin_size_um.is_some() {
            if let Some(mapping_path) = hd_binned::discover_barcode_mapping_parquet(root) {
                let (table, summary) = load_barcode_mapping_parquet(&mapping_path, &effective_cfg)?;
                if let Some(tenx_obj) = metadata_json
                    .get_mut("tenx")
                    .and_then(|v| v.as_object_mut())
                {
                    tenx_obj.insert("barcode_mapping_present".to_string(), json!(true));
                    tenx_obj.insert(
                        "barcode_mapping_path".to_string(),
                        json!(mapping_path.display().to_string()),
                    );
                    tenx_obj.insert("barcode_mapping_rows".to_string(), json!(summary.row_count));
                    tenx_obj.insert(
                        "barcode_mapping_has_cell_id".to_string(),
                        json!(summary.has_cell_id),
                    );
                    tenx_obj.insert(
                        "barcode_mapping_has_grid".to_string(),
                        json!(summary.has_grid),
                    );
                    tenx_obj.insert("barcode_mapping_has_xy".to_string(), json!(summary.has_xy));
                    tenx_obj.insert(
                        "mapping_has_duplicates".to_string(),
                        json!(summary.mapping_has_duplicates),
                    );
                }
                Some(table)
            } else {
                if let Some(tenx_obj) = metadata_json
                    .get_mut("tenx")
                    .and_then(|v| v.as_object_mut())
                {
                    tenx_obj.insert("barcode_mapping_present".to_string(), json!(false));
                }
                None
            }
        } else {
            None
        }
    };

    let metadata_json = canonicalize_json(&metadata_json);

    let dataset = Dataset::from_parts(
        spatial_domain,
        csr,
        feature_build.table,
        metadata_core,
        metadata_json,
    );

    #[cfg(feature = "parquet")]
    let dataset = dataset.with_barcode_mapping(barcode_mapping);

    Ok(dataset)
}
