use std::path::{Path, PathBuf};

use crate::error::SpatialIoError;

#[derive(Debug, Clone)]
pub(crate) struct TenxPaths {
    pub matrix_dir: PathBuf,
    pub matrix_mtx: PathBuf,
    pub barcodes_tsv: PathBuf,
    pub features_tsv: PathBuf,
    pub spatial_csv: PathBuf,
}

pub(crate) fn discover_tenx_paths(root: &Path) -> Result<TenxPaths, SpatialIoError> {
    let filtered = root.join("filtered_feature_bc_matrix");
    let matrix_dir = if filtered.is_dir() {
        filtered
    } else {
        root.to_path_buf()
    };

    let matrix_mtx = choose_existing(&[
        matrix_dir.join("matrix.mtx"),
        matrix_dir.join("matrix.mtx.gz"),
    ])
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing matrix.mtx(.gz) in {}",
            matrix_dir.display()
        ))
    })?;

    let barcodes_tsv = choose_existing(&[
        matrix_dir.join("barcodes.tsv"),
        matrix_dir.join("barcodes.tsv.gz"),
    ])
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing barcodes.tsv(.gz) in {}",
            matrix_dir.display()
        ))
    })?;

    let features_tsv = choose_existing(&[
        matrix_dir.join("features.tsv"),
        matrix_dir.join("features.tsv.gz"),
    ])
    .ok_or_else(|| {
        if choose_existing(&[
            matrix_dir.join("genes.tsv"),
            matrix_dir.join("genes.tsv.gz"),
        ])
        .is_some()
        {
            SpatialIoError::UnsupportedFormat(
                "genes.tsv legacy layout is not supported in Stage 1".to_string(),
            )
        } else {
            SpatialIoError::UnsupportedFormat(format!(
                "missing features.tsv(.gz) in {}",
                matrix_dir.display()
            ))
        }
    })?;

    let spatial_dir = root.join("spatial");
    if !spatial_dir.is_dir() {
        return Err(SpatialIoError::UnsupportedFormat(
            "missing spatial/".to_string(),
        ));
    }

    let spatial_csv = choose_existing(&[
        spatial_dir.join("tissue_positions.csv"),
        spatial_dir.join("tissue_positions_list.csv"),
    ])
    .ok_or_else(|| {
        SpatialIoError::UnsupportedFormat(format!(
            "missing tissue_positions.csv or tissue_positions_list.csv in {}",
            spatial_dir.display()
        ))
    })?;

    Ok(TenxPaths {
        matrix_dir,
        matrix_mtx,
        barcodes_tsv,
        features_tsv,
        spatial_csv,
    })
}

fn choose_existing(paths: &[PathBuf]) -> Option<PathBuf> {
    paths.iter().find(|p| p.is_file()).cloned()
}
