use std::path::{Path, PathBuf};

use crate::error::SpatialIoError;

#[derive(Debug, Clone)]
pub(crate) enum SpatialInput {
    Csv(PathBuf),
    #[cfg(feature = "parquet")]
    Parquet(PathBuf),
}

#[derive(Debug, Clone)]
pub(crate) struct H5Paths {
    pub root_dir: PathBuf,
    pub h5_path: PathBuf,
    pub spatial_input: SpatialInput,
}

pub(crate) fn discover_h5_paths(path: &Path) -> Result<H5Paths, SpatialIoError> {
    let (root_dir, h5_path) = if path.is_dir() {
        let candidate = path.join("feature_slice.h5");
        if !candidate.is_file() {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "missing feature_slice.h5 in {}",
                path.display()
            )));
        }
        (path.to_path_buf(), candidate)
    } else {
        let ext_ok = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("h5"));
        if !ext_ok {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "expected .h5 file or directory with feature_slice.h5: {}",
                path.display()
            )));
        }
        let root = path.parent().ok_or_else(|| {
            SpatialIoError::UnsupportedFormat("cannot determine parent dir for h5 file".to_string())
        })?;
        (root.to_path_buf(), path.to_path_buf())
    };

    let spatial_csv = {
        let spatial_dir = root_dir.join("spatial");
        if spatial_dir.is_dir() {
            [
                spatial_dir.join("tissue_positions.csv"),
                spatial_dir.join("tissue_positions_list.csv"),
            ]
            .into_iter()
            .find(|p| p.is_file())
        } else {
            None
        }
    };

    if let Some(spatial_csv) = spatial_csv {
        return Ok(H5Paths {
            root_dir,
            h5_path,
            spatial_input: SpatialInput::Csv(spatial_csv),
        });
    }

    #[cfg(feature = "parquet")]
    {
        let spatial_parquet = [
            root_dir.join("barcode_mappings.parquet"),
            root_dir.join("barcode_mapping.parquet"),
            root_dir.join("spatial").join("barcode_mappings.parquet"),
            root_dir.join("spatial").join("barcode_mapping.parquet"),
        ]
        .into_iter()
        .find(|p| p.is_file());

        if let Some(spatial_parquet) = spatial_parquet {
            return Ok(H5Paths {
                root_dir,
                h5_path,
                spatial_input: SpatialInput::Parquet(spatial_parquet),
            });
        }
    }

    Err(SpatialIoError::UnsupportedFormat(
        "missing spatial coordinates: provide spatial/tissue_positions.csv (or tissue_positions_list.csv), or barcode_mappings.parquet when parquet feature is enabled"
            .to_string(),
    ))
}
