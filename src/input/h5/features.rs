use crate::error::SpatialIoError;
use crate::input::h5::strings::{read_optional_string_dataset, read_string_dataset_any};
use crate::input::tenx::features::{FeatureBuildResult, FeatureRowRaw, build_feature_table};

pub(crate) fn load_features(file: &hdf5::File) -> Result<FeatureBuildResult, SpatialIoError> {
    let name_ds = file.dataset("/matrix/features/name").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/features/name dataset".to_string())
    })?;
    let feature_type_ds = file.dataset("/matrix/features/feature_type").map_err(|_| {
        SpatialIoError::UnsupportedFormat(
            "missing /matrix/features/feature_type dataset".to_string(),
        )
    })?;

    let names = read_string_dataset_any(&name_ds, "/matrix/features/name")?;
    let feature_types = read_string_dataset_any(&feature_type_ds, "/matrix/features/feature_type")?;

    if names.len() != feature_types.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "feature name and feature_type lengths differ: {} vs {}",
            names.len(),
            feature_types.len()
        )));
    }

    let feature_ids = read_optional_string_dataset(file, "/matrix/features/id")?
        .unwrap_or_else(|| vec![String::new(); names.len()]);
    if feature_ids.len() != names.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "feature id length {} differs from name length {}",
            feature_ids.len(),
            names.len()
        )));
    }

    let rows: Vec<FeatureRowRaw> = names
        .into_iter()
        .zip(feature_types)
        .zip(feature_ids)
        .map(|((gene_name, feature_type), feature_id)| FeatureRowRaw {
            feature_id,
            gene_name,
            feature_type,
        })
        .collect();

    build_feature_table(rows)
}
