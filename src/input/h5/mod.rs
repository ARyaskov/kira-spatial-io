use std::path::Path;

use crate::api::dataset::Dataset;
use crate::config::LoadConfig;
use crate::error::SpatialIoError;

pub(crate) mod barcodes;
pub(crate) mod discover;
pub(crate) mod features;
pub(crate) mod loader;
pub(crate) mod matrix;

pub(crate) fn load_h5_dataset<P: AsRef<Path>>(
    path: P,
    cfg: LoadConfig,
) -> Result<Dataset, SpatialIoError> {
    loader::load_10x_h5(path, cfg)
}
