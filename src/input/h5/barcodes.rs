use std::collections::HashSet;

use crate::error::SpatialIoError;
use crate::input::h5::strings::read_string_dataset_any;

pub(crate) fn load_barcodes(file: &hdf5::File) -> Result<Vec<String>, SpatialIoError> {
    let ds = file.dataset("/matrix/barcodes").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/barcodes dataset".to_string())
    })?;

    let barcodes = read_string_dataset_any(&ds, "/matrix/barcodes")?;

    let mut seen = HashSet::with_capacity(barcodes.len().min(1024));
    for barcode in &barcodes {
        if !seen.insert(barcode.clone()) {
            return Err(SpatialIoError::DuplicateBarcode(barcode.clone()));
        }
    }

    Ok(barcodes)
}
