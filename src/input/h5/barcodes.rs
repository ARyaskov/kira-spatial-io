use std::collections::HashSet;

use hdf5::types::{VarLenAscii, VarLenUnicode};

use crate::error::SpatialIoError;

pub(crate) fn load_barcodes(file: &hdf5::File) -> Result<Vec<String>, SpatialIoError> {
    let ds = file.dataset("/matrix/barcodes").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/barcodes dataset".to_string())
    })?;

    let barcodes = read_string_dataset(&ds, "/matrix/barcodes")?;

    let mut seen = HashSet::with_capacity(barcodes.len());
    for barcode in &barcodes {
        if !seen.insert(barcode.clone()) {
            return Err(SpatialIoError::DuplicateBarcode(barcode.clone()));
        }
    }

    Ok(barcodes)
}

pub(crate) fn read_string_dataset(
    ds: &hdf5::Dataset,
    path: &str,
) -> Result<Vec<String>, SpatialIoError> {
    if let Ok(vals) = ds.read_raw::<VarLenUnicode>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }
    if let Ok(vals) = ds.read_raw::<VarLenAscii>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported string dataset type at {}",
        path
    )))
}
