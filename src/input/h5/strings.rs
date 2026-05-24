//! HDF5 string-dataset readers covering variable-length and fixed-width layouts.

use crate::error::SpatialIoError;

/// Decodes an HDF5 string dataset of unknown storage flavour into owned `String`s.
pub fn read_string_dataset_any(
    ds: &hdf5::Dataset,
    path: &str,
) -> Result<Vec<String>, SpatialIoError> {
    if let Ok(vals) = ds.read_raw::<hdf5::types::VarLenUnicode>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }
    if let Ok(vals) = ds.read_raw::<hdf5::types::VarLenAscii>() {
        return Ok(vals.into_iter().map(|v| v.to_string()).collect());
    }

    macro_rules! try_fixed {
        ($n:expr) => {
            if let Ok(vals) = ds.read_raw::<hdf5::types::FixedAscii<$n>>() {
                return Ok(vals
                    .into_iter()
                    .map(|v| v.as_str().trim_end_matches(char::from(0)).to_string())
                    .collect());
            }
            if let Ok(vals) = ds.read_raw::<hdf5::types::FixedUnicode<$n>>() {
                return Ok(vals
                    .into_iter()
                    .map(|v| v.as_str().trim_end_matches(char::from(0)).to_string())
                    .collect());
            }
            if let Ok(vals) = ds.read_raw::<[u8; $n]>() {
                return Ok(vals
                    .into_iter()
                    .map(|buf| {
                        let end = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
                        std::str::from_utf8(&buf[..end])
                            .unwrap_or("")
                            .trim()
                            .to_string()
                    })
                    .collect());
            }
        };
    }

    try_fixed!(16);
    try_fixed!(32);
    try_fixed!(48);
    try_fixed!(64);
    try_fixed!(96);
    try_fixed!(128);
    try_fixed!(192);
    try_fixed!(256);
    try_fixed!(384);
    try_fixed!(512);
    try_fixed!(1024);

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported string dataset type at {path}"
    )))
}

/// Tries to read an optional string dataset; returns `Ok(None)` when the path is missing.
pub fn read_optional_string_dataset(
    file: &hdf5::File,
    path: &str,
) -> Result<Option<Vec<String>>, SpatialIoError> {
    match file.dataset(path) {
        Ok(ds) => Ok(Some(read_string_dataset_any(&ds, path)?)),
        Err(_) => Ok(None),
    }
}
