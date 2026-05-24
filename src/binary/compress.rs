//! Optional zstd compression for section payloads.

use crate::error::SpatialIoError;

#[cfg(feature = "compression")]
pub fn zstd_compress(input: &[u8], level: i32) -> Result<Vec<u8>, SpatialIoError> {
    zstd::encode_all(input, level).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!("zstd encode error: {e}"))
    })
}

#[cfg(feature = "compression")]
pub fn zstd_decompress(input: &[u8]) -> Result<Vec<u8>, SpatialIoError> {
    zstd::decode_all(input).map_err(|e| {
        SpatialIoError::UnsupportedFormat(format!("zstd decode error: {e}"))
    })
}

#[cfg(not(feature = "compression"))]
pub fn zstd_compress(_input: &[u8], _level: i32) -> Result<Vec<u8>, SpatialIoError> {
    Err(SpatialIoError::UnsupportedFormat(
        "compression feature not enabled".to_string(),
    ))
}

#[cfg(not(feature = "compression"))]
pub fn zstd_decompress(_input: &[u8]) -> Result<Vec<u8>, SpatialIoError> {
    Err(SpatialIoError::UnsupportedFormat(
        "compression feature not enabled: cannot decode a compressed section".to_string(),
    ))
}
