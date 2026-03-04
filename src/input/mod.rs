/// H5 input adapters.
pub mod h5;
#[cfg(feature = "parquet")]
/// Parquet input adapters.
pub mod parquet;
/// 10x directory/MTX input adapters.
pub mod tenx;
/// Shared input utilities.
pub mod util;
