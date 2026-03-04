/// Main immutable dataset type.
pub mod dataset;

/// Feature-slice targeted loaders.
#[cfg(feature = "parquet")]
pub mod feature_slice;
