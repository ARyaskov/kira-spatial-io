/// Coordinate-system enum.
pub mod coord;
/// CSR expression model.
pub mod csr;
/// Feature table model.
pub mod features;
#[cfg(feature = "parquet")]
/// Optional barcode mapping model (parquet feature).
pub mod mapping;
/// Fixed metadata model.
pub mod metadata;
/// Spatial domain model.
pub mod spatial_domain;
