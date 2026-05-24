//! Deterministic IO library for spatial transcriptomics datasets.

#![deny(unsafe_code)]
#![deny(clippy::all)]
#![allow(
    clippy::bool_to_int_with_if,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_lossless,
    clippy::cast_sign_loss,
    clippy::doc_markdown,
    clippy::module_name_repetitions,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::needless_pass_by_value,
    clippy::redundant_closure_for_method_calls,
    clippy::ref_option,
    clippy::similar_names,
    clippy::struct_excessive_bools,
    clippy::too_many_lines,
    clippy::trivially_copy_pass_by_ref,
    clippy::uninlined_format_args
)]

#[cfg(target_endian = "big")]
compile_error!("kira-spatial-io requires a little-endian target");

pub mod api;
pub mod binary;
pub mod config;
pub mod determinism;
pub mod error;
pub mod input;
pub mod model;

pub use api::dataset::Dataset;
#[cfg(all(feature = "parquet", feature = "hdf5"))]
pub use api::feature_slice::{FeatureSliceGene, load_feature_slice_gene};
pub use config::{CompressionPolicy, DuplicatePolicy, LoadConfig};
pub use error::{IoPathExt, SpatialIoError};
pub use model::coord::CoordSystem;
pub use model::csr::{BinsCsr, Indptr};
pub use model::features::{FeatureRow, FeatureTable};
#[cfg(feature = "parquet")]
pub use model::mapping::{BarcodeMappingRow, BarcodeMappingTable};
pub use model::metadata::DatasetMetaCore;
pub use model::spatial_domain::SpatialDomain;
