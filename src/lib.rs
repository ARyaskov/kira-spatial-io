#![forbid(unsafe_code)]
#![deny(clippy::all)]
#![deny(clippy::pedantic)]
#![deny(missing_docs)]
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
    clippy::uninlined_format_args,
    dead_code,
    unused_imports
)]

//! Deterministic IO library for spatial transcriptomics datasets.

/// Public dataset API entry points.
pub mod api;
/// Binary format primitives, reader, writer, and hashing.
pub mod binary;
/// Dataset loading configuration.
pub mod config;
/// Determinism helpers for sorting, float checks, and canonical JSON.
pub mod determinism;
/// Error model used across all IO operations.
pub mod error;
/// Input format loaders.
pub mod input;
/// Core in-memory data model types.
pub mod model;

pub use api::dataset::Dataset;
#[cfg(feature = "parquet")]
pub use api::feature_slice::{FeatureSliceGene, load_feature_slice_gene};
pub use config::LoadConfig;
pub use error::SpatialIoError;
pub use model::coord::CoordSystem;
pub use model::csr::BinsCsr;
pub use model::features::{FeatureRow, FeatureTable};
#[cfg(feature = "parquet")]
pub use model::mapping::{BarcodeMappingRow, BarcodeMappingTable};
pub use model::metadata::DatasetMetaCore;
pub use model::spatial_domain::SpatialDomain;
