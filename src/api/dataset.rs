use std::path::Path;

use serde_json::Value;

use crate::binary::reader::{read_kira_bin, read_kira_bin_with};
use crate::binary::writer::{write_kira_bin, write_kira_bin_with_compression};
use crate::config::{CompressionPolicy, LoadConfig};
use crate::error::SpatialIoError;
use crate::input::tenx::load_10x_mtx;
#[cfg(feature = "hdf5")]
use crate::input::h5::load_h5_dataset;
#[cfg(feature = "parquet")]
use crate::model::mapping::BarcodeMappingTable;
use crate::model::{
    csr::BinsCsr, features::FeatureTable, metadata::DatasetMetaCore, spatial_domain::SpatialDomain,
};

/// Immutable loaded dataset with canonical ordering and validated invariants.
#[derive(Debug)]
pub struct Dataset {
    spatial_domain: SpatialDomain,
    expression_csr: BinsCsr,
    features: FeatureTable,
    metadata_core: DatasetMetaCore,
    metadata_json: Value,
    #[cfg(feature = "parquet")]
    barcode_mapping: Option<BarcodeMappingTable>,
}

impl Dataset {
    /// Opens a 10x-style directory (including Visium HD binned layout).
    pub fn open_10x<P: AsRef<Path>>(path: P, cfg: LoadConfig) -> Result<Self, SpatialIoError> {
        load_10x_mtx(path, cfg)
    }

    /// Opens a 10x H5 matrix layout with required spatial coordinates.
    #[cfg(feature = "hdf5")]
    pub fn open_h5<P: AsRef<Path>>(path: P, cfg: LoadConfig) -> Result<Self, SpatialIoError> {
        load_h5_dataset(path, cfg)
    }

    /// Loads a previously exported `.kira-spatial.bin` file with default config.
    pub fn from_kira_bin<P: AsRef<Path>>(path: P) -> Result<Self, SpatialIoError> {
        read_kira_bin(path)
    }

    /// Loads a `.kira-spatial.bin` file with explicit config (hash skip, memory cap, ...).
    pub fn from_kira_bin_with<P: AsRef<Path>>(
        path: P,
        cfg: &LoadConfig,
    ) -> Result<Self, SpatialIoError> {
        read_kira_bin_with(path, cfg)
    }

    /// Exports the dataset as deterministic `.kira-spatial.bin`.
    pub fn export_kira_bin<P: AsRef<Path>>(&self, path: P) -> Result<(), SpatialIoError> {
        write_kira_bin(path, self)
    }

    /// Exports the dataset with an explicit compression policy.
    pub fn export_kira_bin_with<P: AsRef<Path>>(
        &self,
        path: P,
        compression: CompressionPolicy,
    ) -> Result<(), SpatialIoError> {
        write_kira_bin_with_compression(path, self, compression)
    }

    /// Returns canonical spatial domain.
    pub fn spatial_domain(&self) -> &SpatialDomain {
        &self.spatial_domain
    }

    /// Returns expression matrix in CSR representation.
    pub fn expression_csr(&self) -> &BinsCsr {
        &self.expression_csr
    }

    /// Returns canonical feature table.
    pub fn features(&self) -> &FeatureTable {
        &self.features
    }

    /// Returns fixed metadata core.
    pub fn metadata_core(&self) -> &DatasetMetaCore {
        &self.metadata_core
    }

    /// Returns canonical metadata JSON.
    pub fn metadata_json(&self) -> &Value {
        &self.metadata_json
    }

    #[cfg(feature = "parquet")]
    /// Returns optional parquet barcode mapping table.
    pub fn barcode_mapping(&self) -> Option<&BarcodeMappingTable> {
        self.barcode_mapping.as_ref()
    }

    /// Builds a [`Dataset`] from already-validated components.
    pub fn from_parts(
        spatial_domain: SpatialDomain,
        expression_csr: BinsCsr,
        features: FeatureTable,
        metadata_core: DatasetMetaCore,
        metadata_json: Value,
    ) -> Self {
        Self {
            spatial_domain,
            expression_csr,
            features,
            metadata_core,
            metadata_json,
            #[cfg(feature = "parquet")]
            barcode_mapping: None,
        }
    }

    #[cfg(feature = "parquet")]
    /// Attaches an optional barcode mapping table to the dataset.
    pub fn with_barcode_mapping(
        mut self,
        barcode_mapping: Option<BarcodeMappingTable>,
    ) -> Self {
        self.barcode_mapping = barcode_mapping;
        self
    }
}
