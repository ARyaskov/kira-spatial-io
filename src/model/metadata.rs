use crate::model::coord::CoordSystem;

/// Fixed metadata core serialized in section `MetadataCore`.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct DatasetMetaCore {
    /// Dataset label.
    pub dataset_name: String,
    /// Source format identifier.
    pub source_format: String,
    /// Bin level code.
    pub bin_level: u8,
    /// Number of bins.
    pub n_bins: u32,
    /// Number of genes.
    pub n_genes: u32,
    /// Number of non-zero values.
    pub nnz: u64,
    /// Coordinate system.
    pub coord_system: CoordSystem,
    /// Normalization marker.
    pub normalized: bool,
    /// Canonical dataset hash (BLAKE3 leading 16 bytes). Populated by readers; ignored on write.
    pub dataset_hash: [u8; 16],
}

impl DatasetMetaCore {
    #[allow(clippy::too_many_arguments)]
    /// Creates fixed metadata core.
    pub fn new(
        dataset_name: String,
        source_format: String,
        bin_level: u8,
        n_bins: u32,
        n_genes: u32,
        nnz: u64,
        coord_system: CoordSystem,
        normalized: bool,
        dataset_hash: [u8; 16],
    ) -> Self {
        Self {
            dataset_name,
            source_format,
            bin_level,
            n_bins,
            n_genes,
            nnz,
            coord_system,
            normalized,
            dataset_hash,
        }
    }
}
