use std::path::{Path, PathBuf};

/// Error type for deterministic IO, validation, and format handling.
#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum SpatialIoError {
    /// Dimensional or shape mismatch between linked arrays/sections.
    #[error("dimension mismatch: {0}")]
    DimensionMismatch(String),
    /// Duplicate gene key detected in canonical feature table.
    #[error("duplicate gene: {0}")]
    DuplicateGene(String),
    /// Duplicate barcode key detected where uniqueness is required.
    #[error("duplicate barcode: {0}")]
    DuplicateBarcode(String),
    /// Invalid CSR payload or invariant violation.
    #[error("invalid CSR: {0}")]
    InvalidCsr(String),
    /// Invalid floating-point value or numeric policy violation.
    #[error("invalid float: {0}")]
    InvalidFloat(String),
    /// Requested operation exceeds configured memory budget.
    #[error("memory limit exceeded: {0}")]
    MemoryLimitExceeded(String),
    /// Unsupported input layout or binary format condition.
    #[error("unsupported format: {0}")]
    UnsupportedFormat(String),
    /// Per-section CRC32C check failed (file corruption).
    #[error("CRC mismatch in section {section_id}: expected {expected:#010x}, got {actual:#010x}")]
    CrcMismatch {
        /// Section identifier.
        section_id: u16,
        /// Expected CRC32C from the section table.
        expected: u32,
        /// Actual CRC32C computed from on-disk bytes.
        actual: u32,
    },
    /// Dataset hash mismatch (file corruption or tampering).
    #[error("dataset hash mismatch")]
    HashMismatch,
    /// IO failure with file path context.
    #[error("io error on {path}: {source}")]
    IoAt {
        /// Path involved in the failed operation.
        path: PathBuf,
        /// Underlying IO error.
        #[source]
        source: std::io::Error,
    },
    /// Underlying IO failure without path context.
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

impl SpatialIoError {
    /// Wraps an `io::Error` with the path that caused it.
    pub fn io_at(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Self::IoAt {
            path: path.as_ref().to_path_buf(),
            source,
        }
    }
}

/// Extension trait to ergonomically attach a path to an `io::Result`.
pub trait IoPathExt<T> {
    /// Maps the error variant to [`SpatialIoError::IoAt`].
    fn io_path(self, path: impl AsRef<Path>) -> Result<T, SpatialIoError>;
}

impl<T> IoPathExt<T> for std::io::Result<T> {
    fn io_path(self, path: impl AsRef<Path>) -> Result<T, SpatialIoError> {
        self.map_err(|e| SpatialIoError::io_at(path, e))
    }
}
