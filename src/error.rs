/// Error type for deterministic IO, validation, and format handling.
#[derive(thiserror::Error, Debug)]
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
    /// Underlying IO failure.
    #[error(transparent)]
    Io(#[from] std::io::Error),
}
