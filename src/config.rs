/// Compression policy for the binary writer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CompressionPolicy {
    /// No compression — sections are written as raw bytes.
    None,
    /// zstd compression at the given level (requires the `compression` feature).
    Zstd(i32),
}

/// Policy for handling duplicate `(bin, gene)` triplets encountered during ingest.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum DuplicatePolicy {
    /// Sum the duplicate values deterministically.
    Sum,
    /// Reject the input with a clear error.
    Error,
}

/// Loading options controlling validation strictness, memory budget and writer policy.
#[derive(Clone, Debug)]
pub struct LoadConfig {
    /// Maximum memory budget in megabytes for staged builders.
    pub memory_budget_mb: usize,
    /// Requested HD bin-level code, if applicable.
    pub bin_level: Option<u8>,
    /// Enables strict dimension/consistency checks.
    pub validate_strict: bool,
    /// When `true`, the reader recomputes and verifies the dataset hash on load.
    pub validate_hash: bool,
    /// Duplicate `(bin, gene)` handling on ingest.
    pub duplicate_policy: DuplicatePolicy,
    /// Compression policy for the binary writer.
    pub compression: CompressionPolicy,
}

impl Default for LoadConfig {
    fn default() -> Self {
        Self {
            memory_budget_mb: 8192,
            bin_level: None,
            validate_strict: true,
            validate_hash: true,
            duplicate_policy: DuplicatePolicy::Sum,
            compression: CompressionPolicy::None,
        }
    }
}
