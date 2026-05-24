/// Indptr storage variant — `u32` for compactness when feasible, `u64` otherwise.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Indptr {
    /// Compact 32-bit indptr (used when `nnz <= u32::MAX`).
    U32(Vec<u32>),
    /// Full 64-bit indptr.
    U64(Vec<u64>),
}

impl Indptr {
    /// Returns the entry at `i` as a `u64`.
    pub fn get(&self, i: usize) -> u64 {
        match self {
            Indptr::U32(v) => v[i] as u64,
            Indptr::U64(v) => v[i],
        }
    }

    /// Returns the number of entries (`n_bins + 1`).
    pub fn len(&self) -> usize {
        match self {
            Indptr::U32(v) => v.len(),
            Indptr::U64(v) => v.len(),
        }
    }

    /// Returns `true` when there are no entries.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the last entry as `u64`, or `None` if empty.
    pub fn last(&self) -> Option<u64> {
        match self {
            Indptr::U32(v) => v.last().copied().map(u64::from),
            Indptr::U64(v) => v.last().copied(),
        }
    }

    /// Returns the first entry as `u64`, or `None` if empty.
    pub fn first(&self) -> Option<u64> {
        match self {
            Indptr::U32(v) => v.first().copied().map(u64::from),
            Indptr::U64(v) => v.first().copied(),
        }
    }

    /// Returns `true` if the variant is the compact `u32` form.
    pub fn is_u32(&self) -> bool {
        matches!(self, Indptr::U32(_))
    }

    /// Materializes a `Vec<u64>` copy (allocates).
    pub fn to_u64_vec(&self) -> Vec<u64> {
        match self {
            Indptr::U32(v) => v.iter().copied().map(u64::from).collect(),
            Indptr::U64(v) => v.clone(),
        }
    }

    /// Picks the compact variant when `nnz` fits in `u32`, otherwise `U64`.
    pub fn from_u64(values: Vec<u64>) -> Self {
        if values.last().copied().unwrap_or(0) <= u64::from(u32::MAX) {
            Indptr::U32(values.into_iter().map(|v| v as u32).collect())
        } else {
            Indptr::U64(values)
        }
    }
}

/// Expression matrix in CSR format where rows are bins and columns are genes.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct BinsCsr {
    /// Row pointer array of length `n_bins + 1`.
    pub indptr: Indptr,
    /// Column indices per non-zero value.
    pub indices: Vec<u32>,
    /// Non-zero values aligned with [`indices`](Self::indices).
    pub data: Vec<f32>,
    /// Number of bins (rows).
    pub n_bins: u32,
    /// Number of genes (columns).
    pub n_genes: u32,
    /// Number of non-zero values.
    pub nnz: u64,
    /// Normalization marker (when `true`, negative values are permitted).
    pub normalized: bool,
}

impl BinsCsr {
    /// Constructs a CSR matrix container without additional validation.
    pub fn new(
        indptr: Indptr,
        indices: Vec<u32>,
        data: Vec<f32>,
        n_bins: u32,
        n_genes: u32,
        nnz: u64,
        normalized: bool,
    ) -> Self {
        Self {
            indptr,
            indices,
            data,
            n_bins,
            n_genes,
            nnz,
            normalized,
        }
    }
}
