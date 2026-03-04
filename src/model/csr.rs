/// Expression matrix in CSR format where rows are bins and columns are genes.
#[derive(Clone, Debug)]
pub struct BinsCsr {
    /// Row pointer array of length `n_bins + 1`.
    pub indptr: Vec<u64>,
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
    /// Normalization marker.
    pub normalized: bool,
}

impl BinsCsr {
    /// Constructs a CSR matrix container without additional validation.
    pub fn new(
        indptr: Vec<u64>,
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
