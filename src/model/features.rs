/// Canonical feature row.
#[derive(Clone, Debug)]
pub struct FeatureRow {
    /// Canonical gene id assigned after sorting.
    pub gene_id: u32,
    /// Gene symbol/name used as canonical sort key.
    pub gene_name: String,
    /// Feature type string from source dataset.
    pub feature_type: String,
}

impl FeatureRow {
    /// Creates a canonical feature row.
    pub fn new(gene_id: u32, gene_name: String, feature_type: String) -> Self {
        Self {
            gene_id,
            gene_name,
            feature_type,
        }
    }
}

/// Canonical feature table ordered by `gene_name`.
#[derive(Clone, Debug)]
pub struct FeatureTable {
    /// Feature rows.
    pub rows: Vec<FeatureRow>,
}

impl FeatureTable {
    /// Creates a feature table from rows.
    pub fn new(rows: Vec<FeatureRow>) -> Self {
        Self { rows }
    }

    /// Returns the number of rows.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Returns `true` when no rows are present.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
}
