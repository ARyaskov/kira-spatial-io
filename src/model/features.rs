/// Canonical feature row.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct FeatureRow {
    /// Canonical gene id assigned after sorting.
    pub gene_id: u32,
    /// Stable upstream identifier (Ensembl gene id, vendor accession, etc.).
    pub feature_id: String,
    /// Gene symbol/name used as the canonical sort key.
    pub gene_name: String,
    /// Feature type string from source dataset (e.g. `Gene Expression`).
    pub feature_type: String,
}

impl FeatureRow {
    /// Creates a canonical feature row.
    pub fn new(
        gene_id: u32,
        feature_id: String,
        gene_name: String,
        feature_type: String,
    ) -> Self {
        Self {
            gene_id,
            feature_id,
            gene_name,
            feature_type,
        }
    }
}

/// Canonical feature table ordered by `gene_name`.
#[derive(Clone, Debug)]
#[non_exhaustive]
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
