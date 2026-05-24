use std::path::Path;
use std::io::BufRead;

use crate::error::{IoPathExt, SpatialIoError};
use crate::input::util::open_text_maybe_gz;
use crate::model::features::{FeatureRow, FeatureTable};

/// Raw feature row parsed from `features.tsv`.
#[derive(Debug, Clone)]
pub(crate) struct FeatureRowRaw {
    /// Stable upstream feature identifier (column 0; typically Ensembl gene id).
    pub feature_id: String,
    /// Raw gene name (column 1).
    pub gene_name: String,
    /// Raw feature type string (column 2; defaults to "Gene Expression" if absent).
    pub feature_type: String,
}

/// Canonicalized feature-table build result.
#[derive(Debug, Clone)]
pub(crate) struct FeatureBuildResult {
    /// Canonical feature table.
    pub table: FeatureTable,
    /// Mapping from file-order feature index to canonical `gene_id`.
    pub old_to_new: Vec<u32>,
}

/// Loads raw feature rows from `features.tsv(.gz)`.
pub(crate) fn load_feature_rows(path: &Path) -> Result<Vec<FeatureRowRaw>, SpatialIoError> {
    let reader = open_text_maybe_gz(path)?;
    let mut out = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.io_path(path)?;
        let mut cols = line.split('\t');
        let feature_id = cols.next().unwrap_or_default().to_string();
        let gene_name = cols.next().unwrap_or_default().to_string();
        let feature_type = cols.next().unwrap_or_default().to_string();

        if gene_name.is_empty() {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "features.tsv row {} missing gene_name",
                lineno + 1
            )));
        }

        let feature_type = if feature_type.is_empty() {
            "Gene Expression".to_string()
        } else {
            feature_type
        };

        out.push(FeatureRowRaw {
            feature_id,
            gene_name,
            feature_type,
        });
    }

    Ok(out)
}

/// Builds canonical feature table sorted by `gene_name`.
pub(crate) fn build_feature_table(
    rows: Vec<FeatureRowRaw>,
) -> Result<FeatureBuildResult, SpatialIoError> {
    let mut indexed: Vec<(usize, FeatureRowRaw)> = rows.into_iter().enumerate().collect();
    indexed.sort_by(|a, b| a.1.gene_name.cmp(&b.1.gene_name));

    let mut out = Vec::with_capacity(indexed.len());
    let mut old_to_new = vec![0_u32; indexed.len()];
    let mut last_name: Option<String> = None;

    for (new_id, (old_idx, raw)) in indexed.into_iter().enumerate() {
        if last_name.as_deref().is_some_and(|n| n == raw.gene_name) {
            return Err(SpatialIoError::DuplicateGene(raw.gene_name));
        }
        last_name = Some(raw.gene_name.clone());

        old_to_new[old_idx] = new_id as u32;
        out.push(FeatureRow {
            gene_id: new_id as u32,
            feature_id: raw.feature_id,
            gene_name: raw.gene_name,
            feature_type: raw.feature_type,
        });
    }

    Ok(FeatureBuildResult {
        table: FeatureTable::new(out),
        old_to_new,
    })
}
