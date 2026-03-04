use bitvec::vec::BitVec;

use crate::error::SpatialIoError;
use crate::model::spatial_domain::SpatialDomain;

/// Result of canonical bin sorting.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SortBinsResult {
    /// Mapping from pre-sort index to post-sort index.
    pub old_to_new: Vec<u32>,
}

/// Canonically sorts bins and permutes all spatial SoA fields.
pub fn sort_bins(domain: &mut SpatialDomain) -> Result<SortBinsResult, SpatialIoError> {
    let n = domain.x.len();

    if domain.y.len() != n || domain.bin_id.len() != n || domain.tissue_mask.len() != n {
        return Err(SpatialIoError::DimensionMismatch(
            "spatial domain arrays have inconsistent lengths".to_string(),
        ));
    }

    match (&domain.grid_row, &domain.grid_col) {
        (Some(r), Some(c)) => {
            if r.len() != n || c.len() != n {
                return Err(SpatialIoError::DimensionMismatch(
                    "grid_row/grid_col lengths do not match coordinates".to_string(),
                ));
            }
        }
        (None, None) => {}
        _ => {
            return Err(SpatialIoError::DimensionMismatch(
                "grid_row and grid_col must be both present or both absent".to_string(),
            ));
        }
    }

    let mut perm: Vec<usize> = (0..n).collect();

    if let (Some(row), Some(col)) = (&domain.grid_row, &domain.grid_col) {
        perm.sort_unstable_by_key(|&i| (row[i], col[i]));
    } else {
        perm.sort_unstable_by(|&a, &b| {
            domain.y[a]
                .total_cmp(&domain.y[b])
                .then_with(|| domain.x[a].total_cmp(&domain.x[b]))
        });
    }

    let mut old_to_new = vec![0_u32; n];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        old_to_new[old_idx] = new_idx as u32;
    }

    domain.x = permute_vec(&domain.x, &perm);
    domain.y = permute_vec(&domain.y, &perm);
    domain.bin_id = (0..(n as u32)).collect();

    if let Some(ref mut grid_row) = domain.grid_row {
        *grid_row = permute_vec(grid_row, &perm);
    }
    if let Some(ref mut grid_col) = domain.grid_col {
        *grid_col = permute_vec(grid_col, &perm);
    }

    domain.tissue_mask = permute_bitvec(&domain.tissue_mask, &perm);

    Ok(SortBinsResult { old_to_new })
}

fn permute_vec<T: Copy>(data: &[T], perm: &[usize]) -> Vec<T> {
    perm.iter().map(|&i| data[i]).collect()
}

fn permute_bitvec(data: &BitVec, perm: &[usize]) -> BitVec {
    let mut out = BitVec::with_capacity(data.len());
    for &i in perm {
        out.push(data[i]);
    }
    out
}
