use bitvec::vec::BitVec;

use crate::error::SpatialIoError;
use crate::model::coord::CoordSystem;

/// Canonical spatial domain stored in SoA layout.
#[derive(Clone, Debug)]
pub struct SpatialDomain {
    /// X coordinates for each bin.
    pub x: Vec<f32>,
    /// Y coordinates for each bin.
    pub y: Vec<f32>,
    /// Optional grid row coordinates.
    pub grid_row: Option<Vec<u32>>,
    /// Optional grid column coordinates.
    pub grid_col: Option<Vec<u32>>,
    /// Canonical bin ids after sorting.
    pub bin_id: Vec<u32>,
    /// Tissue membership flags per bin.
    pub tissue_mask: BitVec,
    /// Coordinate system for spatial arrays.
    pub coord_system: CoordSystem,
    /// Bin level code.
    pub bin_level: u8,
}

impl SpatialDomain {
    #[allow(clippy::too_many_arguments)]
    /// Creates a spatial domain after validating consistent lengths.
    pub fn new(
        x: Vec<f32>,
        y: Vec<f32>,
        grid_row: Option<Vec<u32>>,
        grid_col: Option<Vec<u32>>,
        bin_id: Vec<u32>,
        tissue_mask: BitVec,
        coord_system: CoordSystem,
        bin_level: u8,
    ) -> Result<Self, SpatialIoError> {
        let n = x.len();
        if y.len() != n {
            return Err(SpatialIoError::DimensionMismatch(
                "x and y lengths differ".to_string(),
            ));
        }
        if bin_id.len() != n {
            return Err(SpatialIoError::DimensionMismatch(
                "bin_id length does not match coordinates".to_string(),
            ));
        }
        if tissue_mask.len() != n {
            return Err(SpatialIoError::DimensionMismatch(
                "tissue_mask length does not match coordinates".to_string(),
            ));
        }
        if let Some(ref row) = grid_row
            && row.len() != n
        {
            return Err(SpatialIoError::DimensionMismatch(
                "grid_row length does not match coordinates".to_string(),
            ));
        }
        if let Some(ref col) = grid_col
            && col.len() != n
        {
            return Err(SpatialIoError::DimensionMismatch(
                "grid_col length does not match coordinates".to_string(),
            ));
        }

        Ok(Self {
            x,
            y,
            grid_row,
            grid_col,
            bin_id,
            tissue_mask,
            coord_system,
            bin_level,
        })
    }

    /// Returns the number of bins.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Returns `true` when there are no bins.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}
