use std::cmp::Ordering;

use crate::error::SpatialIoError;

/// Ensures a value is finite and non-negative. Use for raw count data.
pub fn ensure_f32_finite_nonneg(v: f32) -> Result<(), SpatialIoError> {
    if !v.is_finite() {
        return Err(SpatialIoError::InvalidFloat(
            "value is not finite".to_string(),
        ));
    }
    if v < 0.0 {
        return Err(SpatialIoError::InvalidFloat(
            "value is negative".to_string(),
        ));
    }
    Ok(())
}

/// Ensures a value is finite. Use for normalized data where negatives are legal.
pub fn ensure_f32_finite(v: f32) -> Result<(), SpatialIoError> {
    if !v.is_finite() {
        return Err(SpatialIoError::InvalidFloat(
            "value is not finite".to_string(),
        ));
    }
    Ok(())
}

/// Total-order float comparison wrapper.
pub fn total_cmp_f32(a: f32, b: f32) -> Ordering {
    a.total_cmp(&b)
}
