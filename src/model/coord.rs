/// Coordinate system used by [`SpatialDomain`](crate::SpatialDomain).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordSystem {
    /// Integer grid coordinates.
    Grid,
    /// Pixel-space coordinates.
    Pixel,
    /// Micron-space coordinates.
    Micron,
}
