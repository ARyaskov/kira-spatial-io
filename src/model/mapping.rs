/// Optional barcode mapping row used for HD metadata enrichment.
#[derive(Clone, Debug)]
pub struct BarcodeMappingRow {
    /// Barcode string.
    pub barcode: String,
    /// Optional cell identifier.
    pub cell_id: Option<u64>,
    /// Optional grid row.
    pub grid_row: Option<u32>,
    /// Optional grid column.
    pub grid_col: Option<u32>,
    /// Optional x coordinate.
    pub x: Option<f32>,
    /// Optional y coordinate.
    pub y: Option<f32>,
}

/// Deterministically sorted barcode mapping table.
#[derive(Clone, Debug)]
pub struct BarcodeMappingTable {
    /// Mapping rows.
    pub rows: Vec<BarcodeMappingRow>,
}

impl BarcodeMappingTable {
    /// Creates a mapping table.
    pub fn new(rows: Vec<BarcodeMappingRow>) -> Self {
        Self { rows }
    }
}
