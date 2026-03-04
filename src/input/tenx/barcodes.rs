use std::collections::HashSet;
use std::io::BufRead;
use std::path::Path;

use crate::error::SpatialIoError;
use crate::input::util::open_text_maybe_gz;

pub fn load_barcodes(path: &Path) -> Result<Vec<String>, SpatialIoError> {
    let mut seen = HashSet::new();
    let mut barcodes = Vec::new();

    let reader = open_text_maybe_gz(path)?;
    for line in reader.lines() {
        let barcode = line?;
        if !seen.insert(barcode.clone()) {
            return Err(SpatialIoError::DuplicateBarcode(barcode));
        }
        barcodes.push(barcode);
    }

    Ok(barcodes)
}
