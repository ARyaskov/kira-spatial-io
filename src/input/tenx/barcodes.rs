use std::collections::HashSet;
use std::io::BufRead;
use std::path::Path;

use crate::error::{IoPathExt, SpatialIoError};
use crate::input::util::open_text_maybe_gz;

/// Loads a `barcodes.tsv(.gz)` file, rejecting duplicates.
pub fn load_barcodes(path: &Path) -> Result<Vec<String>, SpatialIoError> {
    let reader = open_text_maybe_gz(path)?;
    let mut barcodes: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::with_capacity(1024);

    for line in reader.lines() {
        let barcode = line.io_path(path)?;
        if !seen.insert(barcode.clone()) {
            return Err(SpatialIoError::DuplicateBarcode(barcode));
        }
        barcodes.push(barcode);
    }

    Ok(barcodes)
}
