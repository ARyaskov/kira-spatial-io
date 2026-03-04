use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::SpatialIoError;

pub(crate) fn open_text_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, SpatialIoError> {
    let file = File::open(path)?;
    let is_gz = path
        .extension()
        .and_then(OsStr::to_str)
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    if is_gz {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}
