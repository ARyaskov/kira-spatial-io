use std::path::{Path, PathBuf};

use crate::config::LoadConfig;
use crate::error::SpatialIoError;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct HdBinnedLayout {
    pub root: PathBuf,
    pub bins: Vec<HdBinFolder>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct HdBinFolder {
    pub um: u32,
    pub level_code: u8,
    pub path: PathBuf,
}

pub fn discover_hd_binned(root: &Path) -> Result<Option<HdBinnedLayout>, SpatialIoError> {
    let binned_root = root.join("binned_outputs");
    if !binned_root.is_dir() {
        return Ok(None);
    }

    let mut bins = Vec::new();
    for entry in std::fs::read_dir(&binned_root)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }

        let name = entry.file_name();
        let Some(name) = name.to_str() else {
            continue;
        };

        if let Some(um) = parse_bin_folder_um(name)
            && let Some(level_code) = level_code_from_um(um)
        {
            bins.push(HdBinFolder {
                um,
                level_code,
                path: entry.path(),
            });
        }
    }

    bins.sort_unstable_by_key(|b| (b.um, b.path.to_string_lossy().to_string()));

    if bins.is_empty() {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "no supported bin_*um folders found in {}",
            binned_root.display()
        )));
    }

    Ok(Some(HdBinnedLayout {
        root: binned_root,
        bins,
    }))
}

pub fn select_bin_folder(
    layout: &HdBinnedLayout,
    cfg: &LoadConfig,
) -> Result<HdBinFolder, SpatialIoError> {
    if let Some(code) = cfg.bin_level {
        let Some(requested_um) = um_from_level_code(code) else {
            return Err(SpatialIoError::UnsupportedFormat(format!(
                "unsupported bin_level code: {}",
                code
            )));
        };

        return layout
            .bins
            .iter()
            .find(|b| b.um == requested_um)
            .cloned()
            .ok_or_else(|| {
                SpatialIoError::UnsupportedFormat(format!(
                    "requested bin level not found: {}um",
                    requested_um
                ))
            });
    }

    if let Some(bin8) = layout.bins.iter().find(|b| b.um == 8) {
        return Ok(bin8.clone());
    }

    layout.bins.first().cloned().ok_or_else(|| {
        SpatialIoError::UnsupportedFormat("requested bin level not found".to_string())
    })
}

#[cfg(feature = "parquet")]
pub fn discover_barcode_mapping_parquet(root: &Path) -> Option<PathBuf> {
    [
        root.join("barcode_mappings.parquet"),
        root.join("barcode_mapping.parquet"),
    ]
    .into_iter()
    .find(|p| p.is_file())
}

pub fn level_code_from_um(um: u32) -> Option<u8> {
    match um {
        2 => Some(1),
        8 => Some(2),
        16 => Some(3),
        32 => Some(4),
        64 => Some(5),
        _ => None,
    }
}

pub fn um_from_level_code(code: u8) -> Option<u32> {
    match code {
        1 => Some(2),
        2 => Some(8),
        3 => Some(16),
        4 => Some(32),
        5 => Some(64),
        _ => None,
    }
}

fn parse_bin_folder_um(name: &str) -> Option<u32> {
    let prefix = "bin_";
    let suffix = "um";
    if !name.starts_with(prefix) || !name.ends_with(suffix) {
        return None;
    }
    let middle = &name[prefix.len()..name.len() - suffix.len()];
    middle.parse::<u32>().ok()
}
