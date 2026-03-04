use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::binary::format::KIRA_SPATIAL_BIN_VERSION;
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn binary_header_contains_frozen_version_1() {
    let root = temp_root("version_header");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id1\tGene1\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,0,10,10\n",
    )
    .expect("spatial");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open");
    let out = root.join("version.kira-spatial.bin");
    ds.export_kira_bin(&out).expect("export");

    let raw = fs::read(&out).expect("read");
    let version = u16::from_le_bytes([raw[8], raw[9]]);
    assert_eq!(version, KIRA_SPATIAL_BIN_VERSION);

    fs::remove_dir_all(&root).expect("cleanup");
}

fn temp_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
