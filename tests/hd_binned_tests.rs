use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn open_10x_hd_selects_explicit_bin_level() {
    let root = prepare_hd_fixture("hd_explicit");

    let ds = Dataset::open_10x(
        &root,
        LoadConfig {
            bin_level: Some(3),
            ..LoadConfig::default()
        },
    )
    .expect("open_10x explicit level");

    assert_eq!(ds.metadata_core().bin_level, 3);
    assert_eq!(ds.metadata_json()["tenx"]["hd"], true);
    assert_eq!(ds.metadata_json()["tenx"]["bin_size_um"], 16);
    assert_eq!(ds.metadata_json()["tenx"]["bin_level_code"], 3);
    assert_eq!(ds.metadata_json()["source"]["layout"], "visium-hd-binned");

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn open_10x_hd_default_prefers_8um() {
    let root = prepare_hd_fixture("hd_default");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open_10x default");

    assert_eq!(ds.metadata_core().bin_level, 2);
    assert_eq!(ds.metadata_json()["tenx"]["bin_size_um"], 8);
    assert_eq!(ds.metadata_json()["tenx"]["bin_level_code"], 2);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn open_10x_hd_missing_requested_level_errors() {
    let root = temp_dir_path("hd_missing_level");
    write_bin_dataset(&root, 8).expect("write bin 8");

    let err = Dataset::open_10x(
        &root,
        LoadConfig {
            bin_level: Some(1),
            ..LoadConfig::default()
        },
    )
    .expect_err("should fail missing requested level");

    assert!(
        err.to_string().contains("requested bin level not found"),
        "actual error: {}",
        err
    );

    fs::remove_dir_all(&root).expect("cleanup");
}

fn prepare_hd_fixture(tag: &str) -> PathBuf {
    let root = temp_dir_path(tag);
    write_bin_dataset(&root, 2).expect("write bin 2");
    write_bin_dataset(&root, 8).expect("write bin 8");
    write_bin_dataset(&root, 16).expect("write bin 16");
    root
}

fn write_bin_dataset(root: &Path, um: u32) -> std::io::Result<()> {
    let bin_dir = root.join("binned_outputs").join(format!("bin_{}um", um));
    let matrix_dir = bin_dir.join("filtered_feature_bc_matrix");
    let spatial_dir = bin_dir.join("spatial");

    fs::create_dir_all(&matrix_dir)?;
    fs::create_dir_all(&spatial_dir)?;

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n")?;
    fs::write(
        matrix_dir.join("features.tsv"),
        "gene_b\tGeneB\tGene Expression\ngene_a\tGeneA\tGene Expression\n",
    )?;

    let mtx = format!(
        "%%MatrixMarket matrix coordinate integer general\n2 2 3\n2 1 {}\n1 2 {}\n2 2 1\n",
        um,
        um + 1
    );
    fs::write(matrix_dir.join("matrix.mtx"), mtx)?;

    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )?;

    Ok(())
}

fn temp_dir_path(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
