use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use hdf5::types::VarLenUnicode;
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn open_h5_loads_tiny_10x_layout() {
    let root = temp_dir_path("h5_ok");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    let h5_path = root.join("feature_slice.h5");
    write_tiny_h5(&h5_path).expect("write h5");

    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )
    .expect("spatial");

    let ds = Dataset::open_h5(&root, LoadConfig::default()).expect("open_h5");
    assert_eq!(ds.metadata_core().source_format, "10x-h5");
    assert_eq!(ds.metadata_core().n_bins, 2);
    assert_eq!(ds.metadata_core().n_genes, 2);
    assert_eq!(ds.expression_csr().n_bins, 2);
    assert_eq!(ds.expression_csr().n_genes, 2);
    assert_eq!(ds.expression_csr().indptr, vec![0, 1, 3]);
    assert_eq!(ds.expression_csr().indices, vec![0, 0, 1]);
    assert_eq!(ds.expression_csr().data, vec![4.0, 1.0, 7.0]);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn open_h5_errors_on_missing_required_dataset() {
    let root = temp_dir_path("h5_missing");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    let h5_path = root.join("feature_slice.h5");
    write_missing_barcodes_h5(&h5_path).expect("write h5");

    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,1,1,10,10\n",
    )
    .expect("spatial");

    let err = Dataset::open_h5(&root, LoadConfig::default()).expect_err("must fail");
    assert!(err.to_string().contains("missing /matrix/barcodes dataset"));

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn open_h5_respects_low_memory_budget() {
    let root = temp_dir_path("h5_budget");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    let h5_path = root.join("feature_slice.h5");
    write_tiny_h5(&h5_path).expect("write h5");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )
    .expect("spatial");

    let err = Dataset::open_h5(
        &root,
        LoadConfig {
            memory_budget_mb: 0,
            bin_level: None,
            validate_strict: true,
        },
    )
    .expect_err("must fail with low budget");
    assert!(
        err.to_string().contains("memory limit exceeded"),
        "actual error: {}",
        err
    );

    fs::remove_dir_all(&root).expect("cleanup");
}

fn write_tiny_h5(path: &PathBuf) -> hdf5::Result<()> {
    let file = hdf5::File::create(path)?;
    let matrix = file.create_group("matrix")?;
    let features = matrix.create_group("features")?;

    let barcodes: Vec<VarLenUnicode> = ["BC1", "BC2"]
        .iter()
        .map(|s| s.parse().expect("unicode"))
        .collect();
    let names: Vec<VarLenUnicode> = ["GeneB", "GeneA"]
        .iter()
        .map(|s| s.parse().expect("unicode"))
        .collect();
    let feature_types: Vec<VarLenUnicode> = ["Gene Expression", "Gene Expression"]
        .iter()
        .map(|s| s.parse().expect("unicode"))
        .collect();

    matrix
        .new_dataset_builder()
        .with_data(&barcodes)
        .create("barcodes")?;
    features
        .new_dataset_builder()
        .with_data(&names)
        .create("name")?;
    features
        .new_dataset_builder()
        .with_data(&feature_types)
        .create("feature_type")?;

    let shape: [i64; 2] = [2, 2];
    matrix
        .new_dataset_builder()
        .with_data(&shape)
        .create("shape")?;

    let indptr: [i64; 3] = [0, 1, 3];
    let indices: [i32; 3] = [1, 0, 1];
    let data: [f32; 3] = [4.0, 7.0, 1.0];

    matrix
        .new_dataset_builder()
        .with_data(&indptr)
        .create("indptr")?;
    matrix
        .new_dataset_builder()
        .with_data(&indices)
        .create("indices")?;
    matrix
        .new_dataset_builder()
        .with_data(&data)
        .create("data")?;

    Ok(())
}

fn write_missing_barcodes_h5(path: &PathBuf) -> hdf5::Result<()> {
    let file = hdf5::File::create(path)?;
    let matrix = file.create_group("matrix")?;
    let features = matrix.create_group("features")?;

    let names: Vec<VarLenUnicode> = ["GeneA"]
        .iter()
        .map(|s| s.parse().expect("unicode"))
        .collect();
    let feature_types: Vec<VarLenUnicode> = ["Gene Expression"]
        .iter()
        .map(|s| s.parse().expect("unicode"))
        .collect();

    features
        .new_dataset_builder()
        .with_data(&names)
        .create("name")?;
    features
        .new_dataset_builder()
        .with_data(&feature_types)
        .create("feature_type")?;

    let shape: [i64; 2] = [1, 1];
    let indptr: [i64; 2] = [0, 0];
    let indices: [i32; 0] = [];
    let data: [f32; 0] = [];

    matrix
        .new_dataset_builder()
        .with_data(&shape)
        .create("shape")?;
    matrix
        .new_dataset_builder()
        .with_data(&indptr)
        .create("indptr")?;
    matrix
        .new_dataset_builder()
        .with_data(&indices)
        .create("indices")?;
    matrix
        .new_dataset_builder()
        .with_data(&data)
        .create("data")?;

    Ok(())
}

fn temp_dir_path(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
