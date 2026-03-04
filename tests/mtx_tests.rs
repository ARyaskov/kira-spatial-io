use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn unsorted_triplets_become_sorted_by_gene_in_each_row() {
    let root = temp_root("mtx_unsorted");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id1\tA\tGene Expression\nid2\tB\tGene Expression\nid3\tC\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 4\n2 1 5\n1 1 7\n3 2 1\n2 2 3\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,0,10,10\nBC2,1,0,1,10,20\n",
    )
    .expect("spatial");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open");
    assert_eq!(ds.expression_csr().indptr, vec![0, 2, 4]);
    assert_eq!(ds.expression_csr().indices, vec![0, 1, 1, 2]);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn duplicate_gene_entries_are_summed_and_compacted() {
    let root = temp_root("mtx_dups");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id1\tA\tGene Expression\nid2\tB\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 4\n1 1 2\n1 1 3\n2 1 4\n2 2 1\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,0,10,10\nBC2,1,0,1,10,20\n",
    )
    .expect("spatial");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open");
    assert_eq!(ds.expression_csr().indptr, vec![0, 2, 3]);
    assert_eq!(ds.expression_csr().indices, vec![0, 1, 1]);
    assert_eq!(ds.expression_csr().data, vec![5.0, 4.0, 1.0]);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn bin_permutation_mapping_changes_row_assignment() {
    let root = temp_root("mtx_perm");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id1\tA\tGene Expression\nid2\tB\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 8\n2 2 9\n",
    )
    .expect("mtx");

    // Spatial order forces BC2 before BC1 by grid_col.
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,1,10,20\nBC2,1,0,0,10,10\n",
    )
    .expect("spatial");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open");
    assert_eq!(ds.expression_csr().indptr, vec![0, 1, 2]);
    assert_eq!(ds.expression_csr().indices, vec![1, 0]);
    assert_eq!(ds.expression_csr().data, vec![9.0, 8.0]);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn budget_enforcement_triggers_error() {
    let root = temp_root("mtx_budget");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\n").expect("barcodes");
    fs::write(matrix_dir.join("features.tsv"), "id1\tA\tGene Expression\n").expect("features");
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

    let err = Dataset::open_10x(
        &root,
        LoadConfig {
            memory_budget_mb: 0,
            bin_level: None,
            validate_strict: true,
        },
    )
    .expect_err("must fail");
    assert!(err.to_string().contains("memory limit exceeded"));

    fs::remove_dir_all(&root).expect("cleanup");
}

fn temp_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
