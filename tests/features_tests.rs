use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn feature_table_sorted_and_reindexed() {
    let root = temp_root("features_sorted");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id_b\tB\tGene Expression\nid_a\tA\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 1\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,0,10,10\n",
    )
    .expect("spatial");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open");
    assert_eq!(ds.features().rows.len(), 2);
    assert_eq!(ds.features().rows[0].gene_name, "A");
    assert_eq!(ds.features().rows[1].gene_name, "B");
    assert_eq!(ds.features().rows[0].gene_id, 0);
    assert_eq!(ds.features().rows[1].gene_id, 1);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn feature_table_rejects_duplicate_gene_names() {
    let root = temp_root("features_duplicate");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "id1\tG\tGene Expression\nid2\tG\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 1 1\n1 1 1\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC1,1,0,0,10,10\n",
    )
    .expect("spatial");

    let err = Dataset::open_10x(&root, LoadConfig::default()).expect_err("duplicate should fail");
    assert!(err.to_string().contains("duplicate gene: G"));

    fs::remove_dir_all(&root).expect("cleanup");
}

fn temp_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
