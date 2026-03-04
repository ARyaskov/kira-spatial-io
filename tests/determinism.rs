use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::binary::hash::compute_dataset_hash;
use kira_spatial_io::determinism::json::write_canonical_json;
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn cross_run_binary_bytes_are_identical() {
    let root = fixture_root("determinism_bytes");
    write_tiny_mtx_fixture(&root).expect("fixture");

    let ds1 = Dataset::open_10x(&root, LoadConfig::default()).expect("open 1");
    let out1 = root.join("out_1.kira-spatial.bin");
    ds1.export_kira_bin(&out1).expect("export 1");

    let ds2 = Dataset::open_10x(&root, LoadConfig::default()).expect("open 2");
    let out2 = root.join("out_2.kira-spatial.bin");
    ds2.export_kira_bin(&out2).expect("export 2");

    let b1 = fs::read(&out1).expect("read out1");
    let b2 = fs::read(&out2).expect("read out2");
    assert_eq!(b1, b2);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn hash_is_stable_after_reload() {
    let root = fixture_root("determinism_hash");
    write_tiny_mtx_fixture(&root).expect("fixture");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");
    let out = root.join("stable.kira-spatial.bin");
    ds.export_kira_bin(&out).expect("export");

    let loaded = Dataset::from_kira_bin(&out).expect("reload");

    let mut json_bytes = Vec::new();
    write_canonical_json(&mut json_bytes, loaded.metadata_json()).expect("canonical json");

    let recomputed = compute_dataset_hash(
        loaded.spatial_domain(),
        loaded.expression_csr(),
        loaded.features(),
        loaded.metadata_core(),
        &json_bytes,
    )
    .expect("hash");

    assert_eq!(loaded.metadata_core().dataset_hash, recomputed);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[cfg(feature = "parallel")]
#[test]
fn parallel_thread_count_does_not_change_layout_or_hash() {
    use rayon::ThreadPoolBuilder;

    let root = fixture_root("determinism_parallel");
    write_tiny_mtx_fixture(&root).expect("fixture");

    let cfg = LoadConfig::default();

    let ds1 = ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("pool1")
        .install(|| Dataset::open_10x(&root, cfg.clone()))
        .expect("open one thread");

    let dsn = ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("poolN")
        .install(|| Dataset::open_10x(&root, cfg))
        .expect("open many threads");

    assert_eq!(ds1.expression_csr().indptr, dsn.expression_csr().indptr);
    assert_eq!(ds1.expression_csr().indices, dsn.expression_csr().indices);
    assert_eq!(ds1.expression_csr().data, dsn.expression_csr().data);
    assert_eq!(
        ds1.metadata_core().dataset_hash,
        dsn.metadata_core().dataset_hash
    );

    fs::remove_dir_all(&root).expect("cleanup");
}

fn write_tiny_mtx_fixture(root: &Path) -> std::io::Result<()> {
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");

    fs::create_dir_all(&matrix_dir)?;
    fs::create_dir_all(&spatial_dir)?;

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n")?;
    fs::write(
        matrix_dir.join("features.tsv"),
        "g2\tGeneB\tGene Expression\ng1\tGeneA\tGene Expression\n",
    )?;
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 3\n2 1 4\n1 2 5\n2 2 1\n",
    )?;
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )?;

    Ok(())
}

fn fixture_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
