use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use hdf5::types::VarLenUnicode;
use kira_spatial_io::{Dataset, LoadConfig};

fn bench_mtx_ingest(c: &mut Criterion) {
    let root = temp_root("bench_mtx");
    let case = root.join("mtx_case");
    write_sparse_mtx_fixture(&case, 100_000, 5_000, 200_000).expect("fixture");

    let mut group = c.benchmark_group("ingest_mtx");
    group.measurement_time(Duration::from_secs(20));
    group.bench_function(BenchmarkId::new("open_10x", "100k_bins_5k_genes"), |b| {
        b.iter(|| {
            let ds = Dataset::open_10x(&case, LoadConfig::default()).expect("open_10x");
            assert_eq!(ds.metadata_core().source_format, "10x-mtx");
        });
    });
    group.finish();

    fs::remove_dir_all(&root).expect("cleanup");
}

fn bench_h5_ingest(c: &mut Criterion) {
    let root = temp_root("bench_h5");
    let case = root.join("h5_case");
    write_tiny_h5_fixture(&case).expect("fixture");

    let mut group = c.benchmark_group("ingest_h5");
    group.measurement_time(Duration::from_secs(10));
    group.bench_function(BenchmarkId::new("open_h5", "tiny"), |b| {
        b.iter(|| {
            let ds = Dataset::open_h5(&case, LoadConfig::default()).expect("open_h5");
            assert_eq!(ds.metadata_core().source_format, "10x-h5");
        });
    });
    group.finish();

    fs::remove_dir_all(&root).expect("cleanup");
}

fn bench_binary_reload(c: &mut Criterion) {
    let root = temp_root("bench_bin");
    let case = root.join("bin_case");
    write_sparse_mtx_fixture(&case, 100_000, 5_000, 200_000).expect("fixture");

    let ds = Dataset::open_10x(&case, LoadConfig::default()).expect("open_10x");
    let bin = case.join("dataset.kira-spatial.bin");
    ds.export_kira_bin(&bin).expect("export");

    let mut group = c.benchmark_group("reload_bin");
    group.measurement_time(Duration::from_secs(10));
    group.bench_function(
        BenchmarkId::new("from_kira_bin", "100k_bins_5k_genes"),
        |b| {
            b.iter(|| {
                let loaded = Dataset::from_kira_bin(&bin).expect("reload");
                assert_eq!(
                    loaded.metadata_core().source_format,
                    ds.metadata_core().source_format
                );
            });
        },
    );
    group.finish();

    fs::remove_dir_all(&root).expect("cleanup");
}

criterion_group!(
    benches,
    bench_mtx_ingest,
    bench_h5_ingest,
    bench_binary_reload
);
criterion_main!(benches);

fn write_sparse_mtx_fixture(
    root: &Path,
    n_bins: usize,
    n_genes: usize,
    nnz: usize,
) -> std::io::Result<()> {
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&matrix_dir)?;
    fs::create_dir_all(&spatial_dir)?;

    let mut barcodes = String::new();
    for i in 0..n_bins {
        barcodes.push_str(&format!("BC{idx:07}\n", idx = i + 1));
    }
    fs::write(matrix_dir.join("barcodes.tsv"), barcodes)?;

    let mut features = String::new();
    for i in 0..n_genes {
        features.push_str(&format!(
            "id{idx}\tGene{idx:05}\tGene Expression\n",
            idx = i + 1
        ));
    }
    fs::write(matrix_dir.join("features.tsv"), features)?;

    let mut mtx =
        format!("%%MatrixMarket matrix coordinate integer general\n{n_genes} {n_bins} {nnz}\n");
    let mut seed = 0xABCD_EF01_2345_6789_u64;
    for _ in 0..nnz {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let gene = (seed as usize % n_genes) + 1;
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let bin = (seed as usize % n_bins) + 1;
        let value = ((seed >> 32) % 10) + 1;
        mtx.push_str(&format!("{gene} {bin} {value}\n"));
    }
    fs::write(matrix_dir.join("matrix.mtx"), mtx)?;

    let mut spatial = String::from(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n",
    );
    for i in 0..n_bins {
        let row = i / 512;
        let col = i % 512;
        spatial.push_str(&format!(
            "BC{idx:07},1,{row},{col},{py},{px}\n",
            idx = i + 1,
            py = 1000 + row,
            px = 2000 + col
        ));
    }
    fs::write(spatial_dir.join("tissue_positions.csv"), spatial)?;

    Ok(())
}

fn write_tiny_h5_fixture(root: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&spatial_dir)?;

    let h5_path = root.join("feature_slice.h5");
    let file = hdf5::File::create(&h5_path)?;
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

    matrix
        .new_dataset_builder()
        .with_data(&[2_i64, 2_i64])
        .create("shape")?;
    matrix
        .new_dataset_builder()
        .with_data(&[0_i64, 1_i64, 3_i64])
        .create("indptr")?;
    matrix
        .new_dataset_builder()
        .with_data(&[1_i32, 0_i32, 1_i32])
        .create("indices")?;
    matrix
        .new_dataset_builder()
        .with_data(&[4.0_f32, 7.0_f32, 1.0_f32])
        .create("data")?;

    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )?;

    Ok(())
}

fn temp_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
