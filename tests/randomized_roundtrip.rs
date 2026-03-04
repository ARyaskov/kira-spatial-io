use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn randomized_small_roundtrip_and_corruption_rejection() {
    let root = temp_root("randomized_roundtrip");
    fs::create_dir_all(&root).expect("root");

    let mut seed: u64 = 0x1234_5678_9abc_def0;

    for i in 0..8 {
        let ds_root = root.join(format!("case_{i}"));
        write_random_fixture(&ds_root, &mut seed).expect("fixture");

        let ds = Dataset::open_10x(&ds_root, LoadConfig::default()).expect("open");
        let bin = ds_root.join("case.kira-spatial.bin");
        ds.export_kira_bin(&bin).expect("export");

        let loaded = Dataset::from_kira_bin(&bin).expect("read");
        assert_eq!(ds.expression_csr().indptr, loaded.expression_csr().indptr);
        assert_eq!(ds.expression_csr().indices, loaded.expression_csr().indices);
        assert_eq!(ds.expression_csr().data, loaded.expression_csr().data);
        assert_ne!(loaded.metadata_core().dataset_hash, [0_u8; 16]);

        let mut raw = fs::read(&bin).expect("read raw");
        raw[0] ^= 0xFF;
        let corrupted = ds_root.join("corrupted.kira-spatial.bin");
        fs::write(&corrupted, raw).expect("write corrupted");
        assert!(Dataset::from_kira_bin(&corrupted).is_err());
    }

    fs::remove_dir_all(&root).expect("cleanup");
}

fn write_random_fixture(root: &Path, seed: &mut u64) -> std::io::Result<()> {
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&matrix_dir)?;
    fs::create_dir_all(&spatial_dir)?;

    let n_bins = 6usize;
    let n_genes = 5usize;

    let barcodes = (0..n_bins)
        .map(|i| format!("BC{idx:03}", idx = i + 1))
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(matrix_dir.join("barcodes.tsv"), format!("{barcodes}\n"))?;

    let mut features_lines = String::new();
    for i in 0..n_genes {
        features_lines.push_str(&format!(
            "id{idx}\tGene{idx}\tGene Expression\n",
            idx = i + 1
        ));
    }
    fs::write(matrix_dir.join("features.tsv"), features_lines)?;

    let mut triplets = Vec::new();
    for row in 1..=n_genes {
        for col in 1..=n_bins {
            let r = next_rand(seed) % 7;
            if r < 2 {
                let value = ((next_rand(seed) % 5) + 1) as u32;
                triplets.push((row, col, value));
            }
        }
    }
    if triplets.is_empty() {
        triplets.push((1, 1, 1));
    }

    let mut mtx = format!(
        "%%MatrixMarket matrix coordinate integer general\n{n_genes} {n_bins} {}\n",
        triplets.len()
    );
    for (r, c, v) in &triplets {
        mtx.push_str(&format!("{r} {c} {v}\n"));
    }
    fs::write(matrix_dir.join("matrix.mtx"), mtx)?;

    let mut spatial_csv = String::from(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n",
    );
    for i in (0..n_bins).rev() {
        let row = i as u32 / 3;
        let col = i as u32 % 3;
        let tissue = if i % 2 == 0 { 1 } else { 0 };
        spatial_csv.push_str(&format!(
            "BC{idx:03},{tissue},{row},{col},{py},{px}\n",
            idx = i + 1,
            py = 100 + i,
            px = 200 + i
        ));
    }
    fs::write(spatial_dir.join("tissue_positions.csv"), spatial_csv)?;

    Ok(())
}

fn next_rand(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    *seed
}

fn temp_root(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
