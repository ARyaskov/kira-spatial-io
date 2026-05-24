use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::binary::format::HEADER_SIZE;
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn writer_reader_roundtrip_preserves_dataset_content() {
    let root = prepare_fixture("reader_roundtrip");
    let source = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");

    let out = root.join("sample.kira-spatial.bin");
    source.export_kira_bin(&out).expect("export");

    let loaded = Dataset::from_kira_bin(&out).expect("read bin");

    assert_eq!(loaded.spatial_domain().x, source.spatial_domain().x);
    assert_eq!(loaded.spatial_domain().y, source.spatial_domain().y);
    assert_eq!(
        loaded.spatial_domain().grid_row,
        source.spatial_domain().grid_row
    );
    assert_eq!(
        loaded.spatial_domain().grid_col,
        source.spatial_domain().grid_col
    );
    assert_eq!(
        loaded.spatial_domain().bin_id,
        source.spatial_domain().bin_id
    );
    assert_eq!(
        loaded.spatial_domain().tissue_mask,
        source.spatial_domain().tissue_mask
    );
    assert_eq!(
        loaded.spatial_domain().coord_system,
        source.spatial_domain().coord_system
    );
    assert_eq!(
        loaded.spatial_domain().bin_level,
        source.spatial_domain().bin_level
    );

    assert_eq!(
        loaded.expression_csr().indptr.to_u64_vec(),
        source.expression_csr().indptr.to_u64_vec()
    );
    assert_eq!(
        loaded.expression_csr().indices,
        source.expression_csr().indices
    );
    assert_eq!(loaded.expression_csr().data, source.expression_csr().data);
    assert_eq!(
        loaded.expression_csr().n_bins,
        source.expression_csr().n_bins
    );
    assert_eq!(
        loaded.expression_csr().n_genes,
        source.expression_csr().n_genes
    );
    assert_eq!(loaded.expression_csr().nnz, source.expression_csr().nnz);

    assert_eq!(loaded.features().rows.len(), source.features().rows.len());
    for (a, b) in loaded
        .features()
        .rows
        .iter()
        .zip(source.features().rows.iter())
    {
        assert_eq!(a.gene_id, b.gene_id);
        assert_eq!(a.feature_id, b.feature_id, "feature_id (Ensembl) preserved");
        assert_eq!(a.gene_name, b.gene_name);
        assert_eq!(a.feature_type, b.feature_type);
    }

    // The feature_id column (column 0 of features.tsv) must round-trip.
    // Fixture features in canonical (sorted-by-gene_name) order are gene_a -> GeneA, gene_b -> GeneB.
    assert_eq!(loaded.features().rows[0].feature_id, "gene_a");
    assert_eq!(loaded.features().rows[1].feature_id, "gene_b");

    assert_eq!(loaded.metadata_json(), source.metadata_json());
    assert_eq!(
        loaded.metadata_core().dataset_name,
        source.metadata_core().dataset_name
    );
    assert_eq!(
        loaded.metadata_core().source_format,
        source.metadata_core().source_format
    );
    assert_eq!(
        loaded.metadata_core().bin_level,
        source.metadata_core().bin_level
    );
    assert_eq!(loaded.metadata_core().n_bins, source.metadata_core().n_bins);
    assert_eq!(
        loaded.metadata_core().n_genes,
        source.metadata_core().n_genes
    );
    assert_eq!(loaded.metadata_core().nnz, source.metadata_core().nnz);
    assert_eq!(
        loaded.metadata_core().coord_system,
        source.metadata_core().coord_system
    );
    assert_eq!(
        loaded.metadata_core().normalized,
        source.metadata_core().normalized
    );
    assert_ne!(loaded.metadata_core().dataset_hash, [0_u8; 16]);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn reader_rejects_corrupted_magic_and_version() {
    let root = prepare_fixture("reader_corrupt_magic_version");
    let source = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");
    let out = root.join("sample.kira-spatial.bin");
    source.export_kira_bin(&out).expect("export");

    let mut bytes = fs::read(&out).expect("read bytes");
    bytes[0] = b'X';
    let bad_magic = root.join("bad_magic.kira-spatial.bin");
    fs::write(&bad_magic, &bytes).expect("write bad magic");
    let err = Dataset::from_kira_bin(&bad_magic).expect_err("must fail magic");
    assert!(err.to_string().contains("invalid magic"));

    let mut bytes2 = fs::read(&out).expect("read bytes2");
    // Set a definitely-unknown version.
    bytes2[8..10].copy_from_slice(&99_u16.to_le_bytes());
    let bad_version = root.join("bad_version.kira-spatial.bin");
    fs::write(&bad_version, &bytes2).expect("write bad version");
    let err = Dataset::from_kira_bin(&bad_version).expect_err("must fail version");
    assert!(err.to_string().contains("unsupported version"));

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn reader_rejects_section_alignment_and_hash_mismatch() {
    let root = prepare_fixture("reader_corrupt_align_hash");
    let source = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");
    let out = root.join("sample.kira-spatial.bin");
    source.export_kira_bin(&out).expect("export");

    // v2 section entry layout:
    //   0..2   id, 2..4   flags, 4..12  offset, 12..20 length, 20..24 crc32
    let mut bytes = fs::read(&out).expect("read bytes");
    let entry0_offset_field = HEADER_SIZE as usize + 4;
    bytes[entry0_offset_field..entry0_offset_field + 8]
        .copy_from_slice(&65_u64.to_le_bytes());
    let bad_align = root.join("bad_align.kira-spatial.bin");
    fs::write(&bad_align, &bytes).expect("write bad align");
    let err = Dataset::from_kira_bin(&bad_align).expect_err("must fail alignment");
    assert!(err.to_string().contains("not 64-byte aligned"));

    // Hash field in v2 header is at bytes 16..32.
    let mut bytes2 = fs::read(&out).expect("read bytes2");
    bytes2[16] ^= 0xFF;
    let bad_hash = root.join("bad_hash.kira-spatial.bin");
    fs::write(&bad_hash, &bytes2).expect("write bad hash");
    let err = Dataset::from_kira_bin(&bad_hash).expect_err("must fail hash");
    assert!(err.to_string().contains("dataset hash mismatch"));

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn reader_rejects_section_crc_mismatch() {
    let root = prepare_fixture("reader_corrupt_section_crc");
    let source = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");
    let out = root.join("sample.kira-spatial.bin");
    source.export_kira_bin(&out).expect("export");

    let mut bytes = fs::read(&out).expect("read bytes");
    // Section table starts at HEADER_SIZE; first entry id/flags/offset/length/crc occupy
    // bytes 64..88. Read the first section's offset and length, then flip a byte inside
    // its payload — that should trigger a per-section CRC mismatch.
    let offset = {
        let mut arr = [0_u8; 8];
        arr.copy_from_slice(&bytes[HEADER_SIZE as usize + 4..HEADER_SIZE as usize + 12]);
        u64::from_le_bytes(arr) as usize
    };
    let length = {
        let mut arr = [0_u8; 8];
        arr.copy_from_slice(&bytes[HEADER_SIZE as usize + 12..HEADER_SIZE as usize + 20]);
        u64::from_le_bytes(arr) as usize
    };
    assert!(length > 0);
    let target = offset + length / 2;
    bytes[target] ^= 0xAA;

    let bad_crc = root.join("bad_crc.kira-spatial.bin");
    fs::write(&bad_crc, &bytes).expect("write bad crc");
    let err = Dataset::from_kira_bin(&bad_crc).expect_err("must fail crc");
    assert!(err.to_string().contains("CRC mismatch"));

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn reader_rejects_oversized_file_against_budget() {
    let root = prepare_fixture("reader_budget");
    let source = Dataset::open_10x(&root, LoadConfig::default()).expect("open source");
    let out = root.join("sample.kira-spatial.bin");
    source.export_kira_bin(&out).expect("export");

    // Set the budget below the file size to verify the early file-size guard.
    let tiny_cfg = LoadConfig {
        memory_budget_mb: 0,
        ..LoadConfig::default()
    };
    let err = Dataset::from_kira_bin_with(&out, &tiny_cfg).expect_err("must fail budget");
    assert!(err.to_string().contains("memory budget") || err.to_string().contains("memory limit"));

    fs::remove_dir_all(&root).expect("cleanup");
}

fn prepare_fixture(tag: &str) -> PathBuf {
    let root = temp_dir_path(tag);
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&matrix_dir).expect("matrix dir");
    fs::create_dir_all(&spatial_dir).expect("spatial dir");

    fs::write(matrix_dir.join("barcodes.tsv"), "BC1\nBC2\n").expect("barcodes");
    fs::write(
        matrix_dir.join("features.tsv"),
        "gene_b\tGeneB\tGene Expression\ngene_a\tGeneA\tGene Expression\n",
    )
    .expect("features");
    fs::write(
        matrix_dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 3\n2 1 4\n1 2 7\n2 2 1\n",
    )
    .expect("mtx");
    fs::write(
        spatial_dir.join("tissue_positions.csv"),
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\nBC2,1,1,2,10,20\nBC1,0,1,1,11,19\n",
    )
    .expect("spatial");

    root
}

fn temp_dir_path(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}
