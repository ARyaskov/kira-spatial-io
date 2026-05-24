use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::binary::format::{
    HEADER_SIZE, KIRA_SPATIAL_BIN_VERSION, MAGIC, MIN_SECTION_COUNT, SECTION_ENTRY_SIZE,
    SECTION_ID_META_JSON,
};
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn writer_emits_valid_header_and_aligned_sections() {
    let root = temp_dir_path("writer_v2_header");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&matrix_dir).expect("create matrix dir");
    fs::create_dir_all(&spatial_dir).expect("create spatial dir");

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

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open_10x");
    let out = root.join("out.kira-spatial.bin");
    ds.export_kira_bin(&out).expect("export");

    let bytes = fs::read(&out).expect("read output");
    assert!(bytes.len() > HEADER_SIZE as usize);

    // v2 header: magic[0..8] / version[8..10] / section_count[10..12] / reserved[12..16]
    // / hash[16..32] / pad[32..64].
    assert_eq!(&bytes[0..8], &MAGIC);
    assert_eq!(read_u16_le(&bytes, 8), KIRA_SPATIAL_BIN_VERSION);
    let section_count = read_u16_le(&bytes, 10);
    assert!(section_count >= MIN_SECTION_COUNT);

    let mut meta_json_entry_offset = 0_u64;
    let mut meta_json_entry_len = 0_u64;

    // v2 entry: id[0..2] / flags[2..4] / offset[4..12] / length[12..20] / crc[20..24].
    let table_start = HEADER_SIZE as usize;
    for i in 0..section_count as usize {
        let base = table_start + i * SECTION_ENTRY_SIZE as usize;
        let id = read_u16_le(&bytes, base);
        let flags = read_u16_le(&bytes, base + 2);
        let offset = read_u64_le(&bytes, base + 4);
        let length = read_u64_le(&bytes, base + 12);
        let crc = read_u32_le(&bytes, base + 20);

        assert_eq!(offset % 64, 0, "section {id} is not 64-byte aligned");
        assert!(length > 0, "section {id} length must be > 0");
        assert_eq!(flags, 0, "uncompressed write should have zero flags");
        assert!(crc != 0 || length == 0, "non-empty section needs CRC");

        if id == SECTION_ID_META_JSON {
            meta_json_entry_offset = offset;
            meta_json_entry_len = length;
        }
    }

    assert!(meta_json_entry_offset > 0);
    assert!(meta_json_entry_len >= 8);

    let json_section_start = meta_json_entry_offset as usize;
    let json_len = read_u64_le(&bytes, json_section_start) as usize;
    let json_start = json_section_start + 8;
    let json_end = json_start + json_len;
    assert!(json_end <= bytes.len());

    let json_text = std::str::from_utf8(&bytes[json_start..json_end]).expect("utf8 json");
    assert!(json_text.starts_with("{\"source\":{\"dataset_root\":"));
    assert!(json_text.contains("\"layout\":\"tenx-mtx\""));
    assert!(json_text.contains(
        "\"tenx\":{\"bin_level_code\":0,\"bin_size_um\":0,\"duplicate_policy\":\"sum-per-bin-gene\",\"has_spatial\":true,\"hd\":false}"
    ));

    fs::remove_dir_all(&root).expect("cleanup");
}

#[cfg(feature = "compression")]
#[test]
fn writer_zstd_compression_round_trips() {
    use kira_spatial_io::CompressionPolicy;

    let root = temp_dir_path("writer_zstd");
    let matrix_dir = root.join("filtered_feature_bc_matrix");
    let spatial_dir = root.join("spatial");
    fs::create_dir_all(&matrix_dir).expect("create matrix dir");
    fs::create_dir_all(&spatial_dir).expect("create spatial dir");

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

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open_10x");
    let out = root.join("out.kira-spatial.bin");
    ds.export_kira_bin_with(&out, CompressionPolicy::Zstd(3))
        .expect("export zstd");

    let loaded = Dataset::from_kira_bin(&out).expect("read zstd");
    assert_eq!(loaded.metadata_core().n_bins, ds.metadata_core().n_bins);
    assert_eq!(loaded.metadata_core().nnz, ds.metadata_core().nnz);
    assert_eq!(loaded.metadata_core().dataset_hash, ds.metadata_core().dataset_hash);

    fs::remove_dir_all(&root).expect("cleanup");
}

fn temp_dir_path(tag: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time")
        .as_nanos();
    std::env::temp_dir().join(format!("kira_spatial_io_{tag}_{ts}"))
}

fn read_u16_le(bytes: &[u8], offset: usize) -> u16 {
    let mut arr = [0_u8; 2];
    arr.copy_from_slice(&bytes[offset..offset + 2]);
    u16::from_le_bytes(arr)
}

fn read_u32_le(bytes: &[u8], offset: usize) -> u32 {
    let mut arr = [0_u8; 4];
    arr.copy_from_slice(&bytes[offset..offset + 4]);
    u32::from_le_bytes(arr)
}

fn read_u64_le(bytes: &[u8], offset: usize) -> u64 {
    let mut arr = [0_u8; 8];
    arr.copy_from_slice(&bytes[offset..offset + 8]);
    u64::from_le_bytes(arr)
}
