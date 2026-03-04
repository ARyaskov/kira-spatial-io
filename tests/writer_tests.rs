use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kira_spatial_io::binary::format::{
    ENDIAN_LITTLE, HEADER_SIZE, MAGIC, SECTION_COUNT, SECTION_ENTRY_SIZE, SECTION_ID_META_JSON,
};
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn writer_emits_valid_header_and_aligned_sections() {
    let root = temp_dir_path("writer_stage3");
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

    assert_eq!(&bytes[0..8], &MAGIC);
    assert_eq!(read_u16_le(&bytes, 8), 1);
    assert_eq!(bytes[10], ENDIAN_LITTLE);
    assert_eq!(read_u16_le(&bytes, 11), SECTION_COUNT);

    let mut meta_json_entry_offset = 0_u64;
    let mut meta_json_entry_len = 0_u64;

    let table_start = HEADER_SIZE as usize;
    for i in 0..SECTION_COUNT as usize {
        let base = table_start + i * SECTION_ENTRY_SIZE as usize;
        let id = read_u16_le(&bytes, base);
        let offset = read_u64_le(&bytes, base + 2);
        let length = read_u64_le(&bytes, base + 10);

        assert_eq!(offset % 64, 0, "section {id} is not 64-byte aligned");
        assert!(length > 0, "section {id} length must be > 0");

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

fn read_u64_le(bytes: &[u8], offset: usize) -> u64 {
    let mut arr = [0_u8; 8];
    arr.copy_from_slice(&bytes[offset..offset + 8]);
    u64::from_le_bytes(arr)
}
