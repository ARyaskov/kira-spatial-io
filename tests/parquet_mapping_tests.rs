#![cfg(feature = "parquet")]

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use arrow2::array::{Array, Float32Array, UInt32Array, UInt64Array, Utf8Array};
use arrow2::chunk::Chunk;
use arrow2::datatypes::{DataType, Field, Schema};
use arrow2::io::parquet::write::{
    CompressionOptions, Encoding, FileWriter, RowGroupIterator, Version, WriteOptions, transverse,
};
use kira_spatial_io::{Dataset, LoadConfig};

#[test]
fn open_10x_hd_attaches_sorted_parquet_mapping_when_present() {
    let root = temp_dir_path("parquet_mapping_hd");
    write_bin_dataset(&root, 8).expect("write bin fixture");
    write_mapping_parquet(&root.join("barcode_mappings.parquet")).expect("write parquet");

    let ds = Dataset::open_10x(&root, LoadConfig::default()).expect("open_10x");
    let mapping = ds.barcode_mapping().expect("mapping present");

    assert_eq!(mapping.rows.len(), 3);
    assert_eq!(mapping.rows[0].barcode, "BC1");
    assert_eq!(mapping.rows[0].grid_row, Some(1));
    assert_eq!(mapping.rows[1].barcode, "BC1");
    assert_eq!(mapping.rows[1].grid_row, Some(2));
    assert_eq!(mapping.rows[2].barcode, "BC2");

    assert_eq!(ds.metadata_json()["tenx"]["barcode_mapping_present"], true);
    assert_eq!(ds.metadata_json()["tenx"]["barcode_mapping_rows"], 3);
    assert_eq!(ds.metadata_json()["tenx"]["mapping_has_duplicates"], true);

    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn parquet_mapping_respects_low_memory_budget() {
    let root = temp_dir_path("parquet_mapping_budget");
    write_bin_dataset(&root, 8).expect("write bin fixture");
    write_mapping_parquet(&root.join("barcode_mappings.parquet")).expect("write parquet");

    let err = Dataset::open_10x(
        &root,
        LoadConfig {
            memory_budget_mb: 0,
            ..LoadConfig::default()
        },
    )
    .expect_err("should fail with tiny budget");

    assert!(
        err.to_string().contains("memory limit exceeded"),
        "actual error: {err}"
    );

    fs::remove_dir_all(&root).expect("cleanup");
}

fn write_mapping_parquet(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let schema = Schema::from(vec![
        Field::new("barcode", DataType::Utf8, false),
        Field::new("cell_id", DataType::UInt64, true),
        Field::new("row", DataType::UInt32, true),
        Field::new("col", DataType::UInt32, true),
        Field::new("x", DataType::Float32, true),
        Field::new("y", DataType::Float32, true),
    ]);

    let chunk = Chunk::new(vec![
        Box::new(Utf8Array::<i32>::from_slice(["BC2", "BC1", "BC1"])) as Box<dyn Array>,
        Box::new(UInt64Array::from([Some(200), Some(101), Some(100)])) as Box<dyn Array>,
        Box::new(UInt32Array::from([Some(9), Some(2), Some(1)])) as Box<dyn Array>,
        Box::new(UInt32Array::from([Some(3), Some(2), Some(1)])) as Box<dyn Array>,
        Box::new(Float32Array::from([Some(50.0), Some(20.0), Some(10.0)])) as Box<dyn Array>,
        Box::new(Float32Array::from([Some(51.0), Some(21.0), Some(11.0)])) as Box<dyn Array>,
    ]);

    let options = WriteOptions {
        write_statistics: false,
        version: Version::V2,
        compression: CompressionOptions::Uncompressed,
        data_pagesize_limit: None,
    };

    let encodings = schema
        .fields
        .iter()
        .map(|f| transverse(&f.data_type, |_| Encoding::Plain))
        .collect::<Vec<_>>();

    let row_groups =
        RowGroupIterator::try_new(std::iter::once(Ok(chunk)), &schema, options, encodings)?;

    let file = std::fs::File::create(path)?;
    let mut writer = FileWriter::try_new(file, schema, options)?;
    for row_group in row_groups {
        writer.write(row_group?)?;
    }
    writer.end(None)?;
    Ok(())
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

    let mtx = "%%MatrixMarket matrix coordinate integer general\n2 2 3\n2 1 5\n1 2 6\n2 2 1\n";
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
