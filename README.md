# kira-spatial-io

`kira-spatial-io` is a deterministic Rust IO library for spatial transcriptomics datasets.

## Purpose

- Deterministic ingest/export for 10x spatial datasets.
- Canonical in-memory model: SoA spatial domain + CSR expression matrix.
- Stable `.kira-spatial.bin` format for fast validated reload.

## Supported Formats

- 10x MTX/TSV + `spatial/` via `Dataset::open_10x`.
- Visium HD binned outputs (`binned_outputs/bin_*um`) via `Dataset::open_10x`.
- 10x H5 (`feature_slice.h5` and compatible `/matrix/*`) via `Dataset::open_h5`.
- `.kira-spatial.bin` import/export via `Dataset::from_kira_bin` / `Dataset::export_kira_bin`.
- Optional HD barcode mapping parquet (`barcode_mappings.parquet`) with `--features parquet`.

## Determinism Guarantees

- Canonical bin sorting: `(grid_row, grid_col)` or `(y, x)` with total float ordering.
- Canonical feature sorting: `gene_name` ascending; stable remapped `gene_id`.
- Canonical CSR rows: indices strictly increasing; duplicate `(bin,gene)` merged by sum.
- Canonical JSON metadata keys sorted recursively.
- Canonical dataset hash (`BLAKE3-128`) verified on binary read.
- Cross-run binary byte stability for same input + config.

See [docs/DETERMINISM.md](docs/DETERMINISM.md).

## Memory Budget Behavior

`LoadConfig.memory_budget_mb` is enforced in staged builders:

- MTX CSR build checks final-size estimate before allocations.
- H5 CSR build checks final-size estimate and bounds chunk buffers by remaining budget.
- Parquet mapping loader enforces row-envelope estimate and errors with `MemoryLimitExceeded`.

## Feature Flags

| Feature | Default | Purpose |
|---|---:|---|
| `parallel` | no | Enables optional `rayon` integration paths/tests. |
| `parquet` | no | Enables Visium HD barcode mapping parquet ingestion and typed accessor. |

## Versioning Policy

- Crate version follows semver (`1.x.y`).
- `.kira-spatial.bin` layout is frozen for v1.x (`KIRA_SPATIAL_BIN_VERSION = 1`).
- Backward compatibility is guaranteed within major version for stable API and binary layout.

See [docs/API_CONTRACT.md](docs/API_CONTRACT.md) and [docs/FORMAT.md](docs/FORMAT.md).

## Minimal Usage

```rust
use kira_spatial_io::{Dataset, LoadConfig, SpatialIoError};

fn load_export_reload() -> Result<(), SpatialIoError> {
    let cfg = LoadConfig::default();

    let ds = Dataset::open_10x("/data/sample_10x", cfg)?;
    ds.export_kira_bin("/data/sample_10x/sample.kira-spatial.bin")?;

    let reloaded = Dataset::from_kira_bin("/data/sample_10x/sample.kira-spatial.bin")?;
    assert_eq!(ds.metadata_core().dataset_hash, reloaded.metadata_core().dataset_hash);

    Ok(())
}
```
