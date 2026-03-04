# API Contract (v1.0)

This file defines the stable public API surface for `kira-spatial-io` v1.x.

## Stable Public Types

- `Dataset`
- `LoadConfig`
- `SpatialIoError`
- `CoordSystem`
- `SpatialDomain`
- `BinsCsr`
- `FeatureRow`
- `FeatureTable`
- `DatasetMetaCore`
- `BarcodeMappingRow` and `BarcodeMappingTable` (only with `feature = "parquet"`)

## Stable Public Methods

- `Dataset::open_10x`
- `Dataset::open_h5`
- `Dataset::from_kira_bin`
- `Dataset::export_kira_bin`
- `Dataset::spatial_domain`
- `Dataset::expression_csr`
- `Dataset::features`
- `Dataset::metadata_core`
- `Dataset::metadata_json`
- `Dataset::barcode_mapping` (only with `feature = "parquet"`)

## Stability Guarantees

- Semantic compatibility is guaranteed for all stable methods and stable type fields in v1.x.
- New optional fields may be added to `metadata_json` in backward-compatible form.
- Deterministic ordering, hashing, and serialization behavior is considered part of the contract.

## Binary Compatibility Rules

- `.kira-spatial.bin` is frozen at `KIRA_SPATIAL_BIN_VERSION = 1` for all v1.x releases.
- Readers must accept v1 files and reject unsupported versions explicitly.
- Mandatory sections (IDs 1..5) and payload ordering are stable in v1.x.
