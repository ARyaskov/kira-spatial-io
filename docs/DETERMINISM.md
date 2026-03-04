# Determinism Contract

This document defines deterministic behavior requirements for `kira-spatial-io`.

## Ordering Rules

### Bins

- If grid coordinates exist, bins are ordered lexicographically by `(grid_row, grid_col)`.
- Otherwise bins are ordered by `(y, x)` using total float ordering (`f32::total_cmp`).
- The same permutation is applied to all SoA fields (`x`, `y`, `grid_*`, `bin_id`, `tissue_mask`).
- The sort step emits `old_to_new` mapping from barcode file order to canonical sorted bin row id.
- For Visium HD binned layouts, bin folder selection is deterministic: explicit `bin_level` code mapping is applied first; default prefers `8um`, otherwise smallest available `bin_*um`.

### Features / Genes

- Features are sorted by `gene_name` using UTF-8 byte ordering.
- Duplicate `gene_name` values are rejected.
- A deterministic remap `feature_old_to_new` converts MTX/H5 feature indices (file order) to canonical `gene_id`.

### CSR

- MatrixMarket entries `(feature, barcode, value)` are remapped to `(gene_id, sorted_bin_id, value)`.
- H5 sparse entries from `/matrix/{indices,indptr,data}` are remapped with the same `feature_old_to_new` and `barcode_old_to_new` mappings.
- CSR row order always matches sorted bin order.
- Within each row, `indices` are strictly increasing by `gene_id`.
- Duplicate `(row, gene_id)` entries are combined by deterministic summation.
- Invariants are validated:
  - `indptr.len() == n_bins + 1`
  - `indptr[0] == 0`
  - `indptr[last] == nnz`
  - `indptr` monotonic non-decreasing
  - `indices < n_genes`
  - `data` finite and non-negative

### Parquet Barcode Mapping (feature `parquet`)

- Mapping rows are sorted by `barcode` ascending (UTF-8 order), then by `(grid_row, grid_col, cell_id)` with `None` values ordered last.
- Duplicate barcodes are allowed and preserved in stable sorted order.
- Mapping data does not mutate `SpatialDomain` or `BinsCsr`; it enriches metadata and optional typed access only.

## Float Constraints

- No NaN or infinite values where numeric values are required.
- Non-negative constraints are checked at ingestion points where applicable.
- Total ordering comparisons use `f32::total_cmp` only.

## JSON Metadata Canonicalization

- Object keys are sorted lexicographically.
- Array element ordering is preserved unless a field-level contract specifies sorting.
- Numeric serialization strategy is deterministic and explicit (integers only in metadata JSON payloads).

## Hashing

- Dataset hash is derived from canonicalized, contract-defined payload components.
- Hash algorithm: BLAKE3 (with deterministic truncation/representation policy).
