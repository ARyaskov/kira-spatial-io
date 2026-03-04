# .kira-spatial.bin Format

**Binary layout frozen for v1.x. Backward compatible within major version.**

This document defines the binary on-disk layout for `.kira-spatial.bin`.

## Endianness and Alignment

- Endianness: little-endian only.
- Section payload alignment: every section payload starts at a 64-byte aligned file offset.
- Header and section table are not required to be 64-byte aligned.

## File Layout

1. Fixed header (64 bytes)
2. Section table (`section_count` entries)
3. Section payloads in fixed id order, each padded to next 64-byte boundary before payload start

## Header (64 bytes)

- `magic: [u8; 8]` = `KIRASPAT`
- `version: u16` = `1`
- `endian: u8` = `1` (little-endian)
- `section_count: u16` = `5`
- `dataset_hash: [u8; 16]` (BLAKE3-128)
- zero padding to 64 bytes

## Section Table

Each entry is exactly 18 bytes:

- `id: u16`
- `offset: u64`
- `length: u64`

Section ids:

- `1` SpatialDomain
- `2` CSR
- `3` FeatureTable
- `4` MetadataCore
- `5` MetadataJSON

## Section 1: SpatialDomain Payload

- `n_bins: u32`
- `coord_system: u8` (`0=Grid, 1=Pixel, 2=Micron`)
- `bin_level: u8`
- `flags: u16` (bit0 = has_grid)

Arrays:

- `x: [f32; n_bins]`
- `y: [f32; n_bins]`
- if `has_grid`:
  - `grid_row: [u32; n_bins]`
  - `grid_col: [u32; n_bins]`
- `bin_id: [u32; n_bins]`

Tissue mask:

- `n_bits: u64`
- `raw_bytes_len: u64`
- `raw_bytes: [u8; raw_bytes_len]` (`BitVec::as_raw_slice()` bytes in LE word order)

## Section 2: CSR Payload

- `n_bins: u32`
- `n_genes: u32`
- `nnz: u64`
- `normalized: u8`
- `reserved: [u8; 7]`

Arrays:

- `indptr: [u64; n_bins + 1]`
- `indices: [u32; nnz]`
- `data: [f32; nnz]`

## Section 3: FeatureTable Payload

- `n_genes: u32`

For each row in ascending `gene_id`:

- `gene_id: u32`
- `gene_name_len: u32`
- `gene_name_bytes: [u8; gene_name_len]`
- `feature_type_len: u32`
- `feature_type_bytes: [u8; feature_type_len]`

## Section 4: MetadataCore Payload

- `dataset_name_len: u32`, `dataset_name_bytes`
- `source_format_len: u32`, `source_format_bytes`
- `bin_level: u8`
- `n_bins: u32`
- `n_genes: u32`
- `nnz: u64`
- `coord_system: u8` (`0=Grid, 1=Pixel, 2=Micron`)
- `normalized: u8`
- `dataset_hash: [u8; 16]`

## Section 5: MetadataJSON Payload

- `json_len: u64`
- `json_bytes: [u8; json_len]` (canonical JSON, compact)

## Dataset Hash Canonicalization

The header and MetadataCore `dataset_hash` values are BLAKE3 truncated to 16 bytes over canonical payload components:

1. SpatialDomain canonical bytes
2. CSR canonical bytes
3. FeatureTable canonical bytes
4. MetadataCore fields excluding `dataset_hash`
5. Canonical JSON bytes

Numeric values are encoded LE, strings are `u32` length-prefixed UTF-8 bytes, and section identity tags are included in hash stream.

## Reader Validation Rules

- Header validation: magic, version (`1`), endian (`1`), and section table bounds.
- Section table validation: mandatory sections `1..5` must exist exactly once; section ranges must be in-bounds, 64-byte aligned, and non-overlapping.
- SpatialDomain validation: array lengths, canonical `bin_id` sequence, sorted coordinate order, and tissue mask shape.
- CSR validation: strict `indptr` invariants, strictly increasing row indices, index bounds, and finite non-negative values.
- FeatureTable validation: canonical `gene_id`, unique/sorted `gene_name`.
- Metadata validation: cross-section dimension consistency and canonical JSON bytes.
- Hash validation: recomputed dataset hash must match both header hash and metadata core hash.
