# .kira-spatial.bin Format (v2)

This document defines the binary on-disk layout for `.kira-spatial.bin`.

## Endianness and Alignment

- Endianness: **little-endian only**. Building on big-endian targets is a compile error.
- Section payload alignment: every section payload starts at a 64-byte aligned file offset.
- Header and section table are at fixed offsets and are not required to be 64-byte aligned themselves.

## File Layout

1. Fixed header (64 bytes)
2. Section table (`section_count` entries, 24 bytes each)
3. Section payloads in canonical id order, each padded to next 64-byte boundary before payload start

## Header (64 bytes)

| Offset | Field           | Type      | Notes                                 |
|-------:|-----------------|-----------|---------------------------------------|
|   0..8 | `magic`         | `[u8; 8]` | Always `KIRASPAT`                     |
|   8..10 | `version`      | `u16`     | Currently `2`                         |
|  10..12 | `section_count`| `u16`     | `>= MIN_SECTION_COUNT (5)`, `<= 1024` |
|  12..16 | reserved       | 4 bytes   | Must be zero                          |
|  16..32 | `dataset_hash` | `[u8;16]` | BLAKE3 leading 16 bytes (see below)   |
|  32..64 | padding        | 32 bytes  | Must be zero                          |

## Section Table

Each entry is exactly 24 bytes:

| Offset | Field    | Type   | Notes                                                                |
|-------:|----------|--------|----------------------------------------------------------------------|
|   0..2 | `id`     | `u16`  | Section identifier                                                   |
|   2..4 | `flags`  | `u16`  | Bit 0 = section payload bytes are zstd-compressed                    |
|   4..12| `offset` | `u64`  | File offset of payload; must be 64-byte aligned (when length > 0)    |
|  12..20| `length` | `u64`  | On-disk byte length of payload (compressed length when bit 0 is set) |
|  20..24| `crc32`  | `u32`  | CRC-32/IEEE over the on-disk payload bytes                           |

Mandatory section ids: `1` SpatialDomain, `2` CSR, `3` FeatureTable, `4` MetadataCore, `5` MetadataJSON.

Readers **must** ignore unknown section ids (forward compatibility); they must still verify the per-section CRC and detect overlaps.

## Section 1: SpatialDomain Payload

- `n_bins: u32`
- `coord_system: u8` (`0=Grid, 1=Pixel, 2=Micron`)
- `bin_level: u8`
- `flags: u16` (bit 0 = `has_grid`)

Arrays (little-endian, tightly packed):

- `x: [f32; n_bins]`
- `y: [f32; n_bins]`
- if `has_grid`:
  - `grid_row: [u32; n_bins]`
  - `grid_col: [u32; n_bins]`
- `bin_id: [u32; n_bins]` (always `0..n_bins`)

Tissue mask (platform-independent):

- `n_bits: u64` — exactly `n_bins`
- `raw_bytes_len: u64` — must be a multiple of 8
- `raw_bytes: [u64; raw_bytes_len / 8]` — bit `i` is in word `i / 64` at bit position `i % 64` (LSB0)

## Section 2: CSR Payload

- `n_bins: u32`
- `n_genes: u32`
- `nnz: u64`
- `normalized: u8`
- `flags: u8` (bit 0 = `INDPTR_U32`)
- `reserved: [u8; 6]`

Arrays:

- if `INDPTR_U32`: `indptr: [u32; n_bins + 1]`, otherwise `indptr: [u64; n_bins + 1]`
- `indices: [u32; nnz]`
- `data: [f32; nnz]`

When `normalized` is `0`, `data` values must be finite and non-negative. When `1`, only finite is required.

## Section 3: FeatureTable Payload

- `n_genes: u32`

For each row in ascending `gene_id`:

- `gene_id: u32`
- `feature_id_len: u32`, `feature_id_bytes: [u8]` (UTF-8; typically Ensembl or vendor id; may be empty)
- `gene_name_len: u32`, `gene_name_bytes: [u8]` (UTF-8; canonical sort key)
- `feature_type_len: u32`, `feature_type_bytes: [u8]` (UTF-8)

## Section 4: MetadataCore Payload

- `dataset_name_len: u32`, `dataset_name_bytes`
- `source_format_len: u32`, `source_format_bytes`
- `bin_level: u8`
- `n_bins: u32`
- `n_genes: u32`
- `nnz: u64`
- `coord_system: u8`
- `normalized: u8`

The dataset hash lives only in the file header — keeping a single source of truth.

## Section 5: MetadataJSON Payload

- `json_len: u64`
- `json_bytes: [u8; json_len]` (canonical JSON, compact)

Numeric values may be integers or finite floats. Floats are serialized through `serde_json` (Ryū-style), which is deterministic for the same input value. NaN / Infinity are rejected.

## Dataset Hash

The dataset hash is the **leading 16 bytes** of a BLAKE3 digest computed over a fixed sequence of canonical, platform-independent byte streams, in this order:

1. `u16` SECTION_ID_SPATIAL_DOMAIN, then canonical SpatialDomain bytes (tissue mask streamed as `u64` LE words regardless of host pointer width)
2. `u16` SECTION_ID_CSR, then canonical CSR bytes (indptr always streamed as `u64` LE words, independent of the on-disk u32/u64 storage variant)
3. `u16` SECTION_ID_FEATURE_TABLE, then canonical FeatureTable bytes
4. `u16` SECTION_ID_META_CORE, then MetadataCore fields (no `dataset_hash` field exists)
5. `u16` SECTION_ID_META_JSON, then `u64` json length, then canonical JSON bytes

128-bit collision resistance is sufficient for dataset identity and deduplication. The hash is not a cryptographic MAC — use TLS / signatures for tamper detection in transit.

## Compression (optional, requires the `compression` feature)

Each section may be independently zstd-compressed. The writer's `CompressionPolicy::Zstd(level)` applies the same level to every section. Readers detect compression per-section via the `flags` bit and decompress before computing the CRC payload check (the CRC is computed over the *on-disk*, possibly-compressed bytes).

## Reader Validation Rules

- Header validation: magic, version (`2`), `section_count` within `[5, 1024]`, table bounds.
- Section table validation: mandatory sections (ids 1..5) must exist exactly once; section ranges must be in-bounds, 64-byte aligned, and non-overlapping. Unknown section ids are ignored.
- Per-section CRC32/IEEE must match (otherwise `SpatialIoError::CrcMismatch`).
- File size must not exceed `LoadConfig::memory_budget_mb` (early DoS guard).
- Declared sizes (`n_bins`, `nnz`, `json_len`) must fit the memory budget before any `Vec::with_capacity`.
- SpatialDomain validation: array lengths, canonical `bin_id` sequence, sorted coordinate order, tissue mask shape.
- CSR validation: strict `indptr` invariants, strictly increasing row indices, index bounds, finite (and non-negative unless `normalized`) values.
- FeatureTable validation: canonical `gene_id`, unique/sorted `gene_name`.
- Metadata validation: cross-section dimension consistency and canonical JSON bytes.
- Hash validation (when `LoadConfig::validate_hash` is true, the default): recomputed dataset hash must equal the header hash. `SpatialIoError::HashMismatch` otherwise.
