//! Binary format constants and section-table descriptors for `.kira-spatial.bin` v2.

/// File magic prefix.
pub const MAGIC: [u8; 8] = *b"KIRASPAT";

/// Current binary format version.
pub const KIRA_SPATIAL_BIN_VERSION: u16 = 2;
/// Convenience alias for the format version.
pub const VERSION: u16 = KIRA_SPATIAL_BIN_VERSION;

/// Fixed header byte size.
pub const HEADER_SIZE: u64 = 64;
/// Section table entry byte size (v2).
pub const SECTION_ENTRY_SIZE: u64 = 24;
/// Minimum mandatory section count.
pub const MIN_SECTION_COUNT: u16 = 5;
/// Upper bound on section count accepted by the reader (defensive cap).
pub const MAX_SECTION_COUNT: u16 = 1024;

/// Section identifier for the spatial domain payload.
pub const SECTION_ID_SPATIAL_DOMAIN: u16 = 1;
/// Section identifier for the expression CSR payload.
pub const SECTION_ID_CSR: u16 = 2;
/// Section identifier for the feature table payload.
pub const SECTION_ID_FEATURE_TABLE: u16 = 3;
/// Section identifier for the fixed metadata payload.
pub const SECTION_ID_META_CORE: u16 = 4;
/// Section identifier for the canonical JSON metadata payload.
pub const SECTION_ID_META_JSON: u16 = 5;

/// Mandatory section identifiers, in canonical write order.
pub const MANDATORY_SECTION_IDS: [u16; 5] = [
    SECTION_ID_SPATIAL_DOMAIN,
    SECTION_ID_CSR,
    SECTION_ID_FEATURE_TABLE,
    SECTION_ID_META_CORE,
    SECTION_ID_META_JSON,
];

/// Section flag: payload bytes on disk are zstd-compressed.
pub const SECTION_FLAG_ZSTD: u16 = 1 << 0;

/// CSR section flag: indptr is encoded as `u32` instead of `u64`.
pub const CSR_FLAG_INDPTR_U32: u8 = 1 << 0;

/// Parsed in-memory representation of the fixed file header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Header {
    /// Binary format version.
    pub version: u16,
    /// Number of entries in the section table.
    pub section_count: u16,
    /// Canonical dataset hash (BLAKE3 truncated to 16 leading bytes).
    pub dataset_hash: [u8; 16],
}

impl Header {
    /// Creates a current-version header for the supplied dataset hash.
    pub fn new(section_count: u16, dataset_hash: [u8; 16]) -> Self {
        Self {
            version: VERSION,
            section_count,
            dataset_hash,
        }
    }
}

/// A single section-table entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionEntry {
    /// Section identifier.
    pub id: u16,
    /// Section flags.
    pub flags: u16,
    /// Byte offset from file start.
    pub offset: u64,
    /// Section byte length on disk (possibly compressed).
    pub length: u64,
    /// CRC-32/IEEE of on-disk payload bytes.
    pub crc32c: u32,
}

impl SectionEntry {
    /// Builds a section-table entry.
    pub fn new(id: u16, flags: u16, offset: u64, length: u64, crc32c: u32) -> Self {
        Self {
            id,
            flags,
            offset,
            length,
            crc32c,
        }
    }
}
