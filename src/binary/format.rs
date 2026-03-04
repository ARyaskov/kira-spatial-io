/// File magic prefix for `.kira-spatial.bin`.
pub const MAGIC: [u8; 8] = *b"KIRASPAT";
/// Frozen binary format version for all v1.x releases.
pub const KIRA_SPATIAL_BIN_VERSION: u16 = 1;
/// Legacy alias for format version.
pub const VERSION: u16 = KIRA_SPATIAL_BIN_VERSION;
/// Endianness discriminator for little-endian payloads.
pub const ENDIAN_LITTLE: u8 = 1;

/// Fixed header size in bytes.
pub const HEADER_SIZE: u64 = 64;
/// Single section table entry byte size.
pub const SECTION_ENTRY_SIZE: u64 = 18;
/// Mandatory section count for v1.
pub const SECTION_COUNT: u16 = 5;

/// Section identifier for the spatial domain payload.
pub const SECTION_ID_SPATIAL_DOMAIN: u16 = 1;
/// Section identifier for expression CSR payload.
pub const SECTION_ID_CSR: u16 = 2;
/// Section identifier for feature table payload.
pub const SECTION_ID_FEATURE_TABLE: u16 = 3;
/// Section identifier for fixed metadata payload.
pub const SECTION_ID_META_CORE: u16 = 4;
/// Section identifier for canonical JSON metadata payload.
pub const SECTION_ID_META_JSON: u16 = 5;

/// Parsed in-memory representation of the fixed file header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Header {
    /// Binary format version.
    pub version: u16,
    /// Endian marker.
    pub endian: u8,
    /// Number of entries in section table.
    pub section_count: u16,
    /// Canonical dataset hash (BLAKE3-128).
    pub dataset_hash: [u8; 16],
}

impl Header {
    /// Creates a v1 header for the provided dataset hash.
    pub fn new(dataset_hash: [u8; 16]) -> Self {
        Self {
            version: VERSION,
            endian: ENDIAN_LITTLE,
            section_count: SECTION_COUNT,
            dataset_hash,
        }
    }
}

/// A single section-table entry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionEntry {
    /// Section identifier.
    pub id: u16,
    /// Byte offset from file start.
    pub offset: u64,
    /// Section byte length.
    pub length: u64,
}

impl SectionEntry {
    /// Constructs a section-table entry.
    pub fn new(id: u16, offset: u64, length: u64) -> Self {
        Self { id, offset, length }
    }
}
