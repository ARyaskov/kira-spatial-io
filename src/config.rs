/// Loading options controlling validation strictness and memory budget.
#[derive(Clone, Debug)]
pub struct LoadConfig {
    /// Maximum memory budget in megabytes for staged builders.
    pub memory_budget_mb: usize,
    /// Requested HD bin-level code, if applicable.
    pub bin_level: Option<u8>,
    /// Enables strict dimension/consistency checks.
    pub validate_strict: bool,
}

impl Default for LoadConfig {
    fn default() -> Self {
        Self {
            memory_budget_mb: 8192,
            bin_level: None,
            validate_strict: true,
        }
    }
}
