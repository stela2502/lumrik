use clap::Args;
use std::path::PathBuf;

use crate::{Chemistry, Grammar, PrimerDetector};

#[derive(Debug, Clone, Args)]
pub struct PrimerCli {
    /// Preset single-cell chemistry.
    ///
    /// Ignored if --primer-structure is supplied.
    #[arg(long, value_enum, num_args = 1.., default_values_t = [Chemistry::default()])]
    pub chemistry: Vec<Chemistry>,

    /// Custom primer/read structure grammar.
    ///
    /// Overrides --chemistry.
    #[arg(long)]
    pub primer_structure: Option<String>,

    /// Optional line-delimited cell-barcode whitelist.
    ///
    /// For TENX_CELL this replaces the built-in chemistry whitelist. For a
    /// custom grammar, exactly one CELL:N operation is required.
    #[arg(long, value_name = "FILE")]
    pub whitelist: Option<PathBuf>,

    /// Maximum barcode mismatches corrected by --whitelist / whitelist-backed chemistries.
    #[arg(long, default_value_t = 1)]
    pub whitelist_mismatches: u32,

    /// Also search the reverse-complement orientation.
    #[arg(long, default_value_t = true)]
    pub detect_reverse_complement: bool,
}

impl PrimerCli {
    pub fn detector(&self) -> Result<PrimerDetector, String> {
        let mut detector = if let Some(structure) = self.primer_structure.as_deref() {
            PrimerDetector::from_grammar(Grammar::parse("custom", structure)?)?
        } else {
            PrimerDetector::from_chemistries(self.chemistry.iter().copied())?
        };
        if let Some(path) = self.whitelist.as_deref() {
            detector = detector.with_whitelist_path(path, self.whitelist_mismatches)?;
        }

        detector = detector.with_reverse_complement_detection(self.detect_reverse_complement);

        Ok(detector)
    }
}
