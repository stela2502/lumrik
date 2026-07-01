use clap::Args;

use crate::{Chemistry, Grammar, PrimerDetector};

#[derive(Debug, Clone, Args)]
pub struct PrimerCli {
    /// Preset single-cell chemistry.
    ///
    /// Ignored if --primer-structure is supplied.
    #[arg(long, value_enum, default_value_t = Chemistry::default())]
    pub chemistry: Chemistry,

    /// Custom primer/read structure grammar.
    ///
    /// Overrides --chemistry.
    #[arg(long)]
    pub primer_structure: Option<String>,

    /// Also search the reverse-complement orientation.
    #[arg(long, default_value_t = true)]
    pub detect_reverse_complement: bool,
}

impl PrimerCli {
    pub fn detector(&self) -> Result<PrimerDetector, String> {
        let grammar = if let Some(structure) = self.primer_structure.as_deref() {
            Grammar::parse("custom", structure)?
        } else {
            self.chemistry.grammar()?
        };

        let detector = PrimerDetector::from_grammar(grammar)?;

        Ok(detector)
    }
}
