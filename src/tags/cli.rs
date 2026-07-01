use clap::{Args, ValueEnum};
use std::path::PathBuf;
use crate::tags::FastTagMapper;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum TagRead {
    #[value(name = "r1")]
    R1,

    #[value(name = "r2")]
    R2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BuiltinTagSet {
    #[value(name = "none")]
    None,

    #[value(name = "bd-human")]
    BdHuman,

    #[value(name = "bd-mouse")]
    BdMouse,
}

#[derive(Debug, Clone, Args)]
pub struct TagCli {
    #[arg(
        long,
        value_enum,
        default_value_t = BuiltinTagSet::None,
        value_name = "TAGSET",
        help = "Optional built-in DNA tag set to detect. Supported: none, bd-human, bd-mouse."
    )]
    pub tag_set: BuiltinTagSet,

    #[arg(
        long,
        value_name = "TSV/FASTA",
        help = "Optional custom DNA tag file. TSV format: name<TAB>sequence. FASTA is also accepted."
    )]
    pub custom_tag_file: Option<PathBuf>,

    #[arg(
        long,
        value_enum,
        default_value_t = TagRead::R1,
        value_name = "r1|r2",
        help = "Read used for side-channel tag detection."
    )]
    pub search_read: TagRead,

    #[arg(
        long,
        default_value = "feature_tags",
        value_name = "DIR",
        help = "Output directory for feature-tag quantification results."
    )]
    pub feature_tags: PathBuf,

    #[arg(
        long,
        default_value_t = 9,
        value_name = "BP",
        help = "K-mer size used for DNA tag matching."
    )]
    pub tag_kmer_size: usize,

    #[arg(
        long,
        default_value_t = 1,
        value_name = "N",
        help = "Step size while scanning read sequence for tag-matching kmers."
    )]
    pub tag_jump: usize,

    #[arg(
        long,
        default_value_t = 0,
        value_name = "OFFSET",
        help = "Start offset in the selected read for tag-matching kmer scans."
    )]
    pub tag_start: usize,

    #[arg(
        long,
        default_value_t = 3,
        value_name = "N",
        help = "Minimum winning kmer matches required to report a tag."
    )]
    pub tag_min_matches: usize,
}

impl TagCli {
    pub fn enabled(&self) -> bool {
        self.tag_set != BuiltinTagSet::None || self.custom_tag_file.is_some()
    }

    pub fn builtin_name(&self) -> Option<&'static str> {
        match self.tag_set {
            BuiltinTagSet::None => None,
            BuiltinTagSet::BdHuman => Some("bd-human"),
            BuiltinTagSet::BdMouse => Some("bd-mouse"),
        }
    }

    pub fn mapper(&self) -> anyhow::Result<Option<FastTagMapper>> {
        use anyhow::bail;

        let mapper = match (&self.custom_tag_file, self.tag_set) {
            (None, BuiltinTagSet::None) => {
                return Ok(None);
            }

            (Some(path), _) => {
                let ext = path
                    .extension()
                    .and_then(|x| x.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();

                if matches!(ext.as_str(), "fa" | "fasta" | "fna") {
                    FastTagMapper::from_fasta(
                        path,
                        "custom",
                        self.tag_kmer_size,
                    )?
                } else {
                    FastTagMapper::from_tsv(
                        path,
                        "custom",
                        self.tag_kmer_size,
                    )?
                }
            }

            (None, BuiltinTagSet::BdHuman) => {
                FastTagMapper::bd_human(self.tag_kmer_size)?
            }

            (None, BuiltinTagSet::BdMouse) => {
                FastTagMapper::bd_mouse(self.tag_kmer_size)?
            }
        };

        Ok(Some(
            mapper
                .with_max_mismatches(
                    self.tag_min_matches.saturating_sub(1)
                )
        ))
    }
}