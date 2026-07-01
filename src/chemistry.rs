use crate::error::{PrimerError, PrimerResult};
use crate::grammar::Grammar;

use clap::ValueEnum;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum)]
pub enum Chemistry {
    /// 10x Genomics Chromium Single Cell 3' v1.
    ///
    /// Layout: adapter + 16 bp cell barcode + 10 bp UMI + polyT + insert.
    TenxThreePrimeV1,

    /// 10x Genomics Chromium Single Cell 3' v2.
    ///
    /// Layout: adapter + 16 bp cell barcode + 10 bp UMI + polyT + insert.
    TenxThreePrimeV2,

    /// 10x Genomics Chromium Single Cell 3' v3 / v3.1.
    ///
    /// Layout: adapter + 16 bp cell barcode + 12 bp UMI + polyT + insert.
    #[default]
    TenxThreePrimeV3,

    /// 10x Genomics Chromium Single Cell 3' v4 / GEM-X style preset.
    ///
    /// Layout: adapter + 16 bp cell barcode + 12 bp UMI + polyT + insert.
    TenxThreePrimeV4,

    /// 10x Genomics Chromium Single Cell 5'.
    ///
    /// Layout: adapter + 16 bp cell barcode + 10 bp UMI + insert.
    TenxFivePrime,

    /// 10x Genomics Chromium Single Cell Multiome ARC v1.
    ///
    /// Layout: adapter + 16 bp cell barcode + 12 bp UMI + insert.
    TenxMultiomeArcV1,

    /// BD Rhapsody v1 / older layout.
    BdV1,

    /// BD Rhapsody v2 96-cell combinatorial barcode layout.
    BdV2_96,

    /// BD Rhapsody v2 384-cell combinatorial barcode layout.
    BdV2_384,
}

impl Chemistry {
    pub fn parse(raw: &str) -> PrimerResult<Self> {
        match raw {
            "tenx-3p-v1" | "10x-3p-v1" | "chromium-single-cell-3-v1" => Ok(Self::TenxThreePrimeV1),
            "tenx-3p-v2" | "10x-3p-v2" | "tenx-v2" => Ok(Self::TenxThreePrimeV2),
            "tenx-3p-v3" | "10x-3p-v3" | "tenx-v3" => Ok(Self::TenxThreePrimeV3),
            "tenx-3p-v4" | "10x-3p-v4" | "tenx-v4" => Ok(Self::TenxThreePrimeV4),
            "tenx-5p" | "10x-5p" | "tenx-five-prime" => Ok(Self::TenxFivePrime),
            "tenx-multiome-arc-v1" | "10x-multiome-arc-v1" | "tenx-arc-v1" => {
                Ok(Self::TenxMultiomeArcV1)
            }

            "bd-v1" => Ok(Self::BdV1),
            "bd-v2-96" => Ok(Self::BdV2_96),
            "bd-v2-384" => Ok(Self::BdV2_384),

            other => Err(PrimerError::unknown_chemistry(other)),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::TenxThreePrimeV1 => "tenx-3p-v1",
            Self::TenxThreePrimeV2 => "tenx-3p-v2",
            Self::TenxThreePrimeV3 => "tenx-3p-v3",
            Self::TenxThreePrimeV4 => "tenx-3p-v4",
            Self::TenxFivePrime => "tenx-5p",
            Self::TenxMultiomeArcV1 => "tenx-multiome-arc-v1",

            Self::BdV1 => "bd-v1",
            Self::BdV2_96 => "bd-v2-96",
            Self::BdV2_384 => "bd-v2-384",
        }
    }

    pub fn grammar(self) -> PrimerResult<Grammar> {
        match self {
            Self::TenxThreePrimeV1 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:10+POLYT:min=0+INSERT",
            ),

            Self::TenxThreePrimeV2 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:3p-v2+UMI:10+POLYT:min=0+INSERT",
            ),

            Self::TenxThreePrimeV3 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:3p-v3+UMI:12+POLYT:min=0+INSERT",
            ),

            Self::TenxThreePrimeV4 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:3p-v4+UMI:12+POLYT:min=0+INSERT",
            ),

            Self::TenxFivePrime => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:5p+UMI:10+INSERT",
            ),

            Self::TenxMultiomeArcV1 => Grammar::parse(
                self.name(),
                "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:arc-v1+UMI:12+INSERT",
            ),

            Self::BdV1 => Grammar::parse(self.name(), "BD_CELL:v1+POLYT:min=0+INSERT"),

            Self::BdV2_96 => Grammar::parse(
                self.name(),
                "FIXED:ATAGGAAACTCATGGT:mm=2+BD_CELL:v2.96+POLYT:min=0+INSERT",
            ),

            Self::BdV2_384 => Grammar::parse(
                self.name(),
                "FIXED:ATAGGAAACTCATGGT:mm=2+BD_CELL:v2.384+POLYT:min=0+INSERT",
            ),
        }
    }
}
