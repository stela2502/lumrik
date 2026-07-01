use crate::cli::CoverageCli;
use rust_htslib::bam::Record;

/// SAM flag bits (u16) and nice to have
#[allow(unused)]
const FLAG_PAIRED: u16 = 0x1;
#[allow(unused)]
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
#[allow(unused)]
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const FLAG_READ1: u16 = 0x40;
#[allow(unused)]
const FLAG_READ2: u16 = 0x80; // "second in pair"

#[derive(Clone, Copy, Debug)]
pub struct AlignmentPolicy {
    pub min_mapq: u8,

    /// Effective include mask: if 0 => no include constraint.
    sam_flag_include: u16,

    /// Effective exclude mask: bits set here are rejected.
    sam_flag_exclude: u16,
}

impl AlignmentPolicy {
    /// CLI constructor:
    /// - builds computed include/exclude masks from convenience options
    /// - if `cli.sam_flag_include` / `cli.sam_flag_exclude` are provided, they override computed
    ///
    /// IMPORTANT: This matches bamCoverage help defaults: include/exclude default to None (=no filtering).
    pub fn from_cli(cli: &CoverageCli) -> Self {
        let computed_ex = Self::computed_exclude_mask_bamcoverage(cli);
        let computed_in = Self::computed_include_mask_bamcoverage(cli);

        // Superseded by explicit CLI masks (including 0)
        let effective_ex = cli.sam_flag_exclude.unwrap_or(computed_ex);
        let effective_in = cli.sam_flag_include.unwrap_or(computed_in);

        Self {
            min_mapq: cli.min_mapping_quality,
            sam_flag_include: effective_in,
            sam_flag_exclude: effective_ex,
        }
    }

    /// bamCoverage-like computed EXCLUDE mask from your “nice” booleans.
    ///
    /// Since bamCoverage defaults to --samFlagExclude None, the computed default here is 0
    /// unless you add convenience toggles like "ignore_duplicates" etc.
    fn computed_exclude_mask_bamcoverage(cli: &CoverageCli) -> u16 {
        let mut mask = 0u16;

        // If you expose these as "ignore" booleans (bamCoverage has --ignoreDuplicates),
        // you can map them here. If you only have "include_duplicates", invert it:
        if !cli.include_duplicates {
            mask |= FLAG_DUPLICATE;
        }

        // If you *want* to keep the old toggles (include_secondary/supplementary),
        // they can be mapped too. Note: bamCoverage default is not to exclude these via samFlagExclude,
        // but users sometimes do. Keeping them as convenience flags is fine.
        if !cli.include_secondary {
            mask |= FLAG_SECONDARY;
        }
        if !cli.include_supplementary {
            mask |= FLAG_SUPPLEMENTARY;
        }

        // I recommend NOT forcing QCFAIL or UNMAPPED into the mask by default.
        // Unmapped should be handled explicitly in passes_filter (see below).
        // QCFAIL: only exclude if you have an explicit toggle; otherwise leave it alone.

        mask
    }

    /// bamCoverage-like computed INCLUDE mask from your “nice” booleans.
    ///
    /// In bamCoverage, `--samFlagInclude` default is None, meaning no include constraint.
    /// We represent "None" as 0 here.
    fn computed_include_mask_bamcoverage(cli: &CoverageCli) -> u16 {
        let mut mask = 0u16;

        // Your "only_r1" convenience option corresponds to --samFlagInclude 64.
        if cli.only_r1 {
            mask |= FLAG_READ1;
        }

        mask
    }

    /// Check whether a BAM record should be included given this policy.
    /// Semantics:
    /// - Unmapped reads are always excluded (regardless of masks)
    /// - Exclude-mask: if any excluded bit is set -> reject
    /// - Include-mask: if non-zero, all bits must be present -> accept; otherwise no constraint
    /// - MAPQ threshold
    #[inline]
    pub fn passes_filter(&self, rec: &Record) -> bool {
        // Always exclude unmapped for coverage.
        if rec.is_unmapped() {
            return false;
        }

        // MAPQ
        if rec.mapq() < self.min_mapq && rec.mapq() == 255 {
            return false;
        }

        let flag: u16 = rec.flags();

        // Exclude bits
        if (flag & self.sam_flag_exclude) != 0 {
            return false;
        }

        // Include bits (0 means "no include constraint")
        if self.sam_flag_include != 0 && (flag & self.sam_flag_include) != self.sam_flag_include {
            return false;
        }

        true
    }

    /// Optional: expose these for logging/debugging
    pub fn sam_flag_exclude(&self) -> u16 {
        self.sam_flag_exclude
    }
    pub fn sam_flag_include(&self) -> u16 {
        self.sam_flag_include
    }
}
