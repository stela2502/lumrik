//! Stable, human-readable names for modifications of an annotated transcript.
//! Exon numbers refer to the chosen reference transcript in transcriptional order.
use anyhow::{Result, bail};
use crate::{RefBlock, Strand, Transcript};

/// A change relative to a specific annotated transcript. Genomic positions are
/// zero-based, half-open, just like `RefBlock` and `Transcript::exons()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExonModification {
    Delete { exon: usize },
    Add { after_exon: usize, interval: RefBlock },
    AlternativeStart { exon: usize, genomic_start: u32 },
    AlternativeEnd { exon: usize, genomic_end: u32 },
}

/// Construct a descriptive identifier. Changes are sorted so input order does
/// not affect the result. The caller supplies the *reference* transcript name,
/// not the newly constructed transcript's name.
pub fn modified_transcript_name(reference_name: &str, changes: &[ExonModification]) -> Result<String> {
    if reference_name.is_empty() || changes.is_empty() {
        bail!("reference transcript name and modifications must not be empty");
    }
    if !reference_name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.' || b == b'-') {
        bail!("reference transcript name contains unsupported characters");
    }
    let mut changes = changes.to_vec();
    changes.sort_by_key(|change| match *change {
        ExonModification::Delete { exon } => (exon, 0, 0, 0),
        ExonModification::Add { after_exon, interval } => (after_exon, 1, interval.start, interval.end),
        ExonModification::AlternativeStart { exon, genomic_start } => (exon, 2, genomic_start, 0),
        ExonModification::AlternativeEnd { exon, genomic_end } => (exon, 3, genomic_end, 0),
    });
    changes.dedup();
    let suffixes: Vec<String> = changes.into_iter().map(|change| match change {
        ExonModification::Delete { exon } => format!("ex{exon}_del"),
        ExonModification::Add { after_exon, interval } => format!("ex{after_exon}_add_{}_{}", interval.start + 1, interval.end),
        ExonModification::AlternativeStart { exon, genomic_start } => format!("ex{exon}_start_{}", genomic_start + 1),
        ExonModification::AlternativeEnd { exon, genomic_end } => format!("ex{exon}_end_{genomic_end}"),
    }).collect();
    Ok(format!("{reference_name}_mod_{}", suffixes.join("_")))
}

/// Identify a single skipped exon relative to the reference exon chain.
/// Requires the *complete* exon structures, not a partial read alignment.
/// Returns None when the structure cannot be described as exactly one skip.
pub fn single_exon_deletion(reference: &Transcript, candidate: &Transcript) -> Option<ExonModification> {
    if reference.chr_id != candidate.chr_id || reference.strand != candidate.strand || reference.gene_id != candidate.gene_id {
        return None;
    }
    let reference_exons = reference.exons();
    let candidate_exons = candidate.exons();
    if reference_exons.len() != candidate_exons.len() + 1 { return None; }
    for idx in 0..reference_exons.len() {
        if reference_exons[..idx] == candidate_exons[..idx] && reference_exons[idx+1..] == candidate_exons[idx..] {
            let exon = if reference.strand == Strand::Minus { reference_exons.len() - idx } else { idx + 1 };
            return Some(ExonModification::Delete { exon });
        }
    }
    None
}


/// Genomic annotation for a transcript-relative exon modification.
/// Coordinates are 1-based inclusive for human-readable output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicExonEvent {
    pub transcript_exon: usize,
    pub gene_exon: usize,
    pub genomic_start: u32,
    pub genomic_end: u32,
}

/// Gene-wide exon regions are connected components of overlapping exon intervals.
/// Adjacent but non-overlapping exons remain separate. Numbers follow transcription
/// order (descending genomic coordinates on the minus strand).
/// Exact genomic boundaries must be retained because alternative exon boundaries
/// can share the same gene-wide region number.
pub fn gene_exon_regions(transcripts: &[&Transcript]) -> Result<Vec<RefBlock>> {
    let Some(first) = transcripts.first() else { return Ok(Vec::new()); };
    if transcripts.iter().any(|tx| tx.gene_id != first.gene_id || tx.chr_id != first.chr_id || tx.strand != first.strand) {
        bail!("gene exon catalog requires transcripts from one gene, chromosome and strand");
    }
    if first.strand == Strand::Unknown {
        bail!("cannot number gene-wide exons without a known strand");
    }
    let mut blocks: Vec<RefBlock> = transcripts.iter().flat_map(|tx| tx.exons().iter().copied()).collect();
    blocks.sort_by_key(|block| (block.start, block.end));
    let mut regions: Vec<RefBlock> = Vec::new();
    for block in blocks {
        if let Some(last) = regions.last_mut() {
            if block.start < last.end {
                last.end = last.end.max(block.end);
                continue;
            }
        }
        regions.push(block);
    }
    if first.strand == Strand::Minus { regions.reverse(); }
    Ok(regions)
}

/// Resolve a transcript-relative exon to its exact genomic interval and
/// gene-wide region number. The catalog must be constructed from all annotated
/// transcripts of the same gene, not just this reference transcript.
pub fn resolve_gene_exon(reference: &Transcript, transcript_exon: usize, regions: &[RefBlock]) -> Result<GenomicExonEvent> {
    let exons = reference.exons();
    if transcript_exon == 0 || transcript_exon > exons.len() {
        bail!("transcript exon {} outside 1..={}", transcript_exon, exons.len());
    }
    let index = if reference.strand == Strand::Minus { exons.len() - transcript_exon } else { transcript_exon - 1 };
    let exon = exons[index];
    let gene_exon = regions.iter().position(|region| region.contains(exon))
        .ok_or_else(|| anyhow::anyhow!("exon interval not found in gene-wide catalog"))? + 1;
    Ok(GenomicExonEvent {
        transcript_exon, gene_exon,
        genomic_start: exon.start + 1,
        genomic_end: exon.end,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn descriptive_names_are_stable() {
        assert_eq!(modified_transcript_name("TP53_205", &[ExonModification::Delete { exon: 5 }]).unwrap(), "TP53_205_mod_ex5_del");
        let a = ExonModification::Delete { exon: 5 };
        let b = ExonModification::Add { after_exon: 6, interval: RefBlock::new(100, 150) };
        assert_eq!(modified_transcript_name("TP53_205", &[a,b]).unwrap(), modified_transcript_name("TP53_205", &[b,a,a]).unwrap());
        assert_ne!(modified_transcript_name("TP53_205", &[b]).unwrap(), modified_transcript_name("TP53_205", &[ExonModification::Add { after_exon: 6, interval: RefBlock::new(101, 150) }]).unwrap());
    }
    #[test]
    fn detect_exon_skip_on_both_strands() {
        for (strand, expected) in [(Strand::Plus, 2), (Strand::Minus, 2)] {
            let mut reference = Transcript::new(0, 1, "TP53_205", 17, strand);
            let mut candidate = Transcript::new(1, 1, "candidate", 17, strand);
            for (i, block) in [RefBlock::new(10,20), RefBlock::new(30,40), RefBlock::new(50,60)].into_iter().enumerate() {
                reference.add_exon(block);
                if i != 1 { candidate.add_exon(block); }
            }
            reference.finalize(); candidate.finalize();
            assert_eq!(single_exon_deletion(&reference, &candidate), Some(ExonModification::Delete { exon: expected }));
        }
    }
    #[test]
    fn partial_chain_is_not_exon_skip() {
        let mut reference = Transcript::new(0, 1, "ref", 17, Strand::Minus);
        let mut candidate = Transcript::new(1, 1, "candidate", 17, Strand::Minus);
        for block in [RefBlock::new(10,20), RefBlock::new(30,40), RefBlock::new(50,60)] { reference.add_exon(block); }
        candidate.add_exon(RefBlock::new(10,20));
        reference.finalize(); candidate.finalize();
        assert_eq!(single_exon_deletion(&reference, &candidate), None);
    }
}

#[cfg(test)]
mod gene_exon_tests {
    use super::*;
    #[test]
    fn minus_strand_isoforms_share_region_but_keep_exact_coordinates() {
        let mut a = Transcript::new(0, 1, "a", 17, Strand::Minus);
        let mut b = Transcript::new(1, 1, "b", 17, Strand::Minus);
        for block in [RefBlock::new(100, 120), RefBlock::new(200, 220)] { a.add_exon(block); }
        for block in [RefBlock::new(105, 120), RefBlock::new(200, 225)] { b.add_exon(block); }
        a.finalize(); b.finalize();
        let regions = gene_exon_regions(&[&a, &b]).unwrap();
        assert_eq!(regions, vec![RefBlock::new(200,225), RefBlock::new(100,120)]);
        let a_first = resolve_gene_exon(&a, 1, &regions).unwrap();
        let b_first = resolve_gene_exon(&b, 1, &regions).unwrap();
        assert_eq!(a_first.gene_exon, 1);
        assert_eq!(b_first.gene_exon, 1);
        assert_eq!((a_first.genomic_start, a_first.genomic_end), (201,220));
        assert_eq!((b_first.genomic_start, b_first.genomic_end), (201,225));
    }
}
