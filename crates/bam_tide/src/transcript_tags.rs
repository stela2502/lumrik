//! Assign a discovered transcript identifier to a genomic BAM alignment.
use anyhow::{Context, Result};
use rust_htslib::bam::record::{Aux, Record};

/// Set `TX` to the stable, descriptive transcript identifier. Existing `TX`
/// (e.g. the original transcriptome alignment) is replaced deliberately;
/// callers must preserve it separately if it is needed for provenance.
pub fn set_transcript_assignment(record: &mut Record, transcript_name: &str) -> Result<()> {
    anyhow::ensure!(!transcript_name.is_empty(), "empty transcript name");
    if record.aux(b"TX").is_ok() {
        record.remove_aux(b"TX").context("removing previous TX assignment")?;
    }
    record.push_aux(b"TX", Aux::String(transcript_name))
        .context("adding discovered transcript TX assignment")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn writes_and_replaces_transcript_tag() {
        let mut record = Record::new();
        set_transcript_assignment(&mut record, "TP53_205_mod_ex5_del").unwrap();
        assert!(matches!(record.aux(b"TX").unwrap(), Aux::String("TP53_205_mod_ex5_del")));
        set_transcript_assignment(&mut record, "TP53_205_mod_ex5_add_101_150").unwrap();
        assert!(matches!(record.aux(b"TX").unwrap(), Aux::String("TP53_205_mod_ex5_add_101_150")));
    }
}
