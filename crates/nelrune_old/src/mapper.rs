use std::collections::HashMap;

use anyhow::{anyhow, bail, Context, Result};

use bam_tide::quantification::chunk_processor::ChunkProcessor;
use bam_tide::quantification::cli::QuantMode;
use bam_tide::quantification::job::{
    Job,
    JobBuilder,
};
use bam_tide::results::QuantData;

use sc_mapper::{
    MappingCall,
    StreamingMapper,
};

use crate::quant::CHUNK;

use rust_htslib::bam::{Record, HeaderView, record::Aux};
use gtf_splice_index::SpliceIndex;

use crate::fastq::ParsedPair;

pub struct ReadTags {
    pub cell_seq: Vec<u8>,
    pub umi_seq: Vec<u8>,
}

pub fn remember_tags(
    read: &ParsedPair,
    pending: &mut HashMap<String, ReadTags>,
) -> Result<()> {
    if pending
        .insert(
            read.read_id.clone(),
            ReadTags {
                cell_seq:
                    read.cell_seq.clone(),

                umi_seq:
                    read.umi_seq.clone(),
            },
        )
        .is_some()
    {
        bail!(
            "duplicate outstanding read id '{}'",
            read.read_id
        );
    }

    Ok(())
}

pub fn drain_mapper(
    mapper: &mut StreamingMapper,
    pending: &mut HashMap<String, ReadTags>,
    job_builder: &JobBuilder<'_>,
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    while let Some(call) =
        mapper.try_next()?
    {
        consume_mapping_call(
            call,
            pending,
            job_builder,
            processor,
            quant_mode,
            jobs,
            data,
        )?;
    }

    Ok(())
}

pub fn consume_mapping_call(
    call: MappingCall,
    pending: &mut HashMap<String, ReadTags>,
    job_builder: &JobBuilder<'_>,
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    let tags =
        pending
            .remove(&call.read_id)
            .ok_or_else(|| {
                anyhow!(
                    "mapper returned read '{}' without matching primer metadata",
                    call.read_id
                )
            })?;

    let cell =
        std::str::from_utf8(
            &tags.cell_seq,
        )
        .context(
            "cell barcode is not valid ASCII"
        )?;

    let umi =
        std::str::from_utf8(
            &tags.umi_seq,
        )
        .context(
            "UMI is not valid ASCII"
        )?;

    /*
     * sc_mapper may return multiple alignments.
     *
     * For now these are passed through the existing JobBuilder
     * semantics. JobBuilder currently ignores secondary /
     * supplementary records.
     *
     * Proper multimapper-aware quantification belongs in
     * bam_tide, not here.
     */

    for mapper_record in
        call.records.records
    {
        let mut rec =
            mapper_record.into_inner();

        set_aux_string(
            &mut rec,
            b"CB",
            cell,
        )?;

        set_aux_string(
            &mut rec,
            b"UB",
            umi,
        )?;

        if let Some(job) =
            job_builder.build(
                &rec,
                &mut data.report,
            )?
        {
            jobs.push(job);
        }
    }

    if jobs.len() >= CHUNK {
        flush_jobs(
            processor,
            quant_mode,
            jobs,
            data,
        )?;
    }

    Ok(())
}


pub fn flush_jobs(
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    if jobs.is_empty() {
        return Ok(());
    }

    println!(
        "Processing {} mapping jobs",
        jobs.len()
    );

    data.report.stop_file_io_time();

    processor.process_into(
        quant_mode,
        jobs,
        data,
    )?;

    jobs.clear();

    Ok(())
}

fn set_aux_string(
    record: &mut Record,
    tag: &[u8; 2],
    value: &str,
) -> Result<()> {
    /*
     * Mapper output should not already carry Nelrune's
     * barcode tags, but remove first so this is deterministic.
     */

    let _ =
        record.remove_aux(tag);

    record
        .push_aux(
            tag,
            Aux::String(value),
        )
        .with_context(|| {
            format!(
                "adding BAM tag {}",
                String::from_utf8_lossy(tag)
            )
        })?;

    Ok(())
}

pub fn validate_reference_compatibility(
    header: &HeaderView,
    idx: &SpliceIndex,
) -> Result<()> {
    let mut matched = 0usize;
    let mut unmatched = Vec::new();

    for tid in 0..header.target_count() {
        let chr =
            std::str::from_utf8(
                header.tid2name(tid),
            )
            .context(
                "mapper header contains an invalid UTF-8 chromosome name",
            )?;

        if idx.chr_id(chr).is_some() {
            matched += 1;
        } else {
            unmatched.push(chr.to_string());
        }
    }

    if matched == 0 {
        bail!(
            "mapper reference and splice index appear incompatible: \
             none of the {} mapper contigs could be resolved by the splice index",
            header.target_count(),
        );
    }

    eprintln!(
        "Reference check: {matched}/{} mapper contigs are present in the splice index.",
        header.target_count(),
    );

    if !unmatched.is_empty() {
        eprintln!(
            "  {} mapper contigs are not annotated by the splice index.",
            unmatched.len(),
        );
    }

    Ok(())
}
