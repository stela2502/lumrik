use anyhow::{Context, Result, bail};
use bam_tide::{alignment_evidence::read_alignment_evidence, fastq::FastqRecord};
use clap::Parser;
use reference_curator::{Annotation, Producer, ReferenceCurator, Strand};
use rust_htslib::bam;
use sc_mapper::{MappingCall, StreamingMapper, StreamingMapperCli};
use std::{fs, path::PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about = "Map unresolved reference_curator candidates and attach raw genomic alignment evidence")]
struct Args {
    #[arg(long)] store: PathBuf,
    #[arg(long)] reference_id: String,
    #[arg(long)] bam_out: Option<PathBuf>,
    #[command(flatten)] mapper: StreamingMapperCli,
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.mapper.mapper_paired { bail!("reference candidate mapping is single-end; omit --mapper-paired"); }
    let mut curator = ReferenceCurator::open(&args.store)
        .with_context(|| format!("opening curator {}", args.store.display()))?;
    let mut candidates: Vec<_> = curator.unresolved().map(|c| (c.id.clone(), c.sequence.clone())).collect();
    candidates.sort_by(|a,b| a.0.cmp(&b.0));
    if candidates.is_empty() { bail!("curator contains no unresolved candidates"); }

    let bam_path = args.bam_out.clone().unwrap_or_else(|| args.store.with_extension("mapping.bam"));
    if let Some(parent)=bam_path.parent() { fs::create_dir_all(parent)?; }
    let mut mapper = args.mapper.from_cli().context("starting mapper")?;
    let mut writer: Option<bam::Writer> = None;
    let mut calls = Vec::new();
    for (id, seq) in &candidates {
        let read = FastqRecord::new(id.clone(), seq, &vec![40u8; seq.len()]);
        mapper.submit(&read, None)?;
        if writer.is_none() && mapper.header_loaded() { writer=Some(open_writer(&bam_path, &mut mapper)?); }
        while let Some(call)=mapper.try_next()? { calls.push(call); }
    }
    if writer.is_none() { writer=Some(open_writer(&bam_path, &mut mapper)?); }
    while let Some(call)=mapper.try_next()? { calls.push(call); }
    let tail=mapper.finish()?; calls.extend(tail);
    let mut writer=writer.context("mapper produced no BAM header")?;
    for call in calls { write_call(&mut writer, call)?; }
    drop(writer);

    let evidence=read_alignment_evidence(&bam_path)?;
    let producer=Producer { name: format!("{:?}", args.mapper.mapper), version: None,
        parameters: Some(format!("index={} threads={} options={}", args.mapper.mapper_index.display(), args.mapper.mapper_threads, args.mapper.mapper_options.as_deref().unwrap_or(""))) };
    let mut added=0usize;
    for ev in evidence {
        if curator.candidate(&ev.query).is_none() { continue; }
        let mut a=Annotation::new(producer.clone(), ev.target, ev.start, ev.end,
            if ev.reverse { Strand::Reverse } else { Strand::Forward }, ev.query_len);
        a.cigar=Some(ev.cigar); a.mapq=Some(ev.mapq); a.edit_distance=ev.edit_distance;
        a.query_start=ev.query_start; a.query_end=ev.query_end; a.secondary=ev.secondary; a.supplementary=ev.supplementary;
        curator.add_annotation(&ev.query, args.reference_id.clone(), a)?; added+=1;
    }
    curator.save(&args.store)?;
    println!("Reference candidate mapping");
    println!("  unresolved candidates submitted: {}", candidates.len());
    println!("  BAM:                             {}", bam_path.display());
    println!("  mapped alignment records added: {added}");
    println!("  reference curator:              {}", args.store.display());
    Ok(())
}

fn open_writer(path:&PathBuf, mapper:&mut StreamingMapper)->Result<bam::Writer>{
    let header=mapper.header()?.clone();
    bam::Writer::from_path(path,&header,bam::Format::Bam).with_context(||format!("creating {}",path.display()))
}
fn write_call(writer:&mut bam::Writer, call:MappingCall)->Result<()> {
    for rec in call.records.records { writer.write(&rec.into_inner())?; }
    Ok(())
}
