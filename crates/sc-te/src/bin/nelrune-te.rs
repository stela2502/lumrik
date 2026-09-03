use std::collections::HashSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use bam_tide::fastq::FastqRecord;
use clap::Parser;
use int_to_str::int_to_str::IntToStr;
use read_tag_table::ReadTagRecord;
use rust_htslib::bam::{self, HeaderView, Read, record::Aux};
use sc_mapper::{MapperKind, StreamingMapper, StreamingMapperCli};
use sc_te::{TeCollector, TeIndex};
use scdata::Scdata;

#[derive(Debug, Parser)]
#[command(
    name="nelrune-te",
    about="Posterior TE rescue from a queryname-sorted 10x BAM"
)]
struct Cli {
    #[arg(long)] bam: PathBuf,
    #[arg(long="te-index")] te_index: PathBuf,
    #[arg(long="mapper-index")] mapper_index: PathBuf,
    #[arg(long="mapper-bin", default_value="STAR")] mapper_bin: PathBuf,
    #[arg(long, default_value_t=4)] threads: usize,
    #[arg(long, default_value_t=100)] max_multimap: usize,
    /// Reuse complete source-BAM mapping groups up to this NH value. Higher-NH
    /// groups, incomplete groups, and unmapped reads are remapped once with STAR.
    #[arg(long, default_value_t=5)] remap_nh_above: usize,
    #[arg(long, default_value_t=100)] em_iterations: usize,
    #[arg(long, default_value_t=1e-7)] em_epsilon: f64,
    #[arg(long, default_value="CB")] cell_tag: String,
    #[arg(long, default_value="UB")] umi_tag: String,
    #[arg(long, default_value_t=16)] cell_barcode_len: usize,
    #[arg(long)] out: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.cell_tag.len()!=2 || cli.umi_tag.len()!=2 { bail!("--cell-tag and --umi-tag must be exactly two characters"); }
    fs::create_dir_all(&cli.out).with_context(|| format!("creating {}", cli.out.display()))?;
    let mut index = TeIndex::load(&cli.te_index)?;
    eprintln!(
        "[nelrune-te] loaded TE splice index: {} genes, {} transcripts/loci, bin_width={} bp",
        index.annotation_gene_count(),
        index.annotation_transcript_count(),
        index.splice_index().bin_width,
    );

    let mapper_cli = StreamingMapperCli {
        mapper: MapperKind::Star,
        mapper_bin: Some(cli.mapper_bin.clone()),
        mapper_index: cli.mapper_index.clone(),
        mapper_options: Some(format!("--outFilterMultimapNmax {} --winAnchorMultimapNmax {}", cli.max_multimap, cli.max_multimap.saturating_mul(2))),
        mapper_threads: cli.threads.max(1),
        mapper_paired: false,
        mapper_keep_multimappers: true,
    };
    let mut mapper = mapper_cli.from_cli()?;
    let mut collector = TeCollector::new(cli.threads);
    let mut bam = bam::Reader::from_path(&cli.bam).with_context(|| format!("opening {}", cli.bam.display()))?;
    let header = bam.header().to_owned();
    require_queryname_sorted(&header)?;
    let cell_tag: [u8;2] = cli.cell_tag.as_bytes().try_into().unwrap();
    let umi_tag: [u8;2] = cli.umi_tag.as_bytes().try_into().unwrap();

    let mut original_records = 0usize;
    let mut read_groups = 0usize;
    let mut submitted = 0usize;
    let mut reused = 0usize;
    let mut skipped = 0usize;
    let mut group = Vec::<bam::Record>::new();
    let mut current_qname = Vec::<u8>::new();

    for rec in bam.records() {
        let rec = rec?;
        original_records += 1;

        if group.is_empty() {
            current_qname.extend_from_slice(rec.qname());
            group.push(rec);
            continue;
        }

        if rec.qname() == current_qname.as_slice() {
            group.push(rec);
            continue;
        }

        read_groups += 1;
        match process_read_group(
            &group,
            &header,
            &mut index,
            &mut collector,
            &mut mapper,
            &cell_tag,
            &umi_tag,
            cli.remap_nh_above,
        )? {
            1 => reused += 1,
            2 => submitted += 1,
            _ => skipped += 1,
        }
        drain_ready(&mut mapper, &mut collector, &mut index)?;

        group.clear();
        current_qname.clear();
        current_qname.extend_from_slice(rec.qname());
        group.push(rec);
    }

    if !group.is_empty() {
        read_groups += 1;
        match process_read_group(
            &group,
            &header,
            &mut index,
            &mut collector,
            &mut mapper,
            &cell_tag,
            &umi_tag,
            cli.remap_nh_above,
        )? {
            1 => reused += 1,
            2 => submitted += 1,
            _ => skipped += 1,
        }
        drain_ready(&mut mapper, &mut collector, &mut index)?;
    }

    eprintln!(
        "[nelrune-te] scanned {original_records} BAM records in {read_groups} QNAME groups; reused {reused} source-BAM groups; submitted {submitted} hard groups to STAR; skipped {skipped} groups"
    );

    let mapper_header = if submitted > 0 {
        Some(HeaderView::from_header(mapper.header()?))
    } else {
        None
    };
    let remaining = mapper.finish()?;
    if let Some(mh) = mapper_header.as_ref() {
        for call in remaining {
            collector.push_cluster(call.records, mh, &mut index)?;
        }
    }

    let mut result = collector.finish(&index, cli.em_iterations, cli.em_epsilon);
    write_matrix(&cli.out.join("anchor"), &mut result.anchor, &index, cli.cell_barcode_len)?;
    write_matrix(&cli.out.join("rescued_unique"), &mut result.rescued_unique, &index, cli.cell_barcode_len)?;
    write_matrix(&cli.out.join("multi_em"), &mut result.multi_em, &index, cli.cell_barcode_len)?;
    write_matrix(&cli.out.join("multi_anchored_em"), &mut result.multi_anchored_em, &index, cli.cell_barcode_len)?;
    write_matrix(&cli.out.join("multi_unanchored_em"), &mut result.multi_unanchored_em, &index, cli.cell_barcode_len)?;
    let mut f=File::create(cli.out.join("mapping_info.txt"))?;
    write!(f, "{}", result.report)?;
    drop(f);
    let mut combined = result.combined();
    write_matrix(&cli.out.join("combined_em"), &mut combined, &index, cli.cell_barcode_len)?;
    eprintln!("[nelrune-te] wrote {}", cli.out.display());
    Ok(())
}

/// Returns 1 when the complete source-BAM mapping group was reused, 2 when one
/// representative read was submitted to STAR, and 0 when the group was skipped.
fn process_read_group(
    group: &[bam::Record],
    header: &HeaderView,
    index: &mut TeIndex,
    collector: &mut TeCollector,
    mapper: &mut StreamingMapper,
    cell_tag: &[u8; 2],
    umi_tag: &[u8; 2],
    remap_nh_above: usize,
) -> Result<u8> {
    let Some((cell, umi)) = group_tags(group, cell_tag, umi_tag) else {
        return Ok(0);
    };
    let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
    let umi_id = IntToStr::new(umi.as_bytes()).into_u64();

    let mapped: Vec<_> = group
        .iter()
        .filter(|rec| !rec.is_unmapped() && !rec.is_supplementary())
        .collect();
    let reported_nh = group
        .iter()
        .filter_map(|rec| aux_u32(rec, b"NH"))
        .filter(|&nh| nh > 0)
        .max()
        .map(|nh| nh as usize)
        .unwrap_or(mapped.len());

    // We can trust/reuse a low-NH source group only when the BAM actually contains
    // all alignments it says belong to the read. Otherwise remap once rather than
    // silently treating a truncated candidate set as complete.
    let complete_low_nh = !mapped.is_empty()
        && reported_nh <= remap_nh_above
        && mapped.len() >= reported_nh;

    if complete_low_nh {
        let mut candidates = HashSet::new();
        for rec in mapped {
            candidates.extend(index.record_overlaps(rec, header)?);
        }
        let mut candidates: Vec<_> = candidates.into_iter().collect();
        candidates.sort_unstable();
        collector.add_original_candidates(cell_id, umi_id, &candidates);
        return Ok(1);
    }

    let Some(rec) = representative_with_sequence(group) else {
        return Ok(0);
    };
    let qname = String::from_utf8_lossy(rec.qname()).to_string();
    let tag = ReadTagRecord::new(qname.clone(), None, cell.as_bytes(), [], umi.as_bytes(), []);
    let seq = rec.seq().as_bytes();
    let mut fq = FastqRecord::new(tag.extend_qname(&qname), &seq, rec.qual());

    // BAM stores SEQ/QUAL in alignment orientation. Restore the original read
    // orientation before handing a mapped high-NH read back to STAR.
    if rec.is_reverse() {
        fq = fq.revcomp();
    }

    mapper.submit(&fq, None)?;
    Ok(2)
}

fn drain_ready(
    mapper: &mut StreamingMapper,
    collector: &mut TeCollector,
    index: &mut TeIndex,
) -> Result<()> {
    while let Some(call) = mapper.try_next()? {
        let mh = HeaderView::from_header(mapper.header()?);
        collector.push_cluster(call.records, &mh, index)?;
    }
    Ok(())
}

fn group_tags(
    group: &[bam::Record],
    cell_tag: &[u8; 2],
    umi_tag: &[u8; 2],
) -> Option<(String, String)> {
    for rec in group {
        let Some(cell_raw) = aux_string(rec, cell_tag) else {
            continue;
        };
        let Some(umi) = aux_string(rec, umi_tag) else {
            continue;
        };
        let cell = cell_raw
            .split_once('-')
            .map_or(cell_raw.as_str(), |(barcode, _)| barcode)
            .to_string();
        return Some((cell, umi));
    }
    None
}

fn representative_with_sequence(group: &[bam::Record]) -> Option<&bam::Record> {
    group
        .iter()
        .find(|rec| !rec.is_secondary() && !rec.is_supplementary() && rec.seq_len() > 0)
        .or_else(|| group.iter().find(|rec| rec.seq_len() > 0))
}

fn aux_string(rec: &bam::Record, tag: &[u8;2]) -> Option<String> {
    match rec.aux(tag).ok()? { Aux::String(s) => Some(s.to_string()), _ => None }
}

fn aux_u32(rec: &bam::Record, tag: &[u8; 2]) -> Option<u32> {
    match rec.aux(tag).ok()? {
        Aux::U8(x) => Some(x as u32),
        Aux::I8(x) if x >= 0 => Some(x as u32),
        Aux::U16(x) => Some(x as u32),
        Aux::I16(x) if x >= 0 => Some(x as u32),
        Aux::U32(x) => Some(x),
        Aux::I32(x) if x >= 0 => Some(x as u32),
        _ => None,
    }
}

fn require_queryname_sorted(header: &bam::HeaderView) -> Result<()> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let queryname = text.lines().any(|l| {
        l.starts_with("@HD") && l.split('\t').any(|x| x == "SO:queryname")
    });
    if !queryname {
        bail!(
            "nelrune-te requires a queryname-sorted BAM (@HD SO:queryname); use `samtools sort -n` first"
        );
    }
    Ok(())
}

fn write_matrix(out: &Path, data: &mut Scdata, index: &TeIndex, barcode_len: usize) -> Result<()> {
    data.finalize_for_export(0, index);
    data.write_sparse_with_cell_len(&out.to_path_buf(), index, barcode_len)
        .map_err(anyhow::Error::msg)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    fn record(qname: &[u8], flags: u16, nh: Option<u32>) -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set(qname, Some(&CigarString(vec![Cigar::Match(10)])), b"AAAAAAAAAA", &[30; 10]);
        rec.set_flags(flags);
        if let Some(nh) = nh {
            rec.push_aux(b"NH", Aux::U32(nh)).unwrap();
        }
        rec
    }

    #[test]
    fn low_nh_group_is_complete_only_when_all_reported_alignments_are_present() {
        let group = vec![record(b"A", 0, Some(3)), record(b"A", 0x100, Some(3))];
        let mapped = group
            .iter()
            .filter(|rec| !rec.is_unmapped() && !rec.is_supplementary())
            .count();
        let nh = group
            .iter()
            .filter_map(|rec| aux_u32(rec, b"NH"))
            .max()
            .unwrap() as usize;
        assert_eq!(mapped, 2);
        assert_eq!(nh, 3);
        assert!(mapped < nh);
    }

    #[test]
    fn representative_prefers_primary_sequence() {
        let secondary = record(b"A", 0x100, Some(2));
        let primary = record(b"A", 0, Some(2));
        let group = vec![secondary, primary];
        let chosen = representative_with_sequence(&group).unwrap();
        assert!(!chosen.is_secondary());
        assert!(!chosen.is_supplementary());
    }
}
