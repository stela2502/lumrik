use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use bam_tide::fastq::FastqRecord;
use clap::Parser;
use int_to_str::int_to_str::IntToStr;
use read_tag_table::ReadTagRecord;
use rust_htslib::bam::{self, HeaderView, Read, record::Aux};
use sc_mapper::{MapperKind, StreamingMapperCli};
use sc_te::{TeCollector, TeIndex};
use scdata::Scdata;

#[derive(Debug, Parser)]
#[command(name="nelrune-te", about="Posterior TE rescue from a coordinate-sorted 10x BAM")]
struct Cli {
    #[arg(long)] bam: PathBuf,
    #[arg(long="te-index")] te_index: PathBuf,
    #[arg(long="mapper-index")] mapper_index: PathBuf,
    #[arg(long="mapper-bin", default_value="STAR")] mapper_bin: PathBuf,
    #[arg(long, default_value_t=4)] threads: usize,
    #[arg(long, default_value_t=100)] max_multimap: usize,
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
    require_coordinate_sorted(&header)?;
    let cell_tag: [u8;2] = cli.cell_tag.as_bytes().try_into().unwrap();
    let umi_tag: [u8;2] = cli.umi_tag.as_bytes().try_into().unwrap();

    let mut original_records=0usize;
    let mut submitted=0usize;
    for rec in bam.records() {
        let rec = rec?;
        original_records += 1;
        if rec.is_secondary() || rec.is_supplementary() || rec.is_duplicate() { continue; }
        let Some(cell_raw) = aux_string(&rec, &cell_tag) else { continue; };
        let cell = cell_raw.split_once('-').map_or(cell_raw.as_str(), |(barcode, _)| barcode);
        let Some(umi) = aux_string(&rec, &umi_tag) else { continue; };
        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        let umi_id = IntToStr::new(umi.as_bytes()).into_u64();

        if !rec.is_unmapped() {
            let overlaps = index.record_overlaps(&rec, &header)?;
            if !overlaps.is_empty() { collector.add_anchor(cell_id, umi_id, &overlaps); }
            continue;
        }
        if rec.seq_len()==0 { continue; }
        let qname = String::from_utf8_lossy(rec.qname()).to_string();
        let tag = ReadTagRecord::new(qname.clone(), None, cell.as_bytes(), [], umi.as_bytes(), []);
        let fq = FastqRecord::new(tag.extend_qname(&qname), &rec.seq().as_bytes(), rec.qual());
        mapper.submit(&fq, None)?;
        submitted += 1;
        while let Some(call) = mapper.try_next()? {
            let mh = HeaderView::from_header(mapper.header()?);
            collector.push_cluster(call.records, &mh, &mut index)?;
        }
    }
    eprintln!("[nelrune-te] scanned {original_records} BAM records; submitted {submitted} unmapped reads to STAR");
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

fn aux_string(rec: &bam::Record, tag: &[u8;2]) -> Option<String> {
    match rec.aux(tag).ok()? { Aux::String(s) => Some(s.to_string()), _ => None }
}

fn require_coordinate_sorted(header: &bam::HeaderView) -> Result<()> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let coordinate = text.lines().any(|l| l.starts_with("@HD") && l.split('\t').any(|x| x=="SO:coordinate"));
    if !coordinate { bail!("nelrune-te requires a coordinate-sorted BAM (@HD SO:coordinate)"); }
    Ok(())
}

fn write_matrix(out: &Path, data: &mut Scdata, index: &TeIndex, barcode_len: usize) -> Result<()> {
    data.finalize_for_export(0, index);
    data.write_sparse_with_cell_len(&out.to_path_buf(), index, barcode_len)
        .map_err(anyhow::Error::msg)?;
    Ok(())
}
