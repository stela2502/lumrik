use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;
use int_to_dna::TwoBitReader;
use ommverse::Ommverse;
use rust_htslib::bam::{self, Read, record::Aux};
use sc_primer::{BdCellVersion, Chemistry, RhapsodyWhitelist};

const DEFAULT_CELL_ID: u64 = 38_637_011;
const DEFAULT_GEX_GENES: &[&str] = &["Actb", "Gapdh"];
const VDJ_GENES: &[&str] = &[
    "Ighv1-64", "Ighd1-1", "Ighj3", "Igha",
    "Igkv6-15", "Igkj2", "Igkc",
    "Iglv3", "Iglj2",
];
const LOCUS_SPACER: usize = 101;

#[derive(Parser, Debug)]
#[command(
    name = "nelrune materialize-norn-test-data",
    about = "Materialize the tiny real-biological Norn end-to-end fixture"
)]
struct Args {
    /// Ommverse index whose source genome supplies the fixture chromosomes.
    #[arg(long)]
    ommverse_index: PathBuf,

    /// Output directory for genome.fa, genes.gtf, R1.fastq and R2.fastq.
    #[arg(long)]
    out: PathBuf,

    /// Rich one-cell VDJ BAM whose read/UMI evidence is rehydrated as raw BD reads.
    #[arg(long)]
    vdj_bam: PathBuf,

    /// GEX genes to include. Defaults to Actb and Gapdh.
    #[arg(long, num_args = 1..)]
    gex_gene: Vec<String>,

    /// Official one-based BD/Rustody v2.384 cell id used by the fixture.
    #[arg(long, default_value_t = DEFAULT_CELL_ID)]
    cell_id: u64,
}

pub fn run() -> Result<()> {
    let args = Args::parse_from(std::env::args().skip(1));
    materialize(args)
}

fn materialize(args: Args) -> Result<()> {
    let omm = Ommverse::load(&args.ommverse_index)
        .with_context(|| format!("loading Ommverse index {}", args.ommverse_index.display()))?;
    let mut genome = TwoBitReader::open(&omm.genome_twobit)
        .with_context(|| format!("opening source genome {}", omm.genome_twobit.display()))?;

    fs::create_dir_all(&args.out)
        .with_context(|| format!("creating {}", args.out.display()))?;

    let gex_genes: Vec<String> = if args.gex_gene.is_empty() {
        DEFAULT_GEX_GENES.iter().map(|x| (*x).to_owned()).collect()
    } else {
        args.gex_gene.clone()
    };
    let mut requested = gex_genes.clone();
    requested.extend(VDJ_GENES.iter().map(|x| (*x).to_owned()));

    let fasta_path = args.out.join("genome.fa");
    let gtf_path = args.out.join("genes.gtf");
    let mut fasta = BufWriter::new(File::create(&fasta_path)?);
    let mut gtf = BufWriter::new(File::create(&gtf_path)?);

    // One compact synthetic chromosome is deliberate: the 101-N spacer keeps
    // neighboring genes outside the splice index's 100-bp gene overhang.
    let chromosome = "norn_fixture";
    writeln!(fasta, ">{chromosome}")?;
    let mut synthetic_pos = 0u32;
    let mut first = true;
    let mut gex_r2 = Vec::<(String, Vec<u8>)>::new();

    for symbol in &requested {
        let gene = omm.gene(symbol)
            .with_context(|| format!("Ommverse index does not contain required fixture gene {symbol}"))?;
        let tx_id = *gene.transcript_ids().first()
            .with_context(|| format!("required fixture gene {symbol} has no transcript"))?;
        let tx = omm.splice.transcripts.get(tx_id)
            .with_context(|| format!("missing transcript {tx_id} for {symbol}"))?;
        let exons = tx.exons();
        let start = exons.iter().map(|x| x.start).min()
            .with_context(|| format!("required fixture gene {symbol} has no exons"))?;
        let end = exons.iter().map(|x| x.end).max().unwrap();
        let source_chr = omm.splice.chr_names.get(tx.chr_id)
            .with_context(|| format!("missing source chromosome for {symbol}"))?;

        if !first {
            write!(fasta, "{}", "N".repeat(LOCUS_SPACER))?;
            synthetic_pos += LOCUS_SPACER as u32;
        }
        first = false;

        let dna = genome.sequence(source_chr, start, end)?
            .to_string((end - start) as usize);
        write!(fasta, "{dna}")?;

        let gene_start_1 = synthetic_pos + 1;
        let gene_end_1 = synthetic_pos + (end - start);
        let strand = match tx.strand {
            gtf_splice_index::Strand::Minus => "-",
            _ => "+",
        };
        let gene_id = format!("norn_{symbol}");
        let tx_name = tx.primary_name().unwrap_or(symbol);
        writeln!(gtf, "{chromosome}\tLumrik\tgene\t{gene_start_1}\t{gene_end_1}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; gene_name \"{symbol}\";")?;
        writeln!(gtf, "{chromosome}\tLumrik\ttranscript\t{gene_start_1}\t{gene_end_1}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{tx_name}\"; gene_name \"{symbol}\";")?;
        for (i, exon) in exons.iter().enumerate() {
            let s = synthetic_pos + (exon.start - start) + 1;
            let e = synthetic_pos + (exon.end - start);
            writeln!(gtf, "{chromosome}\tLumrik\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{tx_name}\"; gene_name \"{symbol}\"; exon_number \"{}\";", i + 1)?;
        }

        if gex_genes.iter().any(|x| x == symbol) {
            let exon = exons.iter().max_by_key(|x| x.len()).unwrap();
            let read_len = usize::min(75, exon.len() as usize);
            if read_len < 30 {
                bail!("GEX fixture gene {symbol} has no exon long enough for a useful synthetic read");
            }
            let seq = genome.sequence(source_chr, exon.start, exon.start + read_len as u32)?
                .to_string(read_len).into_bytes();
            gex_r2.push((symbol.clone(), seq));
        }

        synthetic_pos += end - start;
    }
    writeln!(fasta)?;
    fasta.flush()?;
    gtf.flush()?;

    let whitelist = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
    let vdj_whitelist = RhapsodyWhitelist::builtin(BdCellVersion::V2_384Vdj);
    let gex_cell_cassette = whitelist.cell_id_to_cassette(args.cell_id)
        .with_context(|| format!("BD v2.384 cell id {} is outside the whitelist", args.cell_id))?;
    let vdj_cell_cassette = vdj_whitelist.cell_id_to_cassette(args.cell_id)
        .with_context(|| format!("BD v2.384-vdj cell id {} is outside the whitelist", args.cell_id))?;
    let gex_grammar = Chemistry::BdV2_384.grammar().map_err(anyhow::Error::msg)?;
    let vdj_grammar = Chemistry::BdV2_384Vdj.grammar().map_err(anyhow::Error::msg)?;

    let mut r1 = BufWriter::new(File::create(args.out.join("R1.fastq"))?);
    let mut r2 = BufWriter::new(File::create(args.out.join("R2.fastq"))?);
    let mut ordinal = 0usize;

    // Enough distinct GEX UMIs to make the single cell explicit even when the
    // Norn test uses a fixed low cell cutoff.
    for (gene_i, (symbol, seq)) in gex_r2.iter().enumerate() {
        for copy in 0..4usize {
            let umi = six_base_umi(gene_i * 4 + copy);
            let primer = gex_grammar.synthesize(&gex_cell_cassette, &umi).map_err(anyhow::Error::msg)?;
            write_pair(&mut r1, &mut r2, ordinal, &format!("gex_{symbol}"), &primer, seq)?;
            ordinal += 1;
        }
    }

    let vdj_bam = args.vdj_bam;
    let mut bam = bam::Reader::from_path(&vdj_bam)
        .with_context(|| format!("opening one-cell VDJ fixture {}", vdj_bam.display()))?;
    let mut vdj_umis = BTreeSet::new();
    let mut vdj_reads = 0usize;
    for rec in bam.records() {
        let rec = rec?;
        if rec.is_unmapped() || rec.seq_len() < 20 {
            continue;
        }
        let umi = match rec.aux(b"UB") {
            Ok(Aux::String(x)) if !x.is_empty() => x.as_bytes(),
            _ => continue,
        };
        let umi6 = normalize_umi6(umi);
        let primer = vdj_grammar.synthesize(&vdj_cell_cassette, &umi6).map_err(anyhow::Error::msg)?;
        let seq = rec.seq().as_bytes();
        write_pair(&mut r1, &mut r2, ordinal, "vdj_fixture", &primer, &seq)?;
        ordinal += 1;
        vdj_reads += 1;
        vdj_umis.insert(String::from_utf8_lossy(&umi6).into_owned());
    }
    if vdj_reads == 0 {
        bail!("one-cell VDJ fixture yielded no mapped reads carrying UB tags");
    }
    r1.flush()?;
    r2.flush()?;

    fs::write(
        args.out.join("expected.txt"),
        format!(
            "cell_id\t{}\ncell_barcode\t{}\nvdj_reads\t{}\nvdj_distinct_synthetic_umis\t{}\nexpected_IGH\tIghv1-64/Ighd1-1/Ighj3/Igha\nexpected_IGK\tIgkv6-15/Igkj2/Igkc\nexpected_IGL\tIglv3/Iglj2\n",
            args.cell_id,
            String::from_utf8_lossy(&whitelist.cell_id_to_seq(args.cell_id).unwrap()),
            vdj_reads,
            vdj_umis.len(),
        ),
    )?;

    eprintln!("[nelrune] materialized Norn fixture in {}", args.out.display());
    eprintln!("[nelrune] source Ommverse assembly: {}", omm.assembly);
    eprintln!("[nelrune] synthetic chromosome length: {} bp", synthetic_pos);
    eprintln!("[nelrune] inter-locus spacer: {} bp", LOCUS_SPACER);
    eprintln!("[nelrune] raw read pairs: {}", ordinal);
    Ok(())
}

fn normalize_umi6(umi: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(6);
    for &b in umi.iter().take(6) {
        out.push(match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
            _ => b'A',
        });
    }
    while out.len() < 6 { out.push(b'A'); }
    out
}

fn six_base_umi(mut n: usize) -> Vec<u8> {
    let mut out = vec![b'A'; 6];
    for base in out.iter_mut().rev() {
        *base = [b'A', b'C', b'G', b'T'][n & 3];
        n >>= 2;
    }
    out
}

fn write_pair(
    r1: &mut impl Write,
    r2: &mut impl Write,
    ordinal: usize,
    label: &str,
    primer: &[u8],
    insert: &[u8],
) -> Result<()> {
    let name = format!("@norn_fixture_{ordinal}_{label}");
    writeln!(r1, "{name}/1")?;
    writeln!(r1, "{}", String::from_utf8_lossy(primer))?;
    writeln!(r1, "+")?;
    writeln!(r1, "{}", "I".repeat(primer.len()))?;
    writeln!(r2, "{name}/2")?;
    writeln!(r2, "{}", String::from_utf8_lossy(insert))?;
    writeln!(r2, "+")?;
    writeln!(r2, "{}", "I".repeat(insert.len()))?;
    Ok(())
}
