//! Reference-backed, deterministic ONT fixture generation and independent FASTQ round-trip.
use std::{fs::{self, File}, io::{Read, Write, BufWriter}, path::Path, time::Instant};
use anyhow::{Result, Context, ensure};

use flate2::read::GzDecoder;
use serde::Deserialize;
use sc_primer::{Chemistry, Grammar};
use gtf_splice_index::{SpliceIndex, IdNameKeys, Transcript, Strand};

#[derive(Deserialize)]
struct Experiment { reference: Reference, cell_barcode: String, chemistry: Option<String>, primer_structure: Option<String>, reads_per_umi: usize, pcr_duplicates: Vec<PcrDuplicate>, seed: u64, haplotypes: Vec<Haplotype> }
#[derive(Deserialize)]
struct PcrDuplicate { haplotype: String, umi: String, fragment: usize }
#[derive(Deserialize)]
struct Reference { fasta: String, gtf: String, vcf: String, gene_id: String }
#[derive(Deserialize)]
struct Haplotype { name: String, allele: usize, umis: Vec<String>, exons: Vec<usize> }
#[derive(Clone)]
struct Variant { pos: usize, reference: u8, alternate: u8, allele: usize, id: String }
struct Fixture { exp: Experiment, chrom: String, genome: Vec<u8>, exons: Vec<(usize,usize)>, variants: Vec<Variant> }
#[derive(Clone)]
struct ReadTruth { name: String, sequence: Vec<u8>, umi: String, haplotype: String, exon: String, junction: String }

fn text_file(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("opening {}",path.display()))?;
    let mut bytes = Vec::new();
    if path.extension().is_some_and(|e| e == "gz") { GzDecoder::new(file).read_to_end(&mut bytes)?; }
    else { std::io::BufReader::new(file).read_to_end(&mut bytes)?; }
    Ok(String::from_utf8(bytes)?)
}
fn load(path: &Path) -> Result<Fixture> {
    let exp: Experiment = serde_yaml::from_str(&text_file(path)?)?;
    let root = path.parent().unwrap_or(Path::new("."));
    let fasta = text_file(&root.join(&exp.reference.fasta))?;
    let mut lines = fasta.lines();
    let chrom = lines.next().context("empty FASTA")?.strip_prefix('>').context("FASTA header missing")?.split_whitespace().next().context("empty FASTA sequence identifier")?.to_owned();
    let genome: Vec<u8> = lines.filter(|l| !l.starts_with('>')).flat_map(|l| l.as_bytes().to_vec()).map(|b|b.to_ascii_uppercase()).collect();
    ensure!(!genome.is_empty(), "empty reference genome");
    // Reuse the production GTF splice index instead of parsing exon attributes here.
    let index = SpliceIndex::from_path(root.join(&exp.reference.gtf), 1000, IdNameKeys::default())?;
    let transcript = index.transcripts.iter().find(|tx| {
        index.genes.get(tx.gene_id).is_some_and(|g| g.names.iter().any(|name| name == &exp.reference.gene_id))
            && index.chr_names.get(tx.chr_id).is_some_and(|name| name == &chrom)
    }).context("GTF splice index contains no transcript for requested gene and chromosome")?;
    ensure!(transcript.strand == Strand::Plus, "fixture requires a positive-strand transcript");
    let blocks = transcript.exons();
    ensure!(blocks.len() == 4, "this regression fixture requires exactly four exons");
    let exons: Vec<(usize, usize)> = blocks.iter().map(|b| (b.start as usize, b.end as usize)).collect();
    for pair in exons.windows(2) {
        let (left,right) = (pair[0].1,pair[1].0);
        ensure!(right >= left+4, "intron too short");
        ensure!(&genome[left..left+2] == b"GT" && &genome[right-2..right] == b"AG", "noncanonical GT-AG splice motif");
    }
    let vcf = text_file(&root.join(&exp.reference.vcf))?;
    let mut variants = Vec::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<_> = line.split('\t').collect();
        ensure!(f.len() >= 10 && f[0] == chrom, "invalid VCF row");
        let pos = f[1].parse::<usize>()? - 1;
        ensure!(pos < genome.len() && f[3].len() == 1 && f[4].len() == 1, "only SNVs supported");
        ensure!(genome[pos] == f[3].as_bytes()[0], "VCF REF mismatch at {}",pos+1);
        let gt: Vec<_> = f[9].split(':').next().unwrap_or("").split('|').collect();
        ensure!(gt.len()==2 && gt.iter().filter(|x| **x=="1").count()==1, "VCF requires phased heterozygous 1|0 or 0|1");
        variants.push(Variant { pos, reference: f[3].as_bytes()[0], alternate:f[4].as_bytes()[0], allele: if gt[0]=="1" {0} else {1}, id:f[2].to_owned() });
    }
    ensure!(variants.len()==3, "expected causal donor + two reporter variants");
    ensure!(variants[0].pos==exons[1].1 && variants[0].reference==b'G' && variants[0].alternate==b'A' && variants[0].allele==0, "expected allele-0 G>A exon2 donor mutation");
    ensure!(variants[1].pos>=exons[2].0 && variants[1].pos<exons[2].1 && variants[1].allele==0, "expected exon3 allele-0 reporter");
    ensure!(variants[2].pos>=exons[3].0 && variants[2].pos<exons[3].1 && variants[2].allele==1, "expected exon4 allele-1 reporter");
    ensure!(exp.cell_barcode.len()==16 && exp.cell_barcode.bytes().all(|b|b"ACGT".contains(&b)), "cell barcode must be 16bp DNA");
    ensure!(exp.reads_per_umi==10 && exp.seed==42, "fixture currently requires 10 reads per UMI and seed 42");
    ensure!(exp.haplotypes.len()==2, "expected two haplotypes");
    for h in &exp.haplotypes {
        ensure!(h.allele<2 && h.umis.len()==2, "expected two UMIs per allele");
        ensure!(h.exons==if h.allele==0 {vec![1,3,4]} else {vec![1,2,3,4]}, "unexpected splice structure");
        for umi in &h.umis { ensure!(umi.len()==12 && umi.bytes().all(|b|b"ACGT".contains(&b)), "UMI must be 12bp DNA"); }
    }
    Ok(Fixture {exp,chrom,genome,exons,variants})
}
fn expected(f: &Fixture, grammar: &Grammar) -> Result<Vec<ReadTruth>> {
    let total = Instant::now();
    let mut reads = Vec::with_capacity(43);
    for hap in &f.exp.haplotypes {
        let hap_started = Instant::now();
        // Apply phased variants to the complete genomic haplotype first.
        // The intronic donor mutation remains genomic and vanishes from mature RNA.
        let mut hap_genome = f.genome.clone();
        for var in &f.variants {
            if var.allele == hap.allele { hap_genome[var.pos] = var.alternate; }
        }
        eprintln!("[nelrune] fixture inspect: allele={} genome+variants={:?}", hap.allele, hap_started.elapsed());
        let tx_started = Instant::now();
        // Use the splice-index transcript model to construct the requested isoform.
        let mut tx = Transcript::new(0, 0, "synthetic_isoform", 0, Strand::Plus);
        for &n in &hap.exons {
            let (start, end) = f.exons[n - 1];
            tx.add_exon(gtf_splice_index::RefBlock::new(start as u32, end as u32));
        }
        tx.finalize();
        // Construct and decode the spliced transcript once per haplotype.
        // Previously cdna() was called again for every exon, repeatedly
        // rebuilding and packing the same complete transcript.
        eprintln!("[nelrune] fixture inspect: allele={} transcript_model={:?}", hap.allele, tx_started.elapsed());
        let cdna_started = Instant::now();
        let cdna = tx.cdna(&hap_genome).map_err(anyhow::Error::msg)?;
        eprintln!("[nelrune] fixture inspect: allele={} cdna={:?} bases={}", hap.allele, cdna_started.elapsed(), cdna.size);
        let decode_started = Instant::now();
        let cdna_bases = cdna.to_string(cdna.size);
        eprintln!("[nelrune] fixture inspect: allele={} decode={:?}", hap.allele, decode_started.elapsed());
        let fragments_started = Instant::now();
        let mut pieces = Vec::with_capacity(hap.exons.len());
        let mut offset = 0;
        for &n in &hap.exons {
            let (start, end) = f.exons[n - 1];
            let len = end - start;
            let end_offset = offset + len;
            let bases = cdna_bases.as_bytes().get(offset..end_offset)
                .context("spliced transcript shorter than its annotated exons")?;
            pieces.push((n, bases.to_vec()));
            offset = end_offset;
        }
        ensure!(offset == cdna_bases.len(), "spliced transcript length differs from exon lengths");
        // Every exon is covered in full, and each junction has a dedicated 100bp observation.
        // Additional reads shift junction windows by one base, without spanning another junction.
        let mut fragments: Vec<(String,String,Vec<u8>)> = Vec::with_capacity(f.exp.reads_per_umi);
        for (n,seq) in &pieces { fragments.push((n.to_string(),".".to_owned(),seq.clone())); }
        for (i,pair) in pieces.windows(2).enumerate() {
            let a=&pair[0].1; let b=&pair[1].1;
            for shift in 0..(if pieces.len()==3 { if i==0 {4} else {3} } else {2}) {
                let left=50+shift;
                let mut seq=a[a.len()-left..].to_vec();
                seq.extend_from_slice(&b[..50-shift]);
                fragments.push((format!("{},{}",pair[0].0,pair[1].0),format!("{}>{}",pair[0].0,pair[1].0),seq));
            }
        }
        ensure!(fragments.len()==f.exp.reads_per_umi, "fragment design must cover exactly {} observations", f.exp.reads_per_umi);
        ensure!(fragments.len()>=f.exp.reads_per_umi, "insufficient fragments");
        eprintln!("[nelrune] fixture inspect: allele={} fragment_design={:?}", hap.allele, fragments_started.elapsed());
        for umi in &hap.umis {
            let primer_started = Instant::now();
            let prefix=grammar.synthesize(f.exp.cell_barcode.as_bytes(),umi.as_bytes())
                .map_err(anyhow::Error::msg)?;
            eprintln!("[nelrune] fixture inspect: allele={} umi={} primer_synthesis={:?}", hap.allele, umi, primer_started.elapsed());
            let reads_started = Instant::now();
            for (i,(exon,junction,fragment)) in fragments.iter().take(f.exp.reads_per_umi).enumerate() {
                let mut sequence=prefix.clone(); sequence.extend_from_slice(fragment);
                reads.push(ReadTruth {name:format!("ont_{}_{}_{}",hap.name,umi,i),sequence,umi:umi.clone(),haplotype:hap.name.clone(),exon:exon.clone(),junction:junction.clone()});
            }
            eprintln!("[nelrune] fixture inspect: allele={} umi={} read_assembly={:?}", hap.allele, umi, reads_started.elapsed());
        }
    }
    ensure!(reads.len()==40,"expected exactly forty original observations");
    ensure!(f.exp.pcr_duplicates.len()==3,"fixture requires exactly three PCR duplicates");
    let duplicates_started = Instant::now();
    let originals = reads.clone();
    let mut selected = std::collections::HashSet::new();
    for (copy_idx, spec) in f.exp.pcr_duplicates.iter().enumerate() {
        ensure!(spec.fragment < f.exp.reads_per_umi, "PCR duplicate fragment out of range");
        ensure!(selected.insert((spec.haplotype.as_str(), spec.umi.as_str(), spec.fragment)),
            "PCR duplicate target specified more than once");
        let original = originals.iter().find(|r|
            r.haplotype == spec.haplotype && r.umi == spec.umi &&
            r.name.ends_with(&format!("_{}", spec.fragment)))
            .context("PCR duplicate references a nonexistent original fragment")?;
        let mut duplicate = original.clone();
        duplicate.name = format!("{}_pcr{}", original.name, copy_idx + 1);
        reads.push(duplicate);
    }
    ensure!(reads.len()==43,"expected 40 observations and three PCR duplicates");
    eprintln!("[nelrune] fixture inspect: pcr_duplicates={:?} total_construct={:?}", duplicates_started.elapsed(), total.elapsed());
    Ok(reads)
}
fn verify(path:&Path, reads:&[ReadTruth]) -> Result<()> {
    let fastq=text_file(path)?;
    let mut lines=fastq.lines();
    for read in reads {
        ensure!(lines.next()==Some(format!("@{}",read.name).as_str()),"read name mismatch: {}",read.name);
        ensure!(lines.next()==Some(std::str::from_utf8(&read.sequence)?),"sequence mismatch: {}",read.name);
        ensure!(lines.next()==Some("+"),"FASTQ separator mismatch");
        ensure!(lines.next().is_some_and(|q|q.len()==read.sequence.len() && q.bytes().all(|b|b==b'I')),"FASTQ quality mismatch");
    }
    ensure!(lines.next().is_none(),"unexpected extra FASTQ records");
    Ok(())
}
fn grammar_for(fixture: &Fixture) -> Result<Grammar> {
    let started = Instant::now();
    let grammar = if let Some(structure) = fixture.exp.primer_structure.as_deref() {
        Grammar::parse("synthetic-ont-fixture", structure).map_err(anyhow::Error::msg)?
    } else {
        let chemistry = Chemistry::parse(fixture.exp.chemistry.as_deref().unwrap_or("tenx-3p-v3"))
            .map_err(anyhow::Error::msg)?;
        chemistry.grammar().map_err(anyhow::Error::msg)?
    };
    eprintln!("[nelrune] fixture inspect: grammar initialization={:?}", started.elapsed());
    Ok(grammar)
}

fn execute_with_grammar(config:&Path,out:&Path,verify_only:bool,grammar:&Grammar)->Result<()> {
    let started = Instant::now();
    let fixture=load(config)?;
    let loaded = Instant::now();
    eprintln!("[nelrune] fixture: reference load {:?}", loaded.duration_since(started));
    let reads=expected(&fixture, grammar)?;
    let synthesized = Instant::now();
    let fq=out.join("ont_reads.fastq");
    if !verify_only {
        fs::create_dir_all(out)?;
        let mut fastq=BufWriter::new(File::create(&fq)?);
        let mut truth=BufWriter::new(File::create(out.join("expected_reads.tsv"))?);
        writeln!(truth,"read_id\tcell_barcode\tumi\thaplotype\texons\tjunction")?;
        for r in &reads {
            writeln!(fastq,"@{}\n{}\n+\n{}",r.name,std::str::from_utf8(&r.sequence)?,"I".repeat(r.sequence.len()))?;
            writeln!(truth,"{}\t{}\t{}\t{}\t{}\t{}",r.name,fixture.exp.cell_barcode,r.umi,r.haplotype,r.exon,r.junction)?;
        }
        fastq.flush()?;truth.flush()?;
        if config != out.join("experiment.yaml") { fs::copy(config,out.join("experiment.yaml"))?; }
        let root=config.parent().unwrap_or(Path::new("."));
        for source in [&fixture.exp.reference.fasta,&fixture.exp.reference.gtf,&fixture.exp.reference.vcf] {
            fs::copy(root.join(source),out.join(source))?;
        }
    }
    verify(&fq,&reads)?;
    eprintln!("[nelrune] fixture timings: load={:?}, construct={:?}, write+verify={:?}", loaded.duration_since(started), synthesized.duration_since(loaded), synthesized.elapsed());
    eprintln!("[nelrune] {}: 43 reads (40 fragments + 3 PCR duplicates) / 4 UMIs / 2 haplotypes; reference {} ({}bp); validated FASTQ", if verify_only {"verified"} else {"generated"},fixture.chrom,fixture.genome.len());
    Ok(())
}
pub(super) fn execute(config:&Path,out:&Path,verify_only:bool)->Result<()> {
    let fixture = load(config)?;
    let grammar = grammar_for(&fixture)?;
    execute_with_grammar(config, out, verify_only, &grammar)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn round_trip_reference_backed_ont_fixture() -> Result<()> {
        let config=Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/nelrune_test_data_generation/experiment.yaml");
        let out=std::env::temp_dir().join(format!("lumrik_ont_fixture_{}",std::process::id()));
        // Both passes use the same grammar: expensive initialization happens only once.
        let fixture = load(&config)?;
        let grammar = grammar_for(&fixture)?;
        execute_with_grammar(&config,&out,false,&grammar)?;
        execute_with_grammar(&out.join("experiment.yaml"),&out,true,&grammar)?;
        // The CLI verify mode checks the FASTQ against reference-derived molecules.
        // Biological recovery belongs to the Nelrune -> mapper -> bam-quant integration test.
        fs::remove_dir_all(&out)?;
        Ok(())
    }
}
