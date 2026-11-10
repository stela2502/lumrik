//! End-to-end regression: YAML -> reference-derived reads -> ONT BAM -> Nelrune -> mapper -> quant.
//! Requires STAR, minimap2 and samtools on PATH. Failed outputs are deliberately retained.
use std::{collections::{HashMap, HashSet}, fs::{self, File}, io::{BufRead, BufReader, Read}, path::{Path, PathBuf}, process::Command};
use flate2::read::GzDecoder;
use gtf_splice_index::AnnotationBuilder;
use rust_htslib::bam::{self, Read as BamRead};

fn root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..") }
fn check_tool(name: &str, arg: &str) {
    assert!(Command::new(name).arg(arg).output().is_ok(),
        "{name} is required for the ONT integration test; install it and retry");
}
fn run(label: &str, command: &mut Command, logs: &Path) {
    let output = command.output().unwrap_or_else(|e| panic!("starting {label}: {e}"));
    fs::write(logs.join(format!("{label}.stdout.txt")), &output.stdout).unwrap();
    fs::write(logs.join(format!("{label}.stderr.txt")), &output.stderr).unwrap();
    assert!(output.status.success(), "{label} failed ({}).\nstdout:\n{}\nstderr:\n{}\nLogs retained in {}",
        output.status, String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr), logs.display());
}
fn write_ont_bam(fastq: &Path, output: &Path) {
    let input = BufReader::new(File::open(fastq).unwrap());
    let mut lines = input.lines();
    let mut header = bam::Header::new();
    header.push_record(bam::header::HeaderRecord::new(b"HD")
        .push_tag(b"VN", "1.6")
        .push_tag(b"SO", "unsorted"));
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam).unwrap();
    let mut count = 0;
    while let Some(name) = lines.next() {
        let name = name.unwrap();
        let seq = lines.next().expect("FASTQ sequence").unwrap();
        assert_eq!(lines.next().unwrap().unwrap(), "+");
        let quality = lines.next().expect("FASTQ quality").unwrap();
        assert_eq!(seq.len(), quality.len());
        let mut record = bam::Record::new();
        record.set(name.strip_prefix('@').expect("FASTQ @ header").as_bytes(),
            None, seq.as_bytes(), &vec![40; seq.len()]);
        record.set_flags(4); // Unaligned raw ONT read; no CB/UB tags pre-assigned.
        writer.write(&record).unwrap();
        count += 1;
    }
    assert_eq!(count, 43, "expected 40 original observations + 3 PCR copies");
}
fn uncompress(source: &Path, dest: &Path) {
    let mut decoder = GzDecoder::new(File::open(source).unwrap());
    let mut output = File::create(dest).unwrap();
    std::io::copy(&mut decoder, &mut output).unwrap();
}
fn summary_number(path: &Path, field: &str) -> u64 {
    let yaml: serde_yaml::Value = serde_yaml::from_reader(File::open(path).unwrap()).unwrap();
    yaml.get(field).and_then(serde_yaml::Value::as_u64)
        .unwrap_or_else(|| panic!("missing numeric field {field} in {}", path.display()))
}
fn count_distinct_truth_molecules(path: &Path) -> usize {
    let text = fs::read_to_string(path).unwrap();
    let mut molecules = HashSet::new();
    for line in text.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 6, "invalid truth TSV: {line}");
        molecules.insert((fields[1].to_owned(), fields[2].to_owned()));
    }
    assert_eq!(text.lines().count() - 1, 43);
    molecules.len()
}
fn inspect_alignment(path: &Path) -> (usize, usize, usize) {
    let mut reader = bam::Reader::from_path(path).unwrap();
    let mut mapped = 0;
    let mut spliced = 0;
    let mut exon_skip = 0;
    for item in reader.records() {
        let rec = item.unwrap();
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() { continue; }
        mapped += 1;
        let mut reference_pos = rec.pos();
        let mut has_splice = false;
        for op in rec.cigar().iter() {
            use bam::record::Cigar;
            match op {
                Cigar::RefSkip(len) => {
                    has_splice = true;
                    // Reference is 0-based: exon 1 ends at 380; exon 3 begins at 1300.
                    if reference_pos == 380 && reference_pos + i64::from(*len) == 1300 {
                        exon_skip += 1;
                    }
                    reference_pos += i64::from(*len);
                }
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::Del(len) => {
                    reference_pos += i64::from(*len);
                }
                _ => {}
            }
        }
        if has_splice { spliced += 1; }
    }
    (mapped, spliced, exon_skip)
}


const MUTANT_TRANSCRIPT: &str = "SYNTH_SPLICE4-001_skip_chrSyntheticMouse1:751-930";

fn read_gzip(path: &Path) -> String {
    let mut text = String::new();
    GzDecoder::new(File::open(path).unwrap())
        .read_to_string(&mut text).unwrap();
    text
}

fn assert_novel_transcript(mapper: &str, output: &Path) -> usize {
    let gtf = fs::read_to_string(output.join("new_transcripts.gtf"))
        .unwrap_or_else(|e| panic!("{mapper}: missing novel-transcript GTF: {e}"));
    let exons: Vec<(u64, u64)> = gtf.lines().filter_map(|line| {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), 9, "{mapper}: malformed novel GTF: {line}");
        if !fields[8].contains(&format!("transcript_id \"{MUTANT_TRANSCRIPT}\"")) {
            return None;
        }
        assert_eq!(fields[0], "chrSyntheticMouse1");
        assert_eq!(fields[2], "exon");
        Some((fields[3].parse().unwrap(), fields[4].parse().unwrap()))
    }).collect();
    assert_eq!(exons, [(201, 380), (1301, 1480), (1851, 2030)],
        "{mapper}: incorrect reconstructed exon-skipping structure");

    // The report uses packed numeric cell/UMI identifiers, not DNA strings.
    // Verify two distinct UMIs, each with ALT at SNP 1 and REF at SNP 2.
    let snp = fs::read_to_string(output.join("new_transcript_snp_evidence.tsv"))
        .unwrap_or_else(|e| panic!("{mapper}: missing transcript-associated SNP evidence: {e}"));
    assert_eq!(snp.lines().next(),
        Some("cell_id\ttranscript_id\tumi\tsnp_id\tref_reads\talt_reads"));
    let mut molecules: HashMap<(u64, u64), HashMap<u64, (u64, u64)>> = HashMap::new();
    for line in snp.lines().skip(1) {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(f.len(), 6, "{mapper}: malformed SNP evidence: {line}");
        if f[1] != MUTANT_TRANSCRIPT { continue; }
        let key = (f[0].parse().unwrap(), f[2].parse().unwrap());
        let snp_id: u64 = f[3].parse().expect("novel transcript must have SNP evidence");
        let counts = (f[4].parse().unwrap(), f[5].parse().unwrap());
        assert!(molecules.entry(key).or_default().insert(snp_id, counts).is_none(),
            "{mapper}: duplicate SNP evidence for one molecule: {line}");
    }
    assert_eq!(molecules.len(), 2, "{mapper}: expected two mutant transcript UMIs");
    assert_eq!(molecules.keys().map(|(cell, _)| cell).collect::<HashSet<_>>().len(), 1,
        "{mapper}: mutant molecules must belong to one cell");
    for (molecule, snps) in &molecules {
        assert_eq!(snps.get(&1), Some(&(0, 1)),
            "{mapper}: mutant molecule {molecule:?} missing exon-3 ALT reporter SNP");
        assert_eq!(snps.get(&2), Some(&(1, 0)),
            "{mapper}: mutant molecule {molecule:?} missing exon-4 REF reporter SNP");
    }

    let support = fs::read_to_string(output.join("new_transcript_supported_regions.tsv"))
        .unwrap_or_else(|e| panic!("{mapper}: missing supported-regions report: {e}"));
    let supported: HashSet<(u64, u64)> = support.lines().skip(1).filter_map(|line| {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(f.len(), 6, "{mapper}: malformed supported region: {line}");
        (f[1] == MUTANT_TRANSCRIPT).then(|| (f[0].parse().unwrap(), f[2].parse().unwrap()))
    }).collect();
    assert_eq!(supported, molecules.keys().copied().collect(),
        "{mapper}: transcript/SNP molecules differ from supported-region molecules");

    // MatrixMarket coordinate format is 1-based: find the actual transcript feature
    // and sum its UMI counts rather than assuming a feature ordering.
    let dir = output.join("new_transcripts_mtx");
    let features = read_gzip(&dir.join("features.tsv.gz"));
    let feature_indices: Vec<usize> = features.lines().enumerate()
        .filter(|(_, line)| line.split('\t').any(|f| f == MUTANT_TRANSCRIPT))
        .map(|(i, _)| i + 1).collect();
    assert_eq!(feature_indices.len(), 1,
        "{mapper}: expected exactly one mutant feature in novel-transcript matrix; features: {features}");
    let matrix = read_gzip(&dir.join("matrix.mtx.gz"));
    let mut rows = matrix.lines().filter(|line| !line.starts_with('%'));
    let dims: Vec<usize> = rows.next().expect("MatrixMarket dimensions")
        .split_whitespace().map(|x| x.parse().unwrap()).collect();
    assert_eq!(dims.len(), 3, "{mapper}: malformed matrix dimensions");
    let mut expression = 0.0f64;
    let mut entries = 0;
    for line in rows {
        let f: Vec<_> = line.split_whitespace().collect();
        assert_eq!(f.len(), 3, "{mapper}: malformed matrix entry: {line}");
        if f[0].parse::<usize>().unwrap() == feature_indices[0] {
            expression += f[2].parse::<f64>().unwrap();
            entries += 1;
        }
    }
    assert_eq!(expression, 2.0, "{mapper}: expected mutant transcript expression of two UMIs");
    assert_eq!(entries, 1, "{mapper}: mutant transcript should be expressed in one cell");
    molecules.len()
}

#[test]
fn ont_reference_backed_star_and_minimap2_recover_molecules() {
    check_tool("STAR", "--version");
    check_tool("minimap2", "--version");
    check_tool("samtools", "--version");

    let source = root().join("tests/data/nelrune_test_data_generation");
    let out = root().join("target/nelrune-ont-reference-integration");
    if out.exists() { fs::remove_dir_all(&out).unwrap(); }
    fs::create_dir_all(&out).unwrap();
    let binary = env!("CARGO_BIN_EXE_nelrune");
    run("generate", Command::new(binary).args(["materialize-norn-test-data", "--config"])
        .arg(source.join("experiment.yaml")).arg("--out").arg(&out), &out);
    assert_eq!(count_distinct_truth_molecules(&out.join("expected_reads.tsv")), 4,
        "fixture must contain four distinct cell/UMI pairs");

    let fasta = out.join("genome.fa");
    let gtf = out.join("genes.gtf");
    let vcf = out.join("variants.vcf.gz");
    assert!(vcf.is_file(), "missing synthetic SNP reference: {}", vcf.display());
    uncompress(&out.join("genome.fa.gz"), &fasta);
    uncompress(&out.join("genes.gtf.gz"), &gtf);
    let index = out.join("reference.splice.idx");
    AnnotationBuilder::new(1000).build_from_path(&gtf).unwrap().save(&index).unwrap();

    let mmi = out.join("genome.mmi");
    run("minimap2-index", Command::new("minimap2").arg("-d").arg(&mmi).arg(&fasta), &out);
    let star = out.join("star");
    fs::create_dir_all(&star).unwrap();
    run("star-index", Command::new("STAR")
        .args(["--runMode", "genomeGenerate", "--runThreadN", "2", "--genomeDir"])
        .arg(&star).arg("--genomeFastaFiles").arg(&fasta)
        .arg("--sjdbGTFfile").arg(&gtf)
        .args(["--genomeSAindexNbases", "3", "--genomeChrBinNbits", "12", "--sjdbOverhang", "99"]), &out);
    assert!(star.join("Genome").is_file());

    // Read the exact same grammar used by the fixture generator from its YAML.
    let experiment: serde_yaml::Value = serde_yaml::from_reader(
        File::open(out.join("experiment.yaml")).unwrap()).unwrap();
    let primer_structure = experiment["primer_structure"].as_str()
        .expect("fixture YAML must define primer_structure for ONT detection");

    let bam = out.join("ont_reads.bam");
    write_ont_bam(&out.join("ont_reads.fastq"), &bam);
    for (mapper, mapper_index) in [("minimap2", &mmi), ("star", &star)] {
        let mapped_out = out.join(mapper);
        fs::create_dir_all(&mapped_out).unwrap();
        run(mapper, Command::new(binary)
            .arg("--bam").arg(&bam)
            .args(["--chemistry", "tenx-three-prime-v3", "--primer-structure", primer_structure,
                "--mapper", mapper, "--mapper-index"])
            .arg(mapper_index)
            .arg("--genome").arg(&fasta)
            .arg("--vcf").arg(&vcf)
            .args(["--mapper-threads", "2", "--index"])
            .arg(&index)
            .args(["--quant-mode", "transcript", "--min-umi-count", "1", "--min-mapq", "0", "--min-transcript-len", "20",
                "--threads", "2", "--no-health-server", "--outpath"])
            .arg(&mapped_out), &out);
        let summary = mapped_out.join("nelrunelog.jaml");
        assert!(summary.is_file(), "{mapper}: missing run summary {}", summary.display());
        assert_eq!(summary_number(&summary, "reads_processed"), 43,
            "{mapper}: ONT input read count differs from fixture");
        assert_eq!(summary_number(&summary, "exonic_umis"), 2,
            "{mapper}: expected two canonical gene-associated molecules in transcript mode");
        let (mapped, spliced, exon_skip) = inspect_alignment(&mapped_out.join("nelrune.mapper.bam"));
        assert!(mapped > 0, "{mapper}: no primary mapped reads");
        assert!(spliced > 0, "{mapper}: no spliced genomic alignments; exon junctions lost");
        assert!(exon_skip > 0,
            "{mapper}: failed to reconstruct mutant exon 1 -> exon 3 junction (chrSyntheticMouse1:380 -> 1300)");
        let novel_umis = assert_novel_transcript(mapper, &mapped_out);
        assert_eq!(summary_number(&summary, "exonic_umis") as usize + novel_umis, 4,
            "{mapper}: canonical and novel transcript UMIs do not account for all four molecules");
        eprintln!("[fixture integration] {mapper}: {mapped} mapped primary records, {spliced} spliced, {exon_skip} exon-skip alignments, 2 canonical + {novel_umis} mutant transcript UMIs with SNP evidence and matrix expression");
    }
    // Keep generated indexes, run summaries and alignments for diagnosis and reruns.
    eprintln!("[fixture integration] results: {}", out.display());
}
