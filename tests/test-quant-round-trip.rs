//test-quant-round-trip.rs
use anyhow::{Context, Result};
use assert_cmd::Command;
use tempfile::tempdir;

use bam_tide::index::GeneFeatureIndex;
use bam_tide::results::QuantData;

use gtf_splice_index::{IdNameKeys, SpliceIndex};
use snp_index::{Genome, SnpIndex, VcfReadOptions};

// Adjust these imports to wherever your generator lives.
use bam_tide::encoder::{SamEncoder, TestDataCli, cli::TruthFeatureMode};

#[test]
#[ignore = "disabled in CI for now"]
fn quant_round_trip_artificial_genome_gene_mode() -> Result<()> {
    let tmp = tempdir()?;

    let keep = std::env::var("KEEP_TMP").is_ok();

    let tmp_path = if keep {
        let p = tmp.keep(); // persist
        println!("Keeping tmp dir at {}", p.display());
        p
    } else {
        tmp.path().to_path_buf()
    };

    let fasta = "tests/data/art_genome.fa";
    let splice_index_path = "tests/data/art_genome_info.gtf.dat";
    let vcf = "tests/data/art_genome_snps.vcf";

    let sam_path = tmp_path.as_path().join("round_trip.sam");
    let truth_out = tmp_path.as_path().join("truth");
    let quant_out = tmp_path.as_path().join("quant");

    let splice_index = SpliceIndex::from_path(
        "tests/data/art_genome_info.gtf",
        100_000,
        IdNameKeys::default(),
    )?;
    splice_index
        .save(&splice_index_path)
        .map_err(|e| anyhow::anyhow!("failed to save index: {e}"))?;
    // ------------------------------------------------------------
    // 1) Generate SAM + truth from fixed artificial genome fixtures
    // ------------------------------------------------------------
    let cli = TestDataCli {
        gtf: splice_index_path.into(),
        fasta: fasta.into(),
        outfile: sam_path.clone(),
        vcf: Some(vcf.into()),

        // ---- defaults (copied from clap) ----
        seq_len: 30,
        reads_per_cell: 70,
        antisense_fraction: 10.0 / 70.0,
        unspliced_fraction: 10.0 / 70.0,
        snp_fraction: 0.20,
        body_error_rate: 0.005,
        end_error_rate: 0.15,
        bad_end_bases: 4,
        seed: 1,

        cell_barcodes: vec![
            "AAACCGCTCTATCCTA-1".to_string(),
            "AAACGAATCACCAGAC-1".to_string(),
            "AAACGAATCCGCGAAT-1".to_string(),
        ],

        genes_per_cell: 100,

        truth_path: Some(truth_out.clone()),
        min_cell_counts: 90,
        truth_feature_mode: TruthFeatureMode::Gene,
    };
    let mut encoder = SamEncoder::new(cli).map_err(anyhow::Error::msg)?;
    let lines = encoder.generate().map_err(anyhow::Error::msg)?;

    std::fs::write(&sam_path, lines.join("\n") + "\n")
        .with_context(|| format!("failed to write generated SAM to {}", sam_path.display()))?;

    assert!(
        sam_path.exists(),
        "SAM was not generated at {}",
        sam_path.display()
    );

    assert!(
        sam_path.exists(),
        "SAM was not generated at {}",
        sam_path.display()
    );

    // ------------------------------------------------------------
    // 2) Run bam-quant on the generated SAM
    // ------------------------------------------------------------
    Command::cargo_bin("bam-quant")?
        .args([
            "--bam",
            sam_path.to_str().context("SAM path is not UTF-8")?,
            "--index",
            splice_index_path,
            "--genome",
            fasta,
            "--vcf",
            vcf,
            "--outpath",
            quant_out.to_str().context("quant path is not UTF-8")?,
            "--min-cell-counts",
            "1",
            "--quant-mode",
            "gene",
            "--split-intronic",
        ])
        .assert()
        .success();

    // ------------------------------------------------------------
    // 3) Recreate indices exactly like SamEncoder does
    // ------------------------------------------------------------
    let genome = Genome::from_fasta(fasta)
        .map_err(|e| anyhow::anyhow!("failed to load FASTA {fasta}: {e}"))?;

    let splice_index = SpliceIndex::load(splice_index_path)
        .map_err(|e| anyhow::anyhow!("failed to load splice index {splice_index_path}: {e}"))?;

    let gene_index = GeneFeatureIndex::new(&splice_index);

    let snp_index = SnpIndex::from_vcf_path(
        vcf,
        genome.chr_names.clone(),
        genome.chr_lengths(),
        10_000,
        &VcfReadOptions::default(),
    )
    .map_err(|e| anyhow::anyhow!("failed to load VCF {vcf}: {e}"))?;

    // ------------------------------------------------------------
    // 4) Read truth + quantification output and compare
    // ------------------------------------------------------------
    let truth = QuantData::from_path(&truth_out, &gene_index, &snp_index)
        .map_err(|e| anyhow::anyhow!("failed to read truth QuantData: {e}"))?;

    let quant = QuantData::from_path(&quant_out, &gene_index, &snp_index)
        .map_err(|e| anyhow::anyhow!("failed to read observed QuantData: {e}"))?;

    truth
        .compare(&quant)
        .map_err(|e| anyhow::anyhow!("round-trip quantification mismatch:\n{e}"))?;

    Ok(())
}
