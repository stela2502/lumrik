mod bam_quant;
mod materialize_norn_test_data;
mod prepare_fastqs;
mod quant;
mod vdj;
use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::process::{Command, Stdio};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{CommandFactory, FromArgMatches};
use rust_htslib::bam;

use bam_tide::FeatureTagCounts;
use bam_tide::illumina_normalizer::cli::{InsertRead, PrimerRead};
use bam_tide::illumina_normalizer::{IlluminaNormalizer, IlluminaNormalizerConfig};
use bam_tide::ont_normalizer::OntNormalizer;
use bam_tide::ont_normalizer::normalizer::OntNormalizerConfig;
use bam_tide::quantification::bam_collector::BamCollector;

use gtf_splice_index::SpliceIndex;
use sc_mapper::{MapperKind, MappingCall, StreamingMapper};
use scdata::QuantData;

use lumrik_status::{public_hostname, spawn_status_server};
use nelrune::cli::Cli;
use nelrune::progress::RunProgress;

fn main() -> Result<()> {
    match std::env::args().nth(1).as_deref() {
        Some("prepare-fastqs") => return prepare_fastqs::run(),
        Some("materialize-norn-test-data") => return materialize_norn_test_data::run(),
        Some("quant") => return quant::run(),
        Some("bam-quant") => return bam_quant::run(),
        Some("vdj") => return vdj::run(),
        _ => {}
    }

    let command = Cli::command();
    let matches = command.clone().get_matches();

    let args = Cli::from_arg_matches(&matches)?;

    args.validate()?;
    Cli::print_overrides(&command, &matches);

    run(args)
}

fn run(args: Cli) -> Result<()> {
    fs::create_dir_all(&args.outpath)
        .with_context(|| format!("creating output directory {}", args.outpath.display()))?;

    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    let mut progress = RunProgress::new();
    progress.start_timer("nelrune/startup");

    progress.open_log(args.outpath.join("nelrune.log"))?;

    let _health_server = if args.no_health_server {
        progress.stage("health server disabled");
        None
    } else {
        let health_addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), args.health_port);

        let server = spawn_status_server(progress.state_handle(), health_addr)?;
        let hostname = public_hostname(args.health_hostname.as_deref());
        let url = format!("http://{}:{}", hostname, server.addr().port());
        progress.set_public_url(url.clone());
        eprintln!("[nelrune] health server: {url}");
        Some(server)
    };

    progress.stage("loading chemistry");
    let mut primer = args
        .primer
        .detector()
        .map_err(anyhow::Error::msg)
        .context("configuring sc_primer chemistry")?;
    if args.barcode_crosscheck {
        primer = primer.with_barcode_crosscheck().map_err(anyhow::Error::msg)?;
        eprintln!("[barcode] crosscheck enabled (diagnostics only; no rescue)");
    }

    let cell_barcode_len = primer.cell_len();

    let mapper_bam = args
        .bam_collector
        .bam_out
        .as_ref()
        .map(|path| path.to_path_buf())
        .unwrap_or_else(|| args.outpath.join("nelrune.mapper.bam"));
    if let Some(parent) = mapper_bam.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("creating mapper BAM directory {}", parent.display()))?;
    }

    let (mut feature_counts, submitted) = if args.mapper.mapper == MapperKind::Minimap2 {
        let bam_input = args.bam.as_deref().ok_or_else(|| {
            anyhow::anyhow!("two-phase minimap2 currently requires ONT --bam input")
        })?;
        let finished = PathBuf::from(format!("{}.finished", mapper_bam.display()));
        if finished.is_file() && mapper_bam.is_file() {
            let metadata = fs::read_to_string(&finished)?;
            let expected = mapping_signature(&args, bam_input)?;
            let compatibility = mapping_signature_compatibility(&metadata, &expected);
            if !compatibility.core_matches {
                eprintln!("[nelrune] mapping completion marker mismatch: {}", finished.display());
                let stored: std::collections::BTreeMap<_, _> = metadata
                    .lines()
                    .filter_map(|line| line.split_once('='))
                    .collect();
                let current: std::collections::BTreeMap<_, _> = expected
                    .lines()
                    .filter_map(|line| line.split_once('='))
                    .collect();
                if metadata.lines().next() != expected.lines().next() {
                    eprintln!("  signature version: stored={:?}, expected={:?}",
                        metadata.lines().next(), expected.lines().next());
                }
                for key in stored.keys().chain(current.keys()).copied().collect::<std::collections::BTreeSet<_>>() {
                    if stored.get(key) != current.get(key) {
                        eprintln!("  {key}:");
                        eprintln!("    stored:   {}", stored.get(key).copied().unwrap_or("<missing>"));
                        eprintln!("    expected: {}", current.get(key).copied().unwrap_or("<missing>"));
                    }
                }
                bail!("mapping completion marker {} does not match current input/parameters; remove the marker to remap", finished.display());
            }
            if compatibility.mapping_limit_is_restrictive {
                eprintln!("[nelrune] WARNING: cached BAM was mapped with a more restrictive --max-reads limit.");
                eprintln!("[nelrune]          stored: {:?}; requested: {:?}", compatibility.stored_limit, compatibility.requested_limit);
                eprintln!("[nelrune]          Reusing existing BAM; results may omit reads not mapped previously.");
                eprintln!("[nelrune]          To regenerate complete mapping, remove {} and {} and rerun Nelrune.", finished.display(), mapper_bam.display());
            }
            check_bam(&mapper_bam)?;
            progress.stage("reusing completed minimap2 BAM");
            // Feature-tag counts from normalization are not persisted here.
            // Reusing a BAM must not silently discard them.
            if !args.additional_features.is_empty() {
                bail!("BAM reuse with --additional-features is not supported yet (normalization feature counts are not persisted)");
            }
            (FeatureTagCounts::default(), 0)
        } else {
            if finished.exists() { fs::remove_file(&finished)?; }
            // Normalization is independent of the mapper reference. Keep its
            // FASTQ and completion marker across remapping runs.
            let fastq = args.outpath.join("normalized.minimap2.fastq");
            let norm_marker = args.outpath.join("normalized.minimap2.fastq.finished");
            let norm_signature = normalization_signature(&args, bam_input)?;
            let reusable = fastq.is_file()
                && norm_marker.is_file()
                && fs::read_to_string(&norm_marker)
                    .map(|stored| stored == norm_signature)
                    .unwrap_or(false);
            let (count, features) = if reusable {
                if !args.additional_features.is_empty() {
                    bail!("cached ONT FASTQ cannot restore additional-feature counts; normalization must be rerun");
                }
                progress.stage("reusing normalized ONT FASTQ");
                eprintln!("[nelrune] reusing normalized FASTQ {}", fastq.display());
                (1, FeatureTagCounts::default())
            } else {
                if norm_marker.exists() {
                    fs::remove_file(&norm_marker)
                        .with_context(|| format!("removing stale normalization marker {}", norm_marker.display()))?;
                }
                progress.stage("normalizing ONT into persistent FASTQ");
                let result = normalize_ont_fastq(&args, &primer, bam_input, &mut progress)?;
                if result.0 == 0 { bail!("normalization emitted no mapper-ready reads"); }
                let tmp = norm_marker.with_extension("finished.partial");
                fs::write(&tmp, &norm_signature)
                    .with_context(|| format!("writing normalization marker {}", tmp.display()))?;
                fs::rename(&tmp, &norm_marker)
                    .with_context(|| format!("publishing normalization marker {}", norm_marker.display()))?;
                result
            };
            progress.stage("minimap2 | samtools: writing mapped BAM");
            map_minimap2_bam(&args, &mapper_bam)?;
            check_bam(&mapper_bam)?;
            let signature = mapping_signature(&args, bam_input)?;
            let marker_tmp = PathBuf::from(format!("{}.finished.partial", mapper_bam.display()));
            fs::write(&marker_tmp, signature)?;
            fs::rename(&marker_tmp, &finished)?;
            (features, count)
        }
    } else {
        progress.stage("starting mapper");
        let mut mapper = args.mapper.from_cli().context("starting streaming mapper")?;
        let mut sink = MapperBamSink::new(mapper_bam.clone());
        progress.stage("waiting for mapper / first batch");
        let (submitted, features) = if let Some(bam_input) = &args.bam {
            run_ont_input(&args, &primer, bam_input, &mut mapper, &mut sink, &mut progress)?
        } else {
            run_illumina_input(&args, &primer, &mut mapper, &mut sink, &mut progress)?
        };
        if submitted == 0 { bail!("normalization emitted no mapper-ready reads"); }
        progress.stage("finishing mapper");
        sink.ensure_open_blocking(&mut mapper)?;
        sink.drain_ready(&mut mapper)?;
        let remaining = mapper.finish().context("finishing streaming mapper")?;
        sink.write_calls(remaining)?;
        sink.finish()?;
        (features, submitted)
    };
    progress.clear_input_file();
    progress.stop_timer("nelrune/startup");
    if submitted > 0 { eprintln!("[nelrune] submitted {submitted} normalization batches"); }

    quantify_mapped_bam(
        &mapper_bam, &args.bam_collector, &args.outpath, args.min_umi_count,
        Some(primer.grammar().clone()), cell_barcode_len, feature_counts, &mut progress,
    )?;

    Ok(())
}

// Deliberately excludes mapper index, mapper threads, and quantification
// settings: none of these affect the normalized ONT read sequences.
fn normalization_signature(args: &Cli, input: &Path) -> Result<String> {
    let meta = fs::metadata(input)
        .with_context(|| format!("reading ONT input metadata {}", input.display()))?;
    let mtime = meta.modified()?.duration_since(std::time::UNIX_EPOCH)?.as_nanos();
    Ok(format!(
        "nelrune-ont-normalization-v1\ninput={}\nsize={}\nmtime={}\nprimer={:?}\nmax_reads={:?}\nmin_transcript_len={}\nadditional_features={:?}\nadditional_feature_min_hits={:?}\n",
        input.canonicalize()?.display(), meta.len(), mtime,
        args.primer, args.bam_collector.max_reads, args.min_transcript_len,
        args.additional_features, args.additional_feature_min_hits,
    ))
}

/// Quantify existing Nelrune alignments without repeating ONT normalization or mapping.
///
/// Encoded QNAME cell/UMI evidence is decoded by bam-tide; the same Bayesian
/// transcript resolution and matrix export are used by normal Nelrune runs.
pub(crate) fn quantify_mapped_bam(
    mapper_bam: &Path,
    config: &bam_tide::quantification::bam_collector::BamCollectorConfig,
    outpath: &Path,
    min_umi_count: Option<usize>,
    grammar: Option<sc_primer::Grammar>,
    cell_barcode_len: usize,
    mut feature_counts: FeatureTagCounts,
    progress: &mut RunProgress,
) -> Result<()> {
    progress.stage("quantifying mapper BAM");
    progress.start_timer("nelrune/quantification");

    let mut collector = BamCollector::from_cli(config.clone())
        .context("configuring BAM collector")?;
    if let Some(grammar) = grammar { collector = collector.with_grammar(grammar); }
    let mut result = collector
        .run_paths_with_progress(std::slice::from_ref(&mapper_bam.to_path_buf()), |data, report| {
            progress.update_quantification_live(data, report);
        })
        .context("collecting BAM quantification")?;
    progress.update_quantification_live(&result.data, &result.report);
    result.new_transcripts.write_reports(&outpath)?;
    let mut data = result.data;
    let mut report = result.report;
    progress.stop_timer("nelrune/quantification");

    let index = SpliceIndex::load(&config.index)
        .with_context(|| {
            format!(
                "reading splice index {} for export",
                config.index.display()
            )
        })?
        .with_match_mode(config.quant_mode.splice_match_mode());

    progress.stage("writing quantification");
    progress.start_timer("nelrune/writing");

    let (retained_cells, calling) = if let Some(min_umi_count) = min_umi_count {
        println!("Applying user-defined cell UMI cutoff: {}", min_umi_count);
        (data.cells_with_min_exonic_umis(min_umi_count), None)
    } else {
        println!("Running sc-beacon barcode-rank knee cell identification and QC...");

        let calling = bam_tide::results::beacon_cell_calling(&data).map_err(anyhow::Error::msg)?;

        std::fs::create_dir_all(&outpath)
            .with_context(|| format!("creating {}", outpath.display()))?;

        calling
            .write_tsv(outpath.join("cell_calling.tsv"))
            .map_err(anyhow::Error::msg)?;

        calling
            .write_qc(outpath.join("qc"))
            .map_err(anyhow::Error::msg)?;

        let retained = calling.retained.clone();
        (retained, Some(calling))
    };
    progress.stage(format!("cell cutoff complete: {} cells retained; preparing matrix indexes", retained_cells.len()));
    let features = index.feature_index();
    let mut indexes: HashMap<String, &dyn scdata::FeatureIndex> = HashMap::new();
    indexes.insert(QuantData::EXPRS.to_string(), &features);
    indexes.insert(QuantData::EXONIC.to_string(), &features);
    indexes.insert(QuantData::INTRONIC.to_string(), &features);
    if let Some(snp) = result.snp.as_ref() {
        indexes.insert(QuantData::SNP_REF.to_string(), &snp.index);
        indexes.insert(QuantData::SNP_ALT.to_string(), &snp.index);
    }
    progress.stage("writing raw and filtered quantification matrices");
    let cell_accounting = data
        .write_raw_and_filtered_for_cells(
            &outpath,
            &retained_cells,
            &indexes,
            Some(cell_barcode_len),
        )
        .map_err(anyhow::Error::msg)
        .context("writing quantification")?;

    // Additional-feature tables have their own feature-type-aware writer.
    // Preserve every observed feature-tag cell under raw/, then retain
    // only the canonical GEX cells under filtered/.
    progress.stage("writing raw additional-feature tables");
    feature_counts
        .finalize_and_write_all(cell_barcode_len, &outpath.join("raw"))
        .context("writing raw additional feature tables")?;
    progress.stage("writing filtered additional-feature tables");
    feature_counts
        .finalize_and_write(
            &retained_cells,
            cell_barcode_len,
            &outpath.join("filtered"),
        )
        .context("writing filtered additional feature tables")?;
    progress.stage("quantification output files written");
    progress.stop_timer("nelrune/writing");

    let bam_records_seen = report.get_issue_count("bam_records_seen");
    let quantified_bam_records = report.get_issue_count("cell_umi from BAM tags");
    let compatible_bam_records = report.get_issue_count("Compatible");
    let unmapped_bam_records = report.get_issue_count("unmapped");
    progress.set_quantification_summary(
        bam_records_seen,
        quantified_bam_records,
        compatible_bam_records,
        unmapped_bam_records,
        retained_cells.len(),
    );

    if let Some(calling) = &calling {
        progress.report_block(
            "Cell calling",
            &format!(
                "method: sc-beacon barcode-rank knee\ncandidate barcodes: {}\ninformative barcodes (>1 UMI): {}\nknee rank: {}\nUMI cutoff: {}\nknee score: {:.6}\nretained cells: {}\nnot called: {}\ncell diagnostics: {}\ncell QC plots: {}\n\n{}",
                calling.fit.candidate_barcodes,
                calling.fit.informative_barcodes,
                calling.fit.knee_rank,
                calling.fit.umi_cutoff,
                calling.fit.score,
                retained_cells.len(),
                cell_accounting.exonic_cells.saturating_sub(retained_cells.len()),
                outpath.join("cell_calling.tsv").display(),
                outpath.join("qc").display(),
                cell_accounting,
            ),
        );
    } else {
        let cutoff = min_umi_count
            .expect("fixed cell calling requires --min-umi-count");
        progress.report_block(
            "Cell calling",
            &format!(
                "method: user-defined UMI cutoff\nUMI cutoff: {}\nretained cells: {}\nnot called: {}\n\n{}",
                cutoff,
                retained_cells.len(),
                cell_accounting.exonic_cells.saturating_sub(retained_cells.len()),
                cell_accounting,
            ),
        );
    }

    // Preserve Nelrune's broad orchestration timings in the final report too.
    report.merge(progress.run_stats());

    progress.report_timings();
    progress.report_block("Quantification report", &report);

    fs::write(outpath.join("nelrune-report.txt"), report.to_string())
        .context("writing nelrune-report.txt")?;

    progress.stage(format!("quantified mapped BAM {}", mapper_bam.display()));

    progress.finish();
    progress.write_final_status(&outpath, "nelrune")?;

    eprintln!(
        "[nelrune] complete: {} reads in {:.1}s ({:.0} reads/s overall, {:.0} reads/s steady after first progress)",
        progress.reads_seen(),
        progress.elapsed().as_secs_f64(),
        progress.average_reads_per_second(),
        progress.processing_reads_per_second(),
    );

    Ok(())
}

fn mapping_signature(args: &Cli, input: &Path) -> Result<String> {
    let input_meta = fs::metadata(input)?;
    let index_meta = fs::metadata(&args.mapper.mapper_index)?;
    let modified = |m: &fs::Metadata| -> Result<u128> {
        Ok(m.modified()?.duration_since(std::time::UNIX_EPOCH)?.as_nanos())
    };
    Ok(format!(
        "nelrune-minimap2-v2\ninput={}\ninput_size={}\ninput_mtime={}\nindex={}\nindex_size={}\nindex_mtime={}\nprimer={:?}\nmapper_options={:?}\nmapper_threads={}\nmax_reads={:?}\nmin_transcript_len={}\nbarcode_crosscheck={}\n",
        input.canonicalize()?.display(), input_meta.len(), modified(&input_meta)?,
        args.mapper.mapper_index.canonicalize()?.display(), index_meta.len(), modified(&index_meta)?,
        args.primer, args.mapper.launch()?.options, args.mapper.mapper_threads,
        args.bam_collector.max_reads, args.min_transcript_len, args.barcode_crosscheck,
    ))
}

// Separate mapping-affecting parameters from diagnostics and read budgets.
// Older v1 completion markers remain readable without forcing remapping.
struct MappingCompatibility {
    core_matches: bool,
    stored_limit: Option<usize>,
    requested_limit: Option<usize>,
    mapping_limit_is_restrictive: bool,
}

fn mapping_signature_compatibility(stored: &str, expected: &str) -> MappingCompatibility {
    fn split(signature: &str) -> Option<(String, Option<usize>)> {
        let mut lines = signature.lines();
        let version = lines.next()?;
        let mut core = String::from(version);
        let mut limit = None;
        let mut found_limit = false;
        for line in lines {
            if let Some(value) = line.strip_prefix("max_reads=") {
                if found_limit { return None; }
                found_limit = true;
                limit = if value == "None" {
                    None
                } else {
                    Some(value.strip_prefix("Some(")?.strip_suffix(')')?.parse().ok()?)
                };
            } else if !line.starts_with("barcode_crosscheck=") {
                core.push('\n');
                core.push_str(line);
            }
        }
        if !found_limit { return None; }
        Some((core, limit))
    }
    let (stored_core, stored_limit) = match split(stored) {
        Some(v) => v,
        None => return MappingCompatibility { core_matches: false, stored_limit: None, requested_limit: None, mapping_limit_is_restrictive: false },
    };
    let (requested_core, requested_limit) = match split(expected) {
        Some(v) => v,
        None => return MappingCompatibility { core_matches: false, stored_limit, requested_limit: None, mapping_limit_is_restrictive: false },
    };
    let restrictive = match (stored_limit, requested_limit) {
        (None, _) => false,
        (Some(_), None) => true,
        (Some(available), Some(requested)) => available < requested,
    };
    MappingCompatibility {
        core_matches: stored_core == requested_core,
        stored_limit,
        requested_limit,
        mapping_limit_is_restrictive: restrictive,
    }
}

#[cfg(test)]
mod mapping_marker_tests {
    use super::mapping_signature_compatibility;

    fn signature(limit: &str, crosscheck: bool) -> String {
        format!("nelrune-minimap2-v2\ninput=example.bam\nmax_reads={limit}\nindex=example.mmi\nbarcode_crosscheck={crosscheck}\n")
    }

    #[test]
    fn diagnostics_do_not_invalidate_mapping() {
        let result = mapping_signature_compatibility(&signature("None", true), &signature("Some(2000000)", false));
        assert!(result.core_matches);
        assert!(!result.mapping_limit_is_restrictive);
    }

    #[test]
    fn restrictive_cached_mapping_warns_without_invalidating() {
        for requested in ["None", "Some(2000000)"] {
            let result = mapping_signature_compatibility(&signature("Some(1000000)", true), &signature(requested, false));
            assert!(result.core_matches);
            assert!(result.mapping_limit_is_restrictive);
        }
    }

    #[test]
    fn changed_index_still_invalidates_mapping() {
        let result = mapping_signature_compatibility(
            &signature("None", false),
            &signature("None", true).replace("index=example.mmi", "index=other.mmi"),
        );
        assert!(!result.core_matches);
    }
}

fn check_bam(path: &Path) -> Result<()> {
    let status = Command::new("samtools").arg("quickcheck").arg("-v").arg(path)
        .status().context("running samtools quickcheck")?;
    if !status.success() { bail!("samtools quickcheck failed for {}", path.display()); }
    Ok(())
}

fn normalize_ont_fastq(
    args: &Cli,
    primer: &sc_primer::PrimerDetector,
    bam_input: &Path,
    progress: &mut RunProgress,
) -> Result<(usize, FeatureTagCounts)> {
    let path = args.outpath.join("normalized.minimap2.fastq");
    let file = fs::File::create(&path)?;
    let mut output = BufWriter::new(file);
    let config = OntNormalizerConfig {
        bam: bam_input.to_path_buf(), dump_failed_reads: args.dump_failed_reads,
        out: args.outpath.join("normalized.fastq"),
        read_tags: args.outpath.join("read_tags.tsv"),
        min_transcript_len: args.min_transcript_len,
        ont_diagnostics: args.ont_diagnostics,
        max_reads: args.bam_collector.max_reads,
        primer: primer.clone().with_ont_fixed_search().map_err(anyhow::Error::msg)?, additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        threads: args.threads, gzip_level: 1, gzip: false,
    };
    let mut normalizer = OntNormalizer::new(config)?;
    let mut batches = 0usize;
    progress.input_file(bam_input.to_string_lossy());
    normalizer.nelrune_run(
        |reads| {
            for (_r1, r2) in reads { write!(output, "{r2}")?; }
            batches += 1;
            Ok(true)
        },
        |stats| progress.update_from_run_stats(stats),
    )?;
    output.flush()?;
    drop(output);
    progress.report_block(&format!("ONT normalization: {}", bam_input.display()), normalizer.stats());
    Ok((batches, normalizer.take_feature_tag_counts()))
}

fn map_minimap2_bam(args: &Cli, output: &Path) -> Result<()> {
    let fastq = args.outpath.join("normalized.minimap2.fastq");
    let partial = PathBuf::from(format!("{}.partial", output.display()));
    if partial.exists() { fs::remove_file(&partial)?; }
    let launch = args.mapper.launch()?;
    let options = launch.options.clone();
    // Match the existing sc-mapper minimap2 defaults without piping FASTQ
    // back and forth through Nelrune.
    let mut cleaned = Vec::new();
    let mut iter = options.into_iter();
    while let Some(opt) = iter.next() {
        if opt == "-t" { let _ = iter.next(); continue; }
        if opt == "-a" || opt == "-c" { continue; }
        cleaned.push(opt);
    }
    let mut options = cleaned;
    options.extend(["-t".into(), launch.threads.to_string(), "-a".into()]);
    if !options.iter().any(|option| option == "--eqx") {
        options.push("--eqx".into());
    }
    options.push(launch.index.to_string_lossy().into_owned());
    options.push(fastq.to_string_lossy().into_owned());
    let mut minimap = Command::new(&launch.mapper_bin)
        .args(&options).stdout(Stdio::piped()).spawn()
        .with_context(|| format!("starting {}", launch.mapper_bin.display()))?;
    let sam = minimap.stdout.take().context("minimap2 stdout unavailable")?;
    let mut samtools = match Command::new("samtools")
        .args(["view", "-b", "-o"]).arg(&partial).arg("-")
        .stdin(Stdio::from(sam)).spawn() {
        Ok(child) => child,
        Err(err) => { let _ = minimap.kill(); let _ = minimap.wait(); return Err(err).context("starting samtools view"); }
    };
    let map_status = minimap.wait().context("waiting for minimap2")?;
    let sam_status = samtools.wait().context("waiting for samtools")?;
    if !map_status.success() || !sam_status.success() {
        bail!("minimap2/samtools failed: minimap2={map_status}, samtools={sam_status}; partial BAM retained at {}", partial.display());
    }
    check_bam(&partial)?;
    fs::rename(&partial, output).context("publishing completed mapping BAM")?;
    Ok(())
}

fn run_illumina_input(
    args: &Cli,
    primer: &sc_primer::PrimerDetector,
    mapper: &mut StreamingMapper,
    sink: &mut MapperBamSink,
    progress: &mut RunProgress,
) -> Result<(usize, FeatureTagCounts)> {
    let mut submitted = 0usize;

    let config = IlluminaNormalizerConfig {
        out: args.outpath.join("normalized.fastq"),
        read_tags: args.outpath.join("read_tags.tsv"),
        primer_read: PrimerRead::R1,
        insert_read: InsertRead::R2,
        primer: primer.clone(),
        additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        min_insert_len: args.min_insert_len,
        max_reads: args.bam_collector.max_reads,
        threads: args.threads,
        gzip_level: 1,
        gzip: false,
    };

    let mut normalizer = IlluminaNormalizer::new(config)?;
    let mut first_progress_seen = false;

    for (r1_path, r2_path) in args.r1.iter().zip(args.r2.iter()) {
        /*
         * Local to this FASTQ pair.
         *
         * --max-reads therefore applies independently to every
         * R1/R2 input pair.
         */

        if first_progress_seen {
            progress.stage(format!("normalizing {}", display_name(r1_path)));
        }

        progress.input_file(r1_path.to_string_lossy());

        normalizer.nelrune_run(
            &r1_path.to_path_buf(),
            &r2_path.to_path_buf(),
            |reads| {
                for (r1, r2) in reads {
                    // Drain mapping results throughout FASTQ submission, not
                    // only after the entire normalization batch. Otherwise
                    // the bounded mapper-output channel can fill, blocking
                    // the stdout reader while the mapper blocks on stdout
                    // and this thread blocks writing FASTQ to its stdin/FIFO.
                    sink.drain_ready(mapper)?;
                    mapper.submit(r2, r1.as_ref())?;
                }

                submitted += 1;

                sink.drain_ready(mapper)?;

                Ok(true)
            },
            |stats| {
                if !first_progress_seen {
                    first_progress_seen = true;
                    progress.stage(format!("normalizing {}", display_name(r1_path)));
                }
                progress.update_from_run_stats(stats);
            },
        )?;

        progress.report_block(
            &format!("Illumina normalization: {}", r1_path.display()),
            normalizer.stats(),
        );
    }

    let feature_counts = normalizer.take_feature_tag_counts();

    Ok((submitted, feature_counts))
}

fn run_ont_input(
    args: &Cli,
    primer: &sc_primer::PrimerDetector,
    bam_input: &Path,
    mapper: &mut StreamingMapper,
    sink: &mut MapperBamSink,
    progress: &mut RunProgress,
) -> Result<(usize, FeatureTagCounts)> {
    progress.stage(format!("preparing {}", display_name(bam_input)));

    progress.input_file(bam_input.to_string_lossy());

    let config = OntNormalizerConfig {
        bam: bam_input.to_path_buf(),
        dump_failed_reads: args.dump_failed_reads,
        ont_diagnostics: args.ont_diagnostics,
        out: args.outpath.join("normalized.fastq"),
        read_tags: args.outpath.join("read_tags.tsv"),
        min_transcript_len: args.min_transcript_len,
        max_reads: args.bam_collector.max_reads,
        primer: primer.clone().with_ont_fixed_search().map_err(anyhow::Error::msg)?,
        additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        threads: args.threads,
        gzip_level: 1,
        gzip: false,
    };

    let mut normalizer = OntNormalizer::new(config)?;

    let mut submitted = 0usize;
    let mut first_progress_seen = false;

    normalizer.nelrune_run(
        |reads| {
            for (r1, r2) in reads {
                // Prevent a full mapper-result channel from stalling the
                // mapper's SAM/BAM stdout while FASTQ input is being sent.
                sink.drain_ready(mapper)?;
                mapper.submit(r2, r1.as_ref())?;
            }

            submitted += 1;

            sink.drain_ready(mapper)?;

            Ok(true)
        },
        |stats| {
            if !first_progress_seen {
                first_progress_seen = true;
                progress.stage(format!("normalizing {}", display_name(bam_input)));
            }
            progress.update_from_run_stats(stats);
        },
    )?;

    progress.report_block(
        &format!("ONT normalization: {}", bam_input.display()),
        normalizer.stats(),
    );

    let feature_counts = normalizer.take_feature_tag_counts();

    Ok((submitted, feature_counts))
}

fn display_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_else(|| path.to_str().unwrap_or("-"))
        .to_string()
}

struct MapperBamSink {
    path: PathBuf,
    writer: Option<bam::Writer>,
}

impl MapperBamSink {
    fn new(path: PathBuf) -> Self {
        Self { path, writer: None }
    }

    fn drain_ready(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        if self.writer.is_none() && mapper.header_loaded() {
            self.open_from_mapper(mapper)?;
        }

        if self.writer.is_none() {
            return Ok(());
        }

        while let Some(call) = mapper.try_next()? {
            self.write_call(call)?;
        }

        Ok(())
    }

    fn ensure_open_blocking(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        if self.writer.is_none() {
            self.open_from_mapper(mapper)?;
        }
        Ok(())
    }

    fn open_from_mapper(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        let header = mapper
            .header()
            .context("waiting for mapper BAM header")?
            .clone();

        let writer = bam::Writer::from_path(&self.path, &header, bam::Format::Bam)
            .with_context(|| format!("creating mapper BAM {}", self.path.display()))?;

        self.writer = Some(writer);
        Ok(())
    }

    fn write_call(&mut self, call: MappingCall) -> Result<()> {
        let writer = self
            .writer
            .as_mut()
            .context("mapper BAM writer not initialized")?;

        for record in call.records.records {
            let record = record.into_inner();
            writer.write(&record).context("writing mapper BAM record")?;
        }

        Ok(())
    }

    fn write_calls(&mut self, calls: Vec<MappingCall>) -> Result<()> {
        for call in calls {
            self.write_call(call)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.writer.take();
        Ok(())
    }
}
