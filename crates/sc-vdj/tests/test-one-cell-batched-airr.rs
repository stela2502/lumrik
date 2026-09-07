use sc_vdj::output::{write_mapping_info, ReportWriter};
use sc_vdj::{NelruneIdentityResolver, VdjIndex, VdjRunner, VdjRunnerConfig};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const CELL: &str = "CTTGGTCTTTTGGTTAATTCTTACCAT";

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("vdj-integration-cell")
}

#[test]
fn one_cell_real_fixture_detects_expected_recombinations() {
    let fixture = fixture_dir();
    let bam = fixture.join("cell.bam");
    let index_path = fixture.join("reference.vdjidx");

    assert!(bam.is_file(), "missing fixture BAM: {}", bam.display());

    assert!(
        index_path.is_file(),
        "missing VDJ index: {}",
        index_path.display()
    );

    let index = VdjIndex::load(&index_path)
        .unwrap_or_else(|e| panic!("loading {}: {e:#}", index_path.display()));

    let mut runner = VdjRunner::new(index, VdjRunnerConfig::default());

    let receptor_records = runner
        .read_bam(&bam, &NelruneIdentityResolver)
        .unwrap_or_else(|e| panic!("reading {}: {e:#}", bam.display()));

    assert!(
        receptor_records > 0,
        "real fixture produced no receptor-overlapping BAM evidence"
    );

    let calls = runner.identify();

    // Always write the current caller output before making biological assertions.
    // This makes a failed integration test useful for diagnosis instead of leaving
    // us staring at a stale output directory from an earlier implementation.
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .join("target")
        .join("sc-vdj-test-output")
        .join("one-cell-batched-airr");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir)
            .unwrap_or_else(|e| panic!("removing stale {}: {e}", out_dir.display()));
    }

    let mut writer = ReportWriter::create(&out_dir, true)
        .unwrap_or_else(|e| panic!("creating output in {}: {e:#}", out_dir.display()));
    let mut by_chain = HashMap::<String, usize>::new();
    let mut recombination_count = 0usize;
    for (cell_id, recombinations) in &calls {
        let cell_name = runner
            .cell_names
            .get(cell_id)
            .map(String::as_str)
            .unwrap_or("unknown");
        writer
            .write_cell(cell_name, recombinations, &runner.index)
            .unwrap_or_else(|e| panic!("writing cell {cell_name}: {e:#}"));
        for recombination in recombinations {
            *by_chain.entry(recombination.chain.to_string()).or_default() += 1;
            recombination_count += 1;
        }
    }
    writer
        .finish()
        .unwrap_or_else(|e| panic!("finishing output in {}: {e:#}", out_dir.display()));
    write_mapping_info(
        out_dir.join("vdj-mapping-info.txt"),
        receptor_records,
        calls.len(),
        recombination_count,
        &by_chain,
    )
    .unwrap_or_else(|e| panic!("writing mapping info in {}: {e:#}", out_dir.display()));

    eprintln!("sc-vdj real-fixture output: {}", out_dir.display());

    assert_eq!(calls.len(), 1, "fixture must resolve to exactly one cell");

    let (cell_id, recombinations) = &calls[0];

    assert_eq!(
        runner.cell_names.get(cell_id).map(String::as_str),
        Some(CELL),
        "wrong cell reconstructed from real fixture BAM"
    );

    assert_eq!(
        recombinations.len(),
        4,
        "expected exactly IGH, IGK and IGL recombinations"
    );

    let expected = HashMap::from([
        ("IGH", ("Ighv1-64", Some("Ighd1-1"), "Ighj3", Some("Igha"))),
        ("IGK", ("Igkv6-15", None, "Igkj2", Some("Igkc"))),
        ("IGL", ("Iglv3", None, "Iglj2", None)),
    ]);

    for recombination in recombinations {
        let chain = recombination.chain.to_string();

        let Some((expected_v, expected_d, expected_j, expected_c)) = expected.get(chain.as_str())
        else {
            panic!("unexpected chain in one-cell fixture: {chain}");
        };

        let v = &runner
            .index
            .segment(recombination.v)
            .unwrap_or_else(|| panic!("missing V segment id {:?}", recombination.v))
            .name;

        let d = recombination
            .d
            .and_then(|id| runner.index.segment(id))
            .map(|segment| segment.name.as_str());

        let j = &runner
            .index
            .segment(recombination.j)
            .unwrap_or_else(|| panic!("missing J segment id {:?}", recombination.j))
            .name;

        let c = recombination
            .constant
            .as_ref()
            .and_then(|constant| runner.index.segment(constant.segment))
            .map(|segment| segment.name.as_str());

        assert_eq!(v, expected_v, "wrong V call for {chain}");

        assert_eq!(d, *expected_d, "wrong D call for {chain}");

        assert_eq!(j, expected_j, "wrong J call for {chain}");

        assert_eq!(c, *expected_c, "wrong C call for {chain}");

        assert!(
            !recombination.observed_rearrangement.is_empty(),
            "empty observed rearrangement for {chain}"
        );

        assert!(
            !recombination.naive_recombination.is_empty(),
            "empty naive recombination for {chain}"
        );

        assert!(
            !recombination.observed_receptor_sequence.is_empty(),
            "empty observed receptor sequence for {chain}"
        );

        let rec_id = recombination.stable_id.to_string();

        let expected_prefix = if recombination.chain.has_d() {
            "HC:"
        } else {
            "LC:"
        };

        assert!(
            rec_id.starts_with(expected_prefix),
            "wrong recombination ID for {chain}: {rec_id}"
        );
    }
}
