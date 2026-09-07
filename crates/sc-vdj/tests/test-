use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

const CELL: &str = "CTTGGTCTTTTGGTTAATTCTTACCAT";

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vdj-integration-cell")
}

fn test_output_dir() -> PathBuf {
    if let Some(target_dir) = std::env::var_os("CARGO_TARGET_DIR") {
        PathBuf::from(target_dir)
            .join("sc-vdj-test-output")
            .join("one-cell-batched-airr")
    } else {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/sc-vdj-test-output/one-cell-batched-airr")
    }
}

fn parsed_batch_flushes(stderr: &str) -> Option<usize> {
    let marker = "batches=";
    let start = stderr.rfind(marker)? + marker.len();
    let digits: String = stderr[start..]
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    (!digits.is_empty()).then(|| digits.parse().ok()).flatten()
}

fn tsv_rows<'a>(text: &'a str) -> (Vec<&'a str>, Vec<HashMap<&'a str, &'a str>>) {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let header: Vec<&str> = lines.next().expect("TSV header").split('\t').collect();
    let rows = lines
        .map(|line| {
            let values: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                values.len(),
                header.len(),
                "TSV row has {} fields but header has {}:\n{}",
                values.len(),
                header.len(),
                line
            );
            header.iter().copied().zip(values).collect()
        })
        .collect();
    (header, rows)
}

#[test]
fn one_cell_fixture_crosses_batch_boundaries_and_writes_expected_airr_calls() {
    let fixture = fixture_dir();
    let exonic = fixture.join("exonic");
    let bam = fixture.join("cell.bam");
    let index = fixture.join("reference.vdjidx");

    assert!(exonic.is_dir(), "missing fixture MEX: {}", exonic.display());
    assert!(bam.is_file(), "missing fixture BAM: {}", bam.display());
    assert!(index.is_file(), "missing VDJ index: {}", index.display());

    // Keep failed test output for inspection. A new iteration removes the
    // previous output first; successful runs clean up at the very end.
    let out = test_output_dir();
    if out.exists() {
        fs::remove_dir_all(&out).expect("remove previous test output");
    }
    fs::create_dir_all(&out).expect("create persistent test output");

    let output = Command::new(env!("CARGO_BIN_EXE_nelrune-vdj"))
        .arg("--exonic")
        .arg(&exonic)
        .arg("--bam")
        .arg(&bam)
        .arg("--index")
        .arg(&index)
        .arg("--out")
        .arg(&out)
        .arg("--routing-batch-size")
        .arg("100")
        .arg("--no-health-server")
        .output()
        .expect("run nelrune-vdj");

    fs::write(out.join("nelrune-vdj.stdout.txt"), &output.stdout)
        .expect("write nelrune-vdj stdout");
    fs::write(out.join("nelrune-vdj.stderr.txt"), &output.stderr)
        .expect("write nelrune-vdj stderr");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "nelrune-vdj failed with {}\nSTDOUT:\n{}\nSTDERR:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        stderr
    );

    let batches = parsed_batch_flushes(&stderr)
        .expect("nelrune-vdj stderr did not report routing batches");
    assert!(
        batches >= 2,
        "batch size 100 did not force repeated compaction; batches={batches}\n{stderr}"
    );

    let mapping = fs::read_to_string(out.join("vdj-mapping-info.txt"))
        .expect("read vdj-mapping-info.txt");
    assert!(mapping.contains("vdj.routing_batch_flushes"));
    assert!(mapping.contains("vdj.compacted_sequences"));
    assert!(mapping.contains("vdj.redundant_segments"));

    // First pin the rich Lumrik calls to the established one-cell result. This
    // catches regressions in the biological caller independently of AIRR
    // serialization.
    let rich = fs::read_to_string(out.join("vdj_calls.tsv"))
        .expect("read vdj_calls.tsv");
    let (_header, rich_rows) = tsv_rows(&rich);
    assert_eq!(rich_rows.len(), 3, "expected exactly IGH, IGK and IGL calls");

    let expected = HashMap::from([
        ("IGH", ("Vdj", "Ighv1-64", "Ighd1-1", "Ighj3", "Igha", "2")),
        ("IGK", ("Vj", "Igkv6-15", "", "Igkj2", "Igkc", "9")),
        ("IGL", ("Vj", "Iglv3", "", "Iglj2", "", "1")),
    ]);

    for row in &rich_rows {
        assert_eq!(row["cell"], CELL);
        let chain = row["chain"];
        let Some((stage, v, d, j, c, umis)) = expected.get(chain) else {
            panic!("unexpected chain in one-cell fixture: {chain}");
        };
        assert_eq!(row["stage"], *stage, "wrong stage for {chain}");
        assert_eq!(row["v"], *v, "wrong V call for {chain}");
        assert_eq!(row["d"], *d, "wrong D call for {chain}");
        assert_eq!(row["j"], *j, "wrong J call for {chain}");
        assert_eq!(row["c"], *c, "wrong C call for {chain}");
        assert_eq!(row["support_umis"], *umis, "wrong UMI support for {chain}");
    }

    // Now pin the AIRR export to exactly the same three rearrangements and
    // require the sequence-level annotations introduced by the AIRR writer.
    let airr = fs::read_to_string(out.join("airr_rearrangements.tsv"))
        .expect("read airr_rearrangements.tsv");
    let (airr_header, airr_rows) = tsv_rows(&airr);
    for required in [
        "sequence_id",
        "sequence",
        "sequence_aa",
        "productive",
        "vj_in_frame",
        "stop_codon",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
        "junction",
        "junction_aa",
        "cdr3",
        "cdr3_aa",
        "locus",
        "cell_id",
        "lumrik_supporting_umis",
        "lumrik_v_cys_anchor",
        "lumrik_j_anchor",
    ] {
        assert!(
            airr_header.contains(&required),
            "AIRR output lacks required/tested column {required}"
        );
    }
    assert_eq!(airr_rows.len(), 3, "AIRR must contain exactly three calls");

    for row in &airr_rows {
        assert_eq!(row["cell_id"], CELL);
        let chain = row["locus"];
        let Some((_stage, v, d, j, c, umis)) = expected.get(chain) else {
            panic!("unexpected AIRR locus in one-cell fixture: {chain}");
        };
        assert_eq!(row["v_call"], *v, "wrong AIRR V call for {chain}");
        assert_eq!(row["d_call"], *d, "wrong AIRR D call for {chain}");
        assert_eq!(row["j_call"], *j, "wrong AIRR J call for {chain}");
        assert_eq!(row["c_call"], *c, "wrong AIRR C call for {chain}");
        assert_eq!(
            row["lumrik_supporting_umis"], *umis,
            "wrong AIRR UMI support for {chain}"
        );

        assert!(!row["sequence_id"].is_empty(), "missing sequence_id for {chain}");
        assert!(!row["sequence"].is_empty(), "missing sequence for {chain}");
        assert!(!row["sequence_aa"].is_empty(), "missing sequence_aa for {chain}");
        assert!(!row["junction"].is_empty(), "missing junction for {chain}");
        assert!(!row["junction_aa"].is_empty(), "missing junction_aa for {chain}");
        assert!(!row["cdr3"].is_empty(), "missing cdr3 for {chain}");
        assert!(!row["cdr3_aa"].is_empty(), "missing cdr3_aa for {chain}");

        for field in [
            "productive",
            "vj_in_frame",
            "stop_codon",
            "lumrik_v_cys_anchor",
            "lumrik_j_anchor",
        ] {
            assert!(
                matches!(row[field], "T" | "F"),
                "{field} must be AIRR T/F for {chain}, got {:?}",
                row[field]
            );
        }

        let junction_aa = row["junction_aa"].as_bytes();
        let cdr3_aa = row["cdr3_aa"].as_bytes();
        assert_eq!(
            junction_aa.len(),
            cdr3_aa.len() + 2,
            "junction_aa must contain exactly the two conserved anchor residues around cdr3_aa for {chain}"
        );
        assert_eq!(junction_aa.first(), Some(&b'C'), "V anchor is not cysteine for {chain}");
        assert!(
            matches!(junction_aa.last(), Some(&b'W') | Some(&b'F')),
            "J anchor is not W/F for {chain}: {}",
            row["junction_aa"]
        );
        assert_eq!(
            &junction_aa[1..junction_aa.len() - 1],
            cdr3_aa,
            "CDR3 amino acid sequence is inconsistent with junction_aa for {chain}"
        );
    }

    fs::remove_dir_all(&out).expect("remove successful test output");
}
