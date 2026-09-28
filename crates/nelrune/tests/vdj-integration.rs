use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

const CELL: &str = "CTTGGTCTTTTGGTTAATTCTTACCAT";

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../sc-vdj/tests/data/vdj-integration-cell")
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

fn field<'a>(row: &'a HashMap<&str, &'a str>, key: &str, context: &str) -> &'a str {
    row.get(key).copied().unwrap_or_else(|| {
        let mut available: Vec<_> = row.keys().copied().collect();
        available.sort_unstable();
        panic!(
            "missing TSV column {key:?} while checking {context}\navailable columns: {}\nrow: {row:#?}",
            available.join(", ")
        );
    })
}

#[test]
fn one_cell_fixture_writes_expected_airr_calls() {
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

    let output = Command::new(env!("CARGO_BIN_EXE_nelrune"))
        .arg("vdj")
        .arg("--exonic")
        .arg(&exonic)
        .arg("--bam")
        .arg(&bam)
        .arg("--index")
        .arg(&index)
        .arg("--out")
        .arg(&out)
        .arg("--bd-cell-version")
        .arg("v2.384")
        .arg("--no-health-server")
        .output()
        .expect("run nelrune vdj");

    fs::write(out.join("nelrune-vdj.stdout.txt"), &output.stdout)
        .expect("write nelrune-vdj stdout");
    fs::write(out.join("nelrune-vdj.stderr.txt"), &output.stderr)
        .expect("write nelrune-vdj stderr");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "nelrune vdj failed with {}\nSTDOUT:\n{}\nSTDERR:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        stderr
    );

    // First pin the rich Lumrik calls to the established one-cell result. This
    // catches regressions in the biological caller independently of AIRR
    // serialization.
    let rich = fs::read_to_string(out.join("vdj_calls.tsv"))
        .expect("read vdj_calls.tsv");
    let (_header, rich_rows) = tsv_rows(&rich);
    assert_eq!(rich_rows.len(), 3, "expected exactly IGH, IGK and IGL calls");

    let expected = HashMap::from([
        ("IGH", ("Vdj", "Ighv1-64", "Ighd1-1", "Ighj3", "Igha", "18")),
        ("IGK", ("Vj", "Igkv6-15", "", "Igkj2", "Igkc", "382")),
        ("IGL", ("Vj", "Iglv3", "", "Iglj2", "", "1")),
    ]);

    for row in &rich_rows {
        assert_eq!(field(row, "cell", "one-cell VDJ fixture"), CELL);
        assert_eq!(
            field(row, "rustody_cell_id", "one-cell VDJ fixture"),
            "38637011",
            "wrong BD/Rustody positional cell id for fixture barcode {CELL}"
        );
        let chain = field(row, "chain", "one-cell VDJ fixture");
        let Some((stage, v, d, j, c, umis)) = expected.get(chain) else {
            panic!("unexpected chain in one-cell fixture: {chain}");
        };
        assert_eq!(field(row, "stage", chain), *stage, "wrong stage for {chain}");
        assert_eq!(field(row, "v", chain), *v, "wrong V call for {chain}");
        assert_eq!(field(row, "d", chain), *d, "wrong D call for {chain}");
        assert_eq!(field(row, "j", chain), *j, "wrong J call for {chain}");
        assert_eq!(field(row, "c", chain), *c, "wrong C call for {chain}");
        assert_eq!(
            field(row, "support_umis", chain),
            *umis,
            "wrong UMI support for {chain}; expected fixture call: stage={stage}, V={v}, D={d:?}, J={j}, C={c:?}, UMIs={umis}; actual row: {row:#?}"
        );
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
        assert_eq!(
            row["lumrik_rustody_cell_id"],
            "38637011",
            "wrong AIRR BD/Rustody positional cell id for fixture barcode {CELL}"
        );
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

        // AIRR sequence/junction amino-acid annotations are only meaningful
        // when the conserved V and J anchors establish a receptor frame.
        // Do not manufacture translations for unresolved/truncated calls.
        for field in ["lumrik_v_cys_anchor", "lumrik_j_anchor"] {
            assert!(
                matches!(row[field], "T" | "F"),
                "{field} must be AIRR T/F for {chain}, got {:?}",
                row[field]
            );
        }
        for field in ["productive", "vj_in_frame", "stop_codon"] {
            assert!(
                row[field].is_empty() || matches!(row[field], "T" | "F"),
                "{field} must be empty (unknown) or AIRR T/F for {chain}, got {:?}",
                row[field]
            );
        }

        let has_v_anchor = row["lumrik_v_cys_anchor"] == "T";
        let has_j_anchor = row["lumrik_j_anchor"] == "T";
        if has_v_anchor && has_j_anchor {
            assert!(!row["sequence_aa"].is_empty(), "missing sequence_aa for anchored {chain}");
            assert!(!row["junction"].is_empty(), "missing junction for anchored {chain}");
            assert!(!row["junction_aa"].is_empty(), "missing junction_aa for anchored {chain}");
            assert!(!row["cdr3"].is_empty(), "missing cdr3 for anchored {chain}");
            assert!(!row["cdr3_aa"].is_empty(), "missing cdr3_aa for anchored {chain}");

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
        } else {
            assert!(row["junction"].is_empty(), "unanchored {chain} must not fabricate junction");
            assert!(row["junction_aa"].is_empty(), "unanchored {chain} must not fabricate junction_aa");
            assert!(row["cdr3"].is_empty(), "unanchored {chain} must not fabricate cdr3");
            assert!(row["cdr3_aa"].is_empty(), "unanchored {chain} must not fabricate cdr3_aa");
        }
    }

    fs::remove_dir_all(&out).expect("remove successful test output");
}
