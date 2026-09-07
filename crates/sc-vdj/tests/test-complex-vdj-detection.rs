use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::Writer;
use sc_vdj::output::{write_mapping_info, ReportWriter};
use sc_vdj::{
    Chain, NelruneIdentityResolver, SegmentKind, Strand, VdjIndex, VdjRunner, VdjRunnerConfig,
    VdjSegment,
};
use std::collections::HashMap;
use std::fs;

fn seg(name: &str, kind: SegmentKind, start: u32, seq: &[u8]) -> VdjSegment {
    VdjSegment {
        id: 0,
        name: name.into(),
        transcript_id: format!("{name}.tx"),
        gene_id: name.into(),
        chain: Chain::Igh,
        kind,
        chromosome: "chr12".into(),
        start,
        end: start + seq.len() as u32,
        strand: Strand::Plus,
        sequence: seq.to_vec(),
    }
}
fn rc(x: &[u8]) -> Vec<u8> {
    x.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[test]
fn complex_fragmented_vdj_with_pn_d_and_constant_is_detected_from_bam() {
    // Deliberately non-repetitive germline pieces.  The biological molecule
    // contains coding-end deletions, V/J P additions, independent N1/N2,
    // a short retained D, and a C continuation.
    let v = b"ACGTCAGTACCTGACCGTATGCACTGATCGTAACGTAGCA";
    let d = b"GATCCGTAGTCA";
    let j = b"TGCATACGGTACCTAGGCTAACGTGACT";
    let c = b"CCGATGTCAGACTGACCTGGTACCGTACGATGCTAGCTAG";
    let index = VdjIndex::from_segments(vec![
        seg("IghvComplex", SegmentKind::V, 100, v),
        seg("IghdComplex", SegmentKind::D, 300, d),
        seg("IghjComplex", SegmentKind::J, 500, j),
        seg("Igha", SegmentKind::C, 700, c),
    ])
    .unwrap();

    let v_del = 3usize;
    let d_del5 = 2usize;
    let d_del3 = 2usize;
    let j_del5 = 3usize;
    let vr = &v[..v.len() - v_del];
    let dr = &d[d_del5..d.len() - d_del3];
    let jr = &j[j_del5..];
    let pv = rc(&vr[vr.len() - 2..]);
    let pj = rc(&jr[..2]);
    let n1 = b"AATGC";
    let n2 = b"GTTGA";
    let mut receptor = Vec::new();
    receptor.extend_from_slice(vr);
    receptor.extend_from_slice(&pv);
    receptor.extend_from_slice(n1);
    receptor.extend_from_slice(dr);
    receptor.extend_from_slice(n2);
    receptor.extend_from_slice(&pj);
    receptor.extend_from_slice(jr);
    let rearr_len = receptor.len();
    receptor.extend_from_slice(&c[..30]);

    // Four overlapping fragments: no BAM record contains the whole VDJ event.
    // Their genomic mapper coordinates nominate V, D, J and C respectively;
    // sequence overlap must reconstruct the receptor before junction inference.
    let ranges = [
        (0, 52, 110u32),
        (34, 72, 302u32),
        (55, 96, 503u32),
        (80, receptor.len(), 702u32),
    ];
    let tmp = tempfile::tempdir().unwrap();
    let bam_path = tmp.path().join("complex.bam");
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr12")
            .push_tag(b"LN", 2000),
    );
    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    for (i, (a, b, pos)) in ranges.into_iter().enumerate() {
        let seq = &receptor[a..b];
        let qual = vec![40u8; seq.len()];
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        let mut rec = bam::Record::new();
        rec.set(format!("fragment{i}").as_bytes(), Some(&cigar), seq, &qual);
        // Synthetic mapper output: bam::Record::new() starts as unmapped.
        rec.set_flags(0);
        rec.set_tid(0);
        rec.set_pos(pos as i64);
        rec.set_mapq(60);
        rec.push_aux(b"CB", Aux::String("ACGTACGTACGTACGT"))
            .unwrap();
        writer.write(&rec).unwrap();
    }
    drop(writer);

    let mut runner = VdjRunner::new(
        index,
        VdjRunnerConfig {
            min_sequence_overlap: 14,
        },
    );
    assert_eq!(
        runner
            .read_bam(&bam_path, &NelruneIdentityResolver)
            .unwrap(),
        4
    );
    let calls = runner.identify();
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].1.len(), 1);
    let r = &calls[0].1[0];
    assert_eq!(r.chain, Chain::Igh);
    assert_eq!(runner.index.segment(r.v).unwrap().name, "IghvComplex");
    assert_eq!(
        runner.index.segment(r.d.unwrap()).unwrap().name,
        "IghdComplex"
    );
    assert_eq!(runner.index.segment(r.j).unwrap().name, "IghjComplex");
    assert_eq!(
        runner
            .index
            .segment(r.constant.as_ref().unwrap().segment)
            .unwrap()
            .name,
        "Igha"
    );
    assert_eq!(r.junction.v_del_3, v_del as u16);
    assert_eq!(r.junction.p_v3, pv);
    assert_eq!(r.junction.n1, n1);
    assert!(r.junction.p_d5.is_empty());
    assert_eq!(r.junction.d_del_5, Some(d_del5 as u16));
    assert_eq!(r.junction.d_retained, dr);
    assert_eq!(r.junction.d_del_3, Some(d_del3 as u16));
    assert!(r.junction.p_d3.is_empty());
    assert_eq!(r.junction.n2, n2);
    assert_eq!(r.junction.p_j5, pj);
    assert_eq!(r.junction.j_del_5, j_del5 as u16);
    assert_eq!(r.observed_rearrangement, &receptor[..rearr_len]);
    assert!(r.observed_receptor_sequence.ends_with(&c[..30]));

    // The compact id is structural: it round-trips the same measurements and
    // excludes constant-region identity even though C was assembled above.
    let decoded = r.stable_id.decode(&runner.index).unwrap();
    assert_eq!(decoded.v, "IghvComplex");
    assert_eq!(decoded.d.as_deref(), Some("IghdComplex"));
    assert_eq!(decoded.j, "IghjComplex");
    assert_eq!(decoded.v_del_3, v_del as u16);
    assert_eq!(decoded.d_del_5, Some(d_del5 as u16));
    assert_eq!(decoded.d_retained_len, Some(dr.len() as u16));
    assert_eq!(decoded.d_del_3, Some(d_del3 as u16));
    assert_eq!(decoded.j_del_5, j_del5 as u16);
    assert!(
        r.stable_id.to_string().len() < 32,
        "compact recombination id unexpectedly bloated: {}",
        r.stable_id
    );

    // Production output contract: the clean caller must still emit the rich
    // Lumrik tables, AIRR table, mapping info and optional sequence FASTAs.
    let out = tmp.path().join("out");
    let mut report = ReportWriter::create(&out, true).unwrap();
    report
        .write_cell("ACGTACGTACGTACGT", &calls[0].1, &runner.index)
        .unwrap();
    report.finish().unwrap();
    let by_chain = HashMap::from([(String::from("IGH"), 1usize)]);
    write_mapping_info(out.join("vdj-mapping-info.txt"), 4, 1, 1, &by_chain).unwrap();
    for name in [
        "vdj_calls.tsv",
        "vdj_receptors.tsv",
        "airr_rearrangements.tsv",
        "vdj-mapping-info.txt",
        "vdj_observed.fasta",
        "vdj_naive.fasta",
    ] {
        let path = out.join(name);
        assert!(
            path.is_file(),
            "missing production output {}",
            path.display()
        );
        assert!(
            fs::metadata(&path).unwrap().len() > 0,
            "empty production output {}",
            path.display()
        );
    }
    let rich = fs::read_to_string(out.join("vdj_calls.tsv")).unwrap();
    assert!(rich.contains("IghvComplex"));
    assert!(rich.contains("IghdComplex"));
    assert!(rich.contains("IghjComplex"));
    assert!(rich.contains("Igha"));
}
