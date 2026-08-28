use fast_tag_mapper::{
    encode_seq_positions_with_int_to_str, FastTagFeatureIndex, FastTagMapper, FeatureEntry,
    HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS, MapStatus, Slot,
};
use mapping_info::MappingInfo;
use scdata::FeatureIndex;

const BD_PREFIX: &[u8] = b"GTTGTCAAGATGCTACCGTTCAGAG";

fn info() -> MappingInfo {
    MappingInfo::new(None, 0.0, 0)
}

#[test]
fn builtin_features_are_real_12_each_and_named() {
    let human = FastTagMapper::human_samples();
    let mouse = FastTagMapper::mouse_samples();

    assert_eq!(HUMAN_SAMPLE_TAGS.len(), 12);
    assert_eq!(MOUSE_SAMPLE_TAGS.len(), 12);
    assert_eq!(human.feature_count(), 12);
    assert_eq!(mouse.feature_count(), 12);

    assert_eq!(human.feature(0).unwrap().id, 1);
    assert_eq!(human.feature(0).unwrap().name, "SampleTag01_hs");
    assert_eq!(human.feature(11).unwrap().id, 12);
    assert_eq!(human.feature(11).unwrap().name, "SampleTag12_hs");

    assert_eq!(mouse.feature(0).unwrap().id, 1);
    assert_eq!(mouse.feature(0).unwrap().name, "SampleTag01_mm");
    assert_eq!(mouse.feature(11).unwrap().id, 12);
    assert_eq!(mouse.feature(11).unwrap().name, "SampleTag12_mm");
}

#[test]
fn map_feature_id_returns_only_scdata_id_if_min_hits_surpassed() {
    let mapper = FastTagMapper::mouse_samples().with_min_hits(2);

    let mut read = BD_PREFIX.to_vec();
    read.extend_from_slice(MOUSE_SAMPLE_TAGS[6]);

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(&read, &mut mi), Some(7));
}

#[test]
fn too_high_min_hits_rejects_otherwise_correct_hit() {
    let mapper = FastTagMapper::mouse_samples().with_min_hits(10_000);

    let mut read = BD_PREFIX.to_vec();
    read.extend_from_slice(MOUSE_SAMPLE_TAGS[6]);

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(&read, &mut mi), None);
}

#[test]
fn real_mouse_sample_read_from_conversation_maps_to_sampletag07_feature_id() {
    let mapper = FastTagMapper::mouse_samples();
    let read = b"GTTGTCAAGATGCTACCGTTCAGAGACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC";

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(read, &mut mi), Some(7));

    let mut mi = info();
    match mapper.map_status(read, &mut mi) {
        MapStatus::Hit {
            feature_id,
            start,
            hits,
            ..
        } => {
            assert_eq!(feature_id, 7);
            assert_eq!(start, 25);
            assert!(hits >= 10);
        }
        other => panic!("expected hit, got {other:?}"),
    }
}

#[test]
fn real_human_sample_read_from_comment_maps_to_sampletag01_feature_id() {
    let mapper = FastTagMapper::human_samples();
    let read = b"GTTGTCAAGATGCTACCGTTCAGAGATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG";

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(read, &mut mi), Some(1));
}

#[test]
fn every_builtin_feature_maps_to_its_own_feature_id() {
    let human = FastTagMapper::human_samples();
    let mouse = FastTagMapper::mouse_samples();

    for (i, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();
        assert_eq!(human.map_feature_id(seq, &mut mi), Some((i + 1) as u64));
    }

    for (i, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();
        assert_eq!(mouse.map_feature_id(seq, &mut mi), Some((i + 1) as u64));
    }
}

#[test]
fn encoded_positions_use_physical_bp_positions_and_return_feature_id() {
    let mapper = FastTagMapper::mouse_samples();

    let mut read = BD_PREFIX.to_vec();
    read.extend_from_slice(MOUSE_SAMPLE_TAGS[6]);

    let encoded = encode_seq_positions_with_int_to_str(&read);

    let mut mi = info();
    assert_eq!(
        mapper.map_encoded_positions_feature_id(encoded, &mut mi),
        Some(7)
    );
}

#[test]
fn duplicate_8mers_are_invalidated_not_used() {
    let mut mapper = FastTagMapper::new().with_min_hits(1);

    mapper.add_feature(b"AAAAAAAACCCCCCCC", FeatureEntry::new(10, "a", "Antibody Capture"));
    mapper.add_feature(b"AAAAAAAAGGGGGGGG", FeatureEntry::new(20, "b", "Antibody Capture"));

    let encoded = fast_tag_mapper::fast_mapper::encode_8mer_with_int_to_str(b"AAAAAAAA")
        .expect("valid kmer");

    assert!(matches!(mapper.slot(encoded), Slot::Duplicate));

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(b"TTTTAAAAAAAATTTT", &mut mi), None);
}

#[test]
fn tie_returns_none_for_hot_api() {
    let mut mapper = FastTagMapper::new().with_min_hits(1);

    mapper.add_feature(b"ACGTACGT", FeatureEntry::new(101, "tag1", "Antibody Capture"));
    mapper.add_feature(b"TGCATGCA", FeatureEntry::new(202, "tag2", "Antibody Capture"));

    let read = b"NNNNACGTACGTNNNNTGCATGCA";

    let mut mi = info();
    assert_eq!(mapper.map_feature_id(read, &mut mi), None);

    let mut mi = info();
    match mapper.map_status(read, &mut mi) {
        MapStatus::Tie { hits, feature_ids } => {
            assert_eq!(hits, 1);
            assert_eq!(feature_ids, vec![101, 202]);
        }
        other => panic!("expected tie, got {other:?}"),
    }
}

#[test]
fn feature_index_uses_feature_entries_not_table_entries() {
    let mapper = FastTagMapper::mouse_samples();
    let index = FastTagFeatureIndex::new(&mapper);

    assert_eq!(index.feature_id("SampleTag07_mm"), Some(7));
    assert_eq!(index.feature_name(7), "SampleTag07_mm");
    assert_eq!(index.ordered_feature_ids(), (1_u64..=12).collect::<Vec<_>>());
    assert_eq!(
        index.to_10x_feature_line(7),
        "SampleTag07_mm\tSampleTag07_mm\tbd_sample_mouse"
    );
}
