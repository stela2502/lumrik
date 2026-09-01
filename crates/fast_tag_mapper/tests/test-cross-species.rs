use fast_tag_mapper::{FastTagMapper, MapStatus, HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS};
use mapping_info::MappingInfo;

fn info() -> MappingInfo {
    MappingInfo::new(None, 0.0, 0)
}

#[test]
fn mouse_mapper_accepts_mouse_and_rejects_human() {
    let mapper = FastTagMapper::mouse_samples();

    // Every mouse sample maps to itself.
    for (idx, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            Some((idx + 1) as u64),
            "Mouse SampleTag{:02} failed",
            idx + 1
        );
    }

    // No human sample maps into mouse space.
    for (idx, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            None,
            "Human SampleTag{:02} cross-reacts with mouse mapper",
            idx + 1
        );
    }
}

#[test]
fn human_mapper_accepts_human_and_rejects_mouse() {
    let mapper = FastTagMapper::human_samples();

    // Every human sample maps to itself.
    for (idx, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            Some((idx + 1) as u64),
            "Human SampleTag{:02} failed",
            idx + 1
        );
    }

    // No mouse sample maps into human space.
    for (idx, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            None,
            "Mouse SampleTag{:02} cross-reacts with human mapper",
            idx + 1
        );
    }
}
#[test]
fn report_max_cross_species_vote_count() {
    let mapper = FastTagMapper::mouse_samples();

    let mut max_votes = 0;

    for seq in HUMAN_SAMPLE_TAGS {
        let mut mi = info();

        if let MapStatus::Tie { hits, .. } = mapper.map_status(seq, &mut mi) {
            max_votes = max_votes.max(hits);
        }

        if let MapStatus::Hit { hits, .. } = mapper.map_status(seq, &mut mi) {
            max_votes = max_votes.max(hits);
        }
    }

    println!("Maximum human->mouse cross reaction vote count: {max_votes}");
}

#[test]
fn report_cross_species_vote_counts() {
    let mapper = FastTagMapper::mouse_samples().with_min_hits(1);

    let mut max_hits = 0;

    for seq in HUMAN_SAMPLE_TAGS {
        let mut mi = info();

        if let MapStatus::Hit { hits, .. } = mapper.map_status(seq, &mut mi) {
            max_hits = max_hits.max(hits);
        }

        if let MapStatus::Tie { hits, .. } = mapper.map_status(seq, &mut mi) {
            max_hits = max_hits.max(hits);
        }
    }

    println!("Max human->mouse cross reaction = {max_hits}");
}
