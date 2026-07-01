# fast_tag_mapper

`fast_tag_mapper` is a small Rust crate for fast identification of short DNA tag sequences from FASTQ/BAM reads.

The original use case is BD Rhapsody sample-tag detection inside `bam_tide`, but the core model is general enough for FASTA/FASTQ-derived tag sets:

- build a table from known tag sequences
- encode every 8 bp window as a 2-bit `u16`
- invalidate non-unique 8-mers
- map read sequences by voting for a unique `(feature, start)` candidate
- return only the matched feature id when the internal `min_hits` threshold is surpassed

The hot-path API is intentionally simple:

```rust
mapper.map_feature_id(seq, &mut mapping_info) -> Option<u64>
```

`Some(feature_id)` is ready to use as a `Scdata` feature id. `None` means no unique, strong enough hit.

## Conceptual model

The crate uses two different concepts that should not be mixed up.

### FeatureEntry

A `FeatureEntry` is the biological or user-facing tag identity.

Examples:

- `SampleTag01_mm`
- `SampleTag07_hs`
- a clean FASTA record name

It owns:

- the `u64` feature id used in `Scdata`
- the clean feature name
- the 10x feature type, usually `Antibody Capture`

### TagEntry

A `TagEntry` is internal lookup-table state.

It represents one 8 bp 2-bit k-mer position inside one feature sequence:

- feature index
- position of that 8-mer inside the feature/tag sequence

Users normally should not interact with `TagEntry`.

## Dependencies

Add the crate and its ecosystem dependencies to your `Cargo.toml`.

```toml
[dependencies]
fast_tag_mapper = { git = "https://github.com/stela2502/fast_tag_mapper" }
mapping_info = { git = "https://github.com/stela2502/mapping_info" }
scdata = { git = "https://github.com/stela2502/scdata" }
```

Internally, `fast_tag_mapper` uses:

```toml
stela_int_to_str = { git = "https://github.com/stela2502/int_to_str" }
```

The crate depends on `stela_int_to_str` because the original `int_to_str` package name was not available on crates.io.

## Built-in BD sample tags

The crate includes the 12 human and 12 mouse BD sample tags.

```rust
use fast_tag_mapper::FastTagMapper;
use mapping_info::MappingInfo;

let mapper = FastTagMapper::mouse_samples().with_min_hits(4);
let mut info = MappingInfo::new(None, 0.0, 0);

let read = b"GTTGTCAAGATGCTACCGTTCAGAGACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC";

let feature_id = mapper.map_feature_id(read, &mut info);

assert_eq!(feature_id, Some(7));
```

For human tags:

```rust
let mapper = FastTagMapper::human_samples().with_min_hits(4);
```

## Why `min_hits = 4`?

The default/minimum safe BD sample-tag threshold was chosen from cross-species tests.

The important regression test is:

- a mouse mapper maps all 12 mouse tags to ids `1..=12`
- a mouse mapper rejects all 12 human tags
- a human mapper maps all 12 human tags to ids `1..=12`
- a human mapper rejects all 12 mouse tags

With `min_hits = 2`, exact cross-species hits occur:

- human `SampleTag10_hs` can cross-detect as mouse `SampleTag06_mm`
- mouse `SampleTag06_mm` can cross-detect as human `SampleTag10_hs`

With `min_hits = 4`, those exact cross-species detections are rejected while exact same-species tags still have many independent 8-mer votes.

For noisy real reads, a higher threshold can be used:

```rust
let mapper = FastTagMapper::mouse_samples().with_min_hits(8);
```

The best threshold depends on how much sequencing error and truncation you expect.

## Main API

### Map a read to a feature id

```rust
pub fn map_feature_id(
    &self,
    seq: &[u8],
    mapping: &mut MappingInfo,
) -> Option<u64>
```

Example:

```rust
if let Some(feature_id) = mapper.map_feature_id(record.seq(), &mut mapping_info) {
    // feature_id can be used directly in Scdata
}
```

The mapper updates `MappingInfo` internally with hit/no-hit/tie counters and timing.

### Debug mapping status

For debugging, use:

```rust
pub fn map_status(
    &self,
    seq: &[u8],
    mapping: &mut MappingInfo,
) -> MapStatus
```

Example:

```rust
use fast_tag_mapper::MapStatus;

match mapper.map_status(read, &mut mapping_info) {
    MapStatus::Hit {
        feature_id,
        feature_index,
        start,
        hits,
    } => {
        eprintln!(
            "feature_id={feature_id} feature_index={feature_index} start={start} hits={hits}"
        );
    }
    MapStatus::NoHit => {
        eprintln!("no usable tag hit");
    }
    MapStatus::Tie { hits, feature_ids } => {
        eprintln!("tie with {hits} hits: {feature_ids:?}");
    }
}
```

The hot path should still use `map_feature_id`.

## Integration with Scdata

Because `map_feature_id` returns `Option<u64>`, integration into a count table should be direct.

Pseudo-example:

```rust
if let Some(sample_feature_id) = sample_mapper.map_feature_id(&r2.seq, &mut stats) {
    feature_tag_counts.increment(cell_id, sample_feature_id, 1);
}
```

The exact method name depends on your `Scdata` API, but the key point is that no name lookup is needed in the hot path.

## FeatureIndex integration

`FastTagFeatureIndex` exposes mapper features to `scdata::FeatureIndex`.

```rust
use fast_tag_mapper::{FastTagFeatureIndex, FastTagMapper};
use scdata::FeatureIndex;

let mapper = FastTagMapper::mouse_samples().with_min_hits(4);
let index = FastTagFeatureIndex::new(&mapper);

assert_eq!(index.feature_id("SampleTag07_mm"), Some(7));
assert_eq!(index.feature_name(7), "SampleTag07_mm");
assert_eq!(
    index.to_10x_feature_line(7),
    "SampleTag07_mm\tSampleTag07_mm\tAntibody Capture"
);
```

`FastTagFeatureIndex` is built from `FeatureEntry`, not from internal lookup-table `TagEntry`.

## Adding custom tags

To build a mapper manually, add feature sequences with clean names and stable ids.

```rust
use fast_tag_mapper::{FastTagMapper, FeatureEntry};

let mut mapper = FastTagMapper::new().with_min_hits(4);

mapper.add_feature(
    b"ACGTACGTACGTACGTACGT",
    FeatureEntry::new(1, "MyTag01", "Antibody Capture"),
);

mapper.add_feature(
    b"TGCATGCATGCATGCATGCA",
    FeatureEntry::new(2, "MyTag02", "Antibody Capture"),
);
```

The returned value from `map_feature_id` will be the `FeatureEntry.id`.

## Encoded-position API

Most users should use `map_feature_id(seq, mapping)`.

If you already encoded read windows elsewhere, use the encoded-position API:

```rust
mapper.map_encoded_positions_feature_id(encoded_positions, &mut mapping_info)
```

Important: encoded positions must be physical base-pair positions in the query read.

Correct:

```text
(24, encoded_8mer)
```

Incorrect:

```text
(3, encoded_8mer)
```

The second form is only a compressed chunk index and will infer the wrong start position.

## Tests worth keeping

The most important tests are not synthetic "does the function run" tests. They should validate the tag set.

### Same-species positive and cross-species negative test

```rust
use fast_tag_mapper::{FastTagMapper, HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS};
use mapping_info::MappingInfo;

fn info() -> MappingInfo {
    MappingInfo::new(None, 0.0, 0)
}

#[test]
fn mouse_mapper_accepts_mouse_and_rejects_human() {
    let mapper = FastTagMapper::mouse_samples().with_min_hits(4);

    for (idx, seq) in MOUSE_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            Some((idx + 1) as u64),
            "Mouse SampleTag{:02} failed",
            idx + 1
        );
    }

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
    let mapper = FastTagMapper::human_samples().with_min_hits(4);

    for (idx, seq) in HUMAN_SAMPLE_TAGS.iter().enumerate() {
        let mut mi = info();

        assert_eq!(
            mapper.map_feature_id(seq, &mut mi),
            Some((idx + 1) as u64),
            "Human SampleTag{:02} failed",
            idx + 1
        );
    }

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
```

### Real BD sample-read test

```rust
#[test]
fn real_mouse_sample_read_maps_to_sampletag07() {
    let mapper = FastTagMapper::mouse_samples().with_min_hits(4);
    let mut mi = MappingInfo::new(None, 0.0, 0);

    let read = b"GTTGTCAAGATGCTACCGTTCAGAGACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC";

    assert_eq!(mapper.map_feature_id(read, &mut mi), Some(7));
}
```

## Design notes

`fast_tag_mapper` deliberately avoids being a general-purpose aligner.

It is optimized for short, known tag sequences where repeated exact 8-mer evidence is enough to identify a clean feature id.

The mapper does not try to rescue ambiguous hits. If the best vote is tied, or if the number of hits is below `min_hits`, the hot API returns `None`.

This behavior is intentional for sample-tag counting, where a wrong count is worse than a missed count.
