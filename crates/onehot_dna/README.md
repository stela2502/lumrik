# onehot_dna

Small fixed-length one-hot DNA encoding for fast barcode and primer matching in Rust.

This crate is meant for short sequences where all candidates have the same length, for example single-cell barcode blocks, sample tags, primer fragments, or other small DNA tokens that need fast exact or fuzzy matching.

The central type is:

```rust
OneHot<const N: usize>
```

It stores `N` DNA bases in a single `u128`, using four bits per base. Therefore `N <= 32`.

## Encoding

Each DNA base is encoded as a four-bit one-hot mask:

| Base | Bits |
|---|---:|
| `A` | `0001` |
| `C` | `0010` |
| `G` | `0100` |
| `T` | `1000` |
| anything else, including `N` | `0000` |

Unknown bases are encoded as zero. They match nothing and therefore count as mismatches. This is intentional for barcode correction: an ambiguous base should not accidentally rescue a hit.

## Install

Add the crate directly from GitHub:

```toml
[dependencies]
onehot_dna = { git = "https://github.com/stela2502/onehot_dna" }
```

Or, for a fixed revision:

```toml
[dependencies]
onehot_dna = { git = "https://github.com/stela2502/onehot_dna", rev = "<commit>" }
```

Then import the types you need:

```rust
use onehot_dna::{OneHot, OneHot9, OneHotSet};
```

## Basic usage

Encode a fixed-length DNA sequence:

```rust
use onehot_dna::OneHot;

let barcode = OneHot::<9>::from_bytes(b"ACGTACGTN")?;

assert_eq!(barcode.to_dna_string(), "ACGTACGTN");
println!("{barcode}");
# Ok::<(), Box<dyn std::error::Error>>(())
```

The length is part of the type. A `OneHot<9>` and a `OneHot<16>` are different types, which helps avoid mixing incompatible barcode classes by accident.

```rust
use onehot_dna::OneHot;

let cell_block = OneHot::<9>::from_bytes(b"ACGTACGTA")?;
let umi_like   = OneHot::<8>::from_bytes(b"ACGTACGT")?;

// cell_block and umi_like cannot be compared directly because their lengths differ.
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Extracting a window from a read

Use `from_window` when the barcode/primer segment is embedded in a larger read:

```rust
use onehot_dna::OneHot;

let read = b"XXACGTACGTNYY";
let block = OneHot::<9>::from_window(read, 2)?;

assert_eq!(block.to_dna_string(), "ACGTACGTN");
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Exact and fuzzy matching

`mismatches` returns the barcode-style Hamming distance between two encoded sequences of the same length:

```rust
use onehot_dna::OneHot;

let a = OneHot::<9>::from_bytes(b"AAAAAAAAA")?;
let b = OneHot::<9>::from_bytes(b"AAAAAAAAC")?;
let n = OneHot::<9>::from_bytes(b"AAAAAAAAN")?;

assert_eq!(a.mismatches(b), 1);
assert_eq!(a.mismatches(n), 1); // N/unknown matches nothing

assert!(a.within(b, 1));
assert!(!a.within(b, 0));
# Ok::<(), Box<dyn std::error::Error>>(())
```

Exact matching is simple bit equality:

```rust
use onehot_dna::OneHot;

let a = OneHot::<9>::from_bytes(b"ACGTACGTA")?;
let b = OneHot::<9>::from_bytes(b"ACGTACGTA")?;

assert!(a.exact_match(b));
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Candidate tables with `OneHotSet`

`OneHotSet<N>` is a small candidate table for fixed-length barcodes or primers.

```rust
use onehot_dna::{OneHot9, OneHotSet};

let candidates = OneHotSet::<9>::from_sequences(&[
    b"AAAAAAAAA".as_slice(),
    b"CCCCCCCCC".as_slice(),
    b"GGGGGGGGG".as_slice(),
])?;

let query = OneHot9::from_bytes(b"AAAAAAAAC")?;
let hit = candidates.best_match(&query, 1);

assert_eq!(hit, Some((0, 1)));
# Ok::<(), Box<dyn std::error::Error>>(())
```

`best_match` returns:

```rust
Option<(usize, u32)>
```

where the tuple is `(candidate_index, mismatch_distance)`.

It returns `None` when:

- no candidate is within `max_mismatches`
- the best hit is tied between multiple candidates

This is useful for barcode correction where a non-unique correction should be rejected rather than guessed.

Example tie rejection:

```rust
use onehot_dna::{OneHot9, OneHotSet};

let candidates = OneHotSet::<9>::from_sequences(&[
    b"AAAAAAAAA".as_slice(),
    b"AAAAAAAAC".as_slice(),
])?;

let query = OneHot9::from_bytes(b"AAAAAAAAN")?;

// The query is one mismatch away from both candidates, so the match is ambiguous.
assert_eq!(candidates.best_match(&query, 1), None);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Estimating a safe correction radius

`OneHotSet` can report the minimum pairwise distance between all candidates:

```rust
use onehot_dna::OneHotSet;

let candidates = OneHotSet::<9>::from_sequences(&[
    b"AAAAAAAAA".as_slice(),
    b"AAAAAAAAC".as_slice(),
    b"CCCCCCCCC".as_slice(),
])?;

let min_distance = candidates.min_pairwise_mismatches();
let radius = candidates.correction_radius();

assert_eq!(min_distance, Some(1));
assert_eq!(radius, Some(0));
# Ok::<(), Box<dyn std::error::Error>>(())
```

The correction radius is computed as:

```text
(min_pairwise_distance - 1) / 2
```

That is the conservative Hamming-style radius where two valid candidates should not claim the same observed sequence.

For real barcode whitelists you may still want to set your own `max_mismatches`, because empirical sequencing error and whitelist design often matter more than the theoretical global radius.

## BD Rhapsody 9 bp blocks

The crate includes a convenience alias:

```rust
pub type OneHot9 = OneHot<9>;
```

This is useful for BD Rhapsody-style C1/C2/C3 barcode blocks:

```rust
use onehot_dna::{OneHot9, OneHotSet};

let c1_table = OneHotSet::<9>::from_sequences(&[
    b"ACGTACGTA".as_slice(),
    b"TGCATGCAT".as_slice(),
])?;

let observed_c1 = OneHot9::from_bytes(b"ACGTACGTT")?;

if let Some((idx, dist)) = c1_table.best_match(&observed_c1, 1) {
    println!("corrected C1 block to candidate {idx} with {dist} mismatch(es)");
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Raw bits

If you need compact storage or precomputed constants, use `bits` and `from_bits`:

```rust
use onehot_dna::OneHot;

let x = OneHot::<9>::from_bytes(b"ACGTACGTA")?;
let packed: u128 = x.bits();

let restored = OneHot::<9>::from_bits(packed);
assert_eq!(x, restored);
# Ok::<(), Box<dyn std::error::Error>>(())
```

`from_bits` does not validate that every four-bit nibble contains exactly one legal base bit. It is intended for trusted packed values.

## Error handling

Construction can fail when the sequence length does not match `N`, or when `N > 32`:

```rust
use onehot_dna::{OneHot, OneHotError};

let err = OneHot::<9>::from_bytes(b"ACGT").unwrap_err();

assert_eq!(
    err,
    OneHotError::WrongLength {
        expected: 9,
        observed: 4,
    }
);
```

## Design notes

This crate is intentionally small:

- no allocation for individual `OneHot<N>` values
- compact `u128` storage
- const-generic fixed length
- simple mismatch counting
- unknown bases treated as non-matching
- deterministic unique-best-hit correction via `OneHotSet::best_match`

It is best suited for short fixed-length DNA tokens, not for general sequence alignment.

## Run tests

```bash
cargo test
```

## License

Add your license information here.
