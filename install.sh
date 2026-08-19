#!/usr/bin/env bash
set -e

cargo install --path crates/bam_tide --force
cargo install --path crates/sc_primer --force
cargo install --path crates/sc-mapper --force
cargo install --path crates/gtf_splice_index --force
cargo install --path crates/fast_tag_mapper --force
cargo install --path crates/read-tag-table --force
cargo install --path crates/snp-index --force
cargo install --path crates/sc-beacon --force

