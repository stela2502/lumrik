#!/usr/bin/env bash
set -e

echo "Installing bam_tide..."
cargo install --path .

DATA_DIR="$(cargo run --quiet --bin show_data_dir)"
echo "Installing whitelist_map.bin to $DATA_DIR"
mkdir -p "$DATA_DIR"
cp data/whitelist_map.bin "$DATA_DIR"
echo "Done."