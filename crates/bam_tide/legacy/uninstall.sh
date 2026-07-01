#!/usr/bin/env bash
set -e

cargo uninstall bam_tide
DATA_DIR="$(cargo run --quiet --bin show_data_dir)"
echo "Removing data directory: $DATA_DIR"
rm -rf "$DATA_DIR"
echo "Uninstalled."
