#!/bin/bash

# Optional first argument: chromosomal region (e.g. "chr1:1000-2000")
region=$1
cause=$2

if [ -z "$region" ]; then
    cargo build -r
    target/release/bam2bigwig -b testData/test.bam -o testData/output/test.rust.bw
    bigwigtobedgraph testData/output/test.rust.bw testData/output/test.rust.bedgraph
    grep -v "0.0\$" testData/output/test.rust.bedgraph > testData/output/test.rust.bedgraph.no0
    ls testData/output/test.*0
    diff testData/output/test.*0 -y
else
    cargo build
    samtools view -hb testData/test.bam > testData/subset.bam
    echo $cause
    target/debug/bam2bedgraph -b testData/subset.bam -o testData/output/subset.rust.bedgraph
    echo $cause
    cat testData/output/subset.rust.bedgraph | grep -v "0$"
fi
