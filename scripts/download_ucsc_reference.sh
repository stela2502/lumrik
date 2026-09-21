#!/usr/bin/env bash
set -euo pipefail

BASE="https://hgdownload.soe.ucsc.edu"

GENOME=""
OUTPATH=""

usage() {
    cat <<EOF
Usage:
  $0 --genome <UCSC genome> --outpath <directory>

Examples:
  $0 --genome mm39 --outpath /data1/UCSC
  $0 --genome hg38 --outpath /data1/UCSC
  $0 --genome rn7  --outpath /data1/UCSC

Output:
  <outpath>/<genome>/
      Genome/
      Genes/
      Protein/
      download-report.txt
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --outpath)
            OUTPATH="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$GENOME" ]]; then
    echo "ERROR: --genome is required" >&2
    exit 1
fi

if [[ -z "$OUTPATH" ]]; then
    echo "ERROR: --outpath is required" >&2
    exit 1
fi

ROOT="$OUTPATH/$GENOME"
GENOME_DIR="$ROOT/Genome"
GENES_DIR="$ROOT/Genes"
PROTEIN_DIR="$ROOT/Protein"
REPORT="$ROOT/download-report.txt"

mkdir -p "$GENOME_DIR" "$GENES_DIR" "$PROTEIN_DIR"

: > "$REPORT"

log() {
    echo "$*"
    echo "$*" >> "$REPORT"
}

exists() {
    wget --quiet --spider "$1"
}

download() {
    local url="$1"
    local output="$2"

    log "  source: $url"
    log "  target: $output"

    wget \
        --continue \
        --output-document="$output" \
        "$url"
}

probe_and_download() {
    local label="$1"
    local url="$2"
    local output="$3"
    local requirement="$4"

    log ""
    log "Checking $label:"
    log "  $url"

    if exists "$url"; then
        log "  AVAILABLE"
        download "$url" "$output"
    else
        if [[ "$requirement" == "required" ]]; then
            log "  WARNING: REQUIRED resource unavailable"
        else
            log "  INFO: optional resource unavailable"
        fi
    fi
}

log "UCSC reference download"
log "======================="
log "Genome:  $GENOME"
log "Output:  $ROOT"
log "UCSC:    $BASE"
log ""

#
# Genome
#

log "============================================================"
log " Genome"
log "============================================================"

BIGZIPS="$BASE/goldenPath/$GENOME/bigZips"
LATEST="$BIGZIPS/latest"

find_genome_resource() {
    local filename="$1"

    if exists "$LATEST/$filename"; then
        echo "$LATEST/$filename"
        return 0
    fi

    if exists "$BIGZIPS/$filename"; then
        echo "$BIGZIPS/$filename"
        return 0
    fi

    return 1
}

download_genome_resource() {
    local label="$1"
    local filename="$2"
    local requirement="$3"

    log ""
    log "Checking $label:"
    log "  $LATEST/$filename"
    log "  $BIGZIPS/$filename"

    local url=""

    if url="$(find_genome_resource "$filename")"; then
        log "  AVAILABLE"
        download "$url" "$GENOME_DIR/$filename"
    else
        if [[ "$requirement" == "required" ]]; then
            log "  WARNING: REQUIRED resource unavailable"
        else
            log "  INFO: optional resource unavailable"
        fi
    fi
}

download_genome_resource \
    "2bit genome" \
    "$GENOME.2bit" \
    required

download_genome_resource \
    "chromosome sizes" \
    "$GENOME.chrom.sizes" \
    required

download_genome_resource \
    "chromosome aliases" \
    "$GENOME.chromAlias.txt" \
    optional

#
# Genes
#

log ""
log "============================================================"
log " Genes / transcripts"
log "============================================================"

GENES="$BIGZIPS/genes"

probe_and_download \
    "UCSC knownGene GTF" \
    "$GENES/$GENOME.knownGene.gtf.gz" \
    "$GENES_DIR/$GENOME.knownGene.gtf.gz" \
    required

probe_and_download \
    "NCBI RefSeq GTF" \
    "$GENES/$GENOME.ncbiRefSeq.gtf.gz" \
    "$GENES_DIR/$GENOME.ncbiRefSeq.gtf.gz" \
    optional

probe_and_download \
    "refGene GTF" \
    "$GENES/refGene.gtf.gz" \
    "$GENES_DIR/refGene.gtf.gz" \
    optional

#
# Protein / UniProt
#

log ""
log "============================================================"
log " Protein / UniProt"
log "============================================================"

UNIPROT_ROOT="$BASE/goldenPath/archive/$GENOME/uniprot"

log ""
log "Discovering UCSC UniProt releases:"
log "  $UNIPROT_ROOT/"

UNIPROT_RELEASE="$(
    wget -qO- "$UNIPROT_ROOT/" 2>/dev/null |
        grep -oE 'href="[0-9]{4}_[0-9]{2}/"' |
        sed -E 's/^href="//; s:/"$::' |
        sort -V |
        tail -n 1 || true
)"

if [[ -n "$UNIPROT_RELEASE" ]]; then

    UNIPROT="$UNIPROT_ROOT/$UNIPROT_RELEASE"

    log "  newest release: $UNIPROT_RELEASE"
    log "  using: $UNIPROT/"
    log ""

    probe_and_download \
        "UniProt version information" \
        "$UNIPROT/version.txt" \
        "$PROTEIN_DIR/version.txt" \
        required

    probe_and_download \
        "UniProt track definitions" \
        "$UNIPROT/trackDb.txt" \
        "$PROTEIN_DIR/trackDb.txt" \
        required

    probe_and_download \
        "protein mapping information" \
        "$UNIPROT/protMapInfo.tsv" \
        "$PROTEIN_DIR/protMapInfo.tsv" \
        required

    probe_and_download \
        "UniProt lift information" \
        "$UNIPROT/liftInfo.json" \
        "$PROTEIN_DIR/liftInfo.json" \
        optional

    probe_and_download \
        "full protein sequences" \
        "$UNIPROT/unipFullSeq.bb" \
        "$PROTEIN_DIR/unipFullSeq.bb" \
        required

    probe_and_download \
        "protein domains" \
        "$UNIPROT/unipDomain.bb" \
        "$PROTEIN_DIR/unipDomain.bb" \
        optional

    probe_and_download \
        "transmembrane regions" \
        "$UNIPROT/unipLocTransMemb.bb" \
        "$PROTEIN_DIR/unipLocTransMemb.bb" \
        optional

    probe_and_download \
        "signal peptides" \
        "$UNIPROT/unipLocSignal.bb" \
        "$PROTEIN_DIR/unipLocSignal.bb" \
        optional

    probe_and_download \
        "cytoplasmic regions" \
        "$UNIPROT/unipLocCytopl.bb" \
        "$PROTEIN_DIR/unipLocCytopl.bb" \
        optional

    probe_and_download \
        "extracellular regions" \
        "$UNIPROT/unipLocExtra.bb" \
        "$PROTEIN_DIR/unipLocExtra.bb" \
        optional

    probe_and_download \
        "modified residues" \
        "$UNIPROT/unipModif.bb" \
        "$PROTEIN_DIR/unipModif.bb" \
        optional

    probe_and_download \
        "disulfide bonds" \
        "$UNIPROT/unipDisulfBond.bb" \
        "$PROTEIN_DIR/unipDisulfBond.bb" \
        optional

    probe_and_download \
        "protein repeats" \
        "$UNIPROT/unipRepeat.bb" \
        "$PROTEIN_DIR/unipRepeat.bb" \
        optional

    probe_and_download \
        "protein chains" \
        "$UNIPROT/unipChain.bb" \
        "$PROTEIN_DIR/unipChain.bb" \
        optional

    probe_and_download \
        "sequence conflicts" \
        "$UNIPROT/unipConflict.bb" \
        "$PROTEIN_DIR/unipConflict.bb" \
        optional

    probe_and_download \
        "regions of interest" \
        "$UNIPROT/unipInterest.bb" \
        "$PROTEIN_DIR/unipInterest.bb" \
        optional

    probe_and_download \
        "protein mutations" \
        "$UNIPROT/unipMut.bb" \
        "$PROTEIN_DIR/unipMut.bb" \
        optional

    probe_and_download \
        "other protein annotations" \
        "$UNIPROT/unipOther.bb" \
        "$PROTEIN_DIR/unipOther.bb" \
        optional

    probe_and_download \
        "splice annotations" \
        "$UNIPROT/unipSplice.bb" \
        "$PROTEIN_DIR/unipSplice.bb" \
        optional

    probe_and_download \
        "structure annotations" \
        "$UNIPROT/unipStruct.bb" \
        "$PROTEIN_DIR/unipStruct.bb" \
        optional

    probe_and_download \
        "Swiss-Prot genomic alignments" \
        "$UNIPROT/unipAliSwissprot.bb" \
        "$PROTEIN_DIR/unipAliSwissprot.bb" \
        optional

    probe_and_download \
        "TrEMBL genomic alignments" \
        "$UNIPROT/unipAliTrembl.bb" \
        "$PROTEIN_DIR/unipAliTrembl.bb" \
        optional

    probe_and_download \
        "UniProt-to-genome chain" \
        "$UNIPROT/unipToGenome.over.chain.gz" \
        "$PROTEIN_DIR/unipToGenome.over.chain.gz" \
        optional

    probe_and_download \
        "UniProt-to-genome PSL" \
        "$UNIPROT/unipToGenomeLift.psl.gz" \
        "$PROTEIN_DIR/unipToGenomeLift.psl.gz" \
        optional

    #
    # The transcript mapping bigBeds contain release-specific hashes in
    # their filenames, so discover them instead of guessing the hash.
    #

    log ""
    log "Discovering release-specific transcript/protein mapping bigBeds..."

    INDEX="$(
        wget -qO- "$UNIPROT/" 2>/dev/null || true
    )"

    for family in knownGene ensGene; do
        for source in swissprot trembl; do

            filename="$(
                printf '%s\n' "$INDEX" |
                    grep -oE "${family}_[0-9a-f]+\.${source}\.bb" |
                    sort -u |
                    tail -n 1 || true
            )"

            if [[ -n "$filename" ]]; then
                log ""
                log "Found $family / $source mapping:"
                log "  $filename"

                download \
                    "$UNIPROT/$filename" \
                    "$PROTEIN_DIR/$filename"
            else
                log ""
                log "INFO: no $family / $source mapping bigBed found"
            fi
        done
    done

else
    log "  WARNING: no UCSC UniProt release found for $GENOME"
    log "  Checked: $UNIPROT_ROOT/"
    log "  Genome and transcript resources remain usable."
fi

#
# Summary
#

log ""
log "============================================================"
log " Downloaded reference"
log "============================================================"
log ""

while IFS= read -r file; do
    size="$(du -h "$file" | cut -f1)"
    log "$size  $file"
done < <(
    find "$ROOT" -type f ! -path "$REPORT" | sort
)

log ""
log "Report:"
log "  $REPORT"
log ""
log "Done."
