#!/bin/bash
# Build the combined host-symbiont reference used to map Panama reads.
#
# Usage:
#   bash 00_build_combined_reference.sh HOST.fna SYMBIONT.fa [OUTPUT.fa]
#
# HOST.fna is the Acropora cervicornis jaAcrCerv1.1 assembly
# (NCBI accession GCA_964034985.1). SYMBIONT.fa contains the S. fitti
# reference sequences used in the study.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 HOST.fna SYMBIONT.fa [OUTPUT.fa]" >&2
    exit 1
fi

HOST_FASTA=$1
SYMBIONT_FASTA=$2
OUTPUT_FASTA=${3:-GCA_964034985.1_Acer_S.fitti.fa}

ml gcc/9.2.0 bwa/0.7.17

cat "$HOST_FASTA" "$SYMBIONT_FASTA" > "$OUTPUT_FASTA"
bwa index "$OUTPUT_FASTA"

