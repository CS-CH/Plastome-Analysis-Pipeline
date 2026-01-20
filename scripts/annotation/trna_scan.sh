#!/usr/bin/env bash

###############################################################################
# tRNA annotation using tRNAscan-SE v2.0 (Organelle mode, intron detection, sorted)
#
# This script annotates plastid tRNA genes with tRNAscan-SE using:
# - Organelle mode (-O)
# - Intron detection (-I)
# - Max sensitivity (--max)
# - Multi-threading (--thread)
# - Sorted output
#
# Input:
#   data_raw/plastome/*.fasta
#
# Output:
#   data_processed/trna/*_trnascan_result.out
#
# Author: [Your Name]
###############################################################################

set -euo pipefail

# ======================
# Paths
# ======================

INPUT_DIR="data_raw/plastome"
OUTPUT_DIR="data_processed/trna"

mkdir -p "${OUTPUT_DIR}"

THREADS=4

# ======================
# Dependency check
# ======================

if ! command -v tRNAscan-SE &> /dev/null
then
    echo "Error: tRNAscan-SE not found in PATH."
    exit 1
fi

# ======================
# Run tRNAscan-SE for each plastome
# ======================

for fasta in ${INPUT_DIR}/*.fasta
do
    base=$(basename "${fasta}" .fasta)

    echo "Processing: ${base}"

    # Pipe 'o' to enable sorted output
    echo "o" | tRNAscan-SE \
        -O -I --max --thread ${THREADS} \
        "${fasta}" \
        -o "${OUTPUT_DIR}/${base}_trnascan_result.out"

done

echo "tRNA annotation finished (Organelle mode, intron detection, sorted)."
echo "Results saved in ${OUTPUT_DIR}"
