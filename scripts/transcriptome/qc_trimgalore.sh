#!/usr/bin/env bash
set -ue

R1=$1
R2=$2
OUT_DIR=$3
THREADS=4

mkdir -p ${OUT_DIR}

trim_galore \
    --paired \
    --quality 20 \
    --length 50 \
    --stringency 3 \
    --fastqc \
    --cores ${THREADS} \
    ${R1} ${R2} \
    -o ${OUT_DIR}

echo "Trim Galore filtering finished. Clean RNA-seq reads in ${OUT_DIR}"
