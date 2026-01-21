#!/usr/bin/env bash
set -ue

R1=$1
R2=$2
OUT_DIR=$3
THREADS=4

mkdir -p ${OUT_DIR}

fastp \
    -i ${R1} \
    -I ${R2} \
    -o ${OUT_DIR}/$(basename ${R1}) \
    -O ${OUT_DIR}/$(basename ${R2}) \
    -3 -5 -e 20 -l 50 \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --trim_poly_x \
    --qualified_quality_phred 20 \
    --n_base_limit 5 \
    --unqualified_percent_limit 40 \
    --thread ${THREADS} \
    --html ${OUT_DIR}/fastp_report.html \
    --json ${OUT_DIR}/fastp_report.json
