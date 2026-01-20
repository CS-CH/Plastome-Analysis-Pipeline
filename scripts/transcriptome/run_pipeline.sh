#!/usr/bin/env bash
set -ue

###############################################################################
# Full strand-specific RNA-seq pipeline
# Input: raw paired-end FASTQ files in data_raw/fastq/
# Output: strand-specific FPKM tables in data_processed/fpkm/
# Requirements:
#   - TopHat2 or HISAT2
#   - samtools
#   - featureCounts (Subread)
#   - Python 3 with pandas/numpy
###############################################################################

FASTQ_DIR="data_raw/fastq"
REF_INDEX="data_raw/reference/reference_index"
GTF="data_raw/reference/reference.gtf"
BAM_DIR="data_processed/bam"
COUNT_DIR="data_processed/counts"
FPKM_DIR="data_processed/fpkm"
THREADS=4

mkdir -p "${BAM_DIR}" "${COUNT_DIR}" "${FPKM_DIR}"

# ----------------------------
# 1. Loop over paired-end FASTQ
# ----------------------------
for R1 in ${FASTQ_DIR}/*_R1.fastq
do
    SAMPLE=$(basename ${R1} _R1.fastq)
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq"
    OUT_DIR="${BAM_DIR}/${SAMPLE}"

    echo "=== Mapping sample: ${SAMPLE} ==="

    # TopHat2 mapping
    tophat2 -p ${THREADS} \
        --library-type fr-firststrand \
        -G ${GTF} \
        -o ${OUT_DIR} \
        ${REF_INDEX} \
        ${R1} ${R2}

    BAM="${OUT_DIR}/accepted_hits.bam"

    # ----------------------------
    # 2. Split BAM into forward/reverse strands
    # ----------------------------
    echo "Splitting BAM for strand-specific counts"
    bash scripts/transcriptome/split.sh ${BAM}

    FWD_BAM="fwd.bam"
    REV_BAM="rev.bam"

    # ----------------------------
    # 3. Count reads per gene with featureCounts
    # ----------------------------
    echo "Counting reads per gene (forward strand)"
    featureCounts -a ${GTF} -o ${COUNT_DIR}/counts_${SAMPLE}_fwd.txt -s 1 -p -T ${THREADS} ${FWD_BAM}

    echo "Counting reads per gene (reverse strand)"
    featureCounts -a ${GTF} -o ${COUNT_DIR}/counts_${SAMPLE}_rev.txt -s 2 -p -T ${THREADS} ${REV_BAM}

    # ----------------------------
    # 4. Calculate FPKM with Python
    # ----------------------------
    echo "Calculating FPKM"
    python scripts/transcriptome/fpkm_calculation.py \
        --fwd "${COUNT_DIR}/counts_${SAMPLE}_fwd.txt" \
        --rev "${COUNT_DIR}/counts_${SAMPLE}_rev.txt" \
        --out_fwd "${FPKM_DIR}/fpkm_${SAMPLE}_fwd.tsv" \
        --out_rev "${FPKM_DIR}/fpkm_${SAMPLE}_rev.tsv"

done

echo "=== All samples processed. FPKM tables are in ${FPKM_DIR} ==="
