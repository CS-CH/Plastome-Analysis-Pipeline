#!/usr/bin/env bash
set -ue

R1=$1
R2=$2
REF_FASTA=$3
OUT_DIR=$4
THREADS=4

mkdir -p ${OUT_DIR}

# Build bowtie2 index if not exist
if [ ! -f ${REF_FASTA}.1.bt2 ]; then
    bowtie2-build ${REF_FASTA} ${REF_FASTA}
fi

# Map reads to assembled plastome
bowtie2 -x ${REF_FASTA} -1 ${R1} -2 ${R2} -p ${THREADS} -S ${OUT_DIR}/mapped.sam

# Convert SAM to BAM and sort
samtools view -bS ${OUT_DIR}/mapped.sam | samtools sort -o ${OUT_DIR}/mapped.bam
samtools index ${OUT_DIR}/mapped.bam

echo "Bowtie2 mapping finished. BAM and index are in ${OUT_DIR}"
