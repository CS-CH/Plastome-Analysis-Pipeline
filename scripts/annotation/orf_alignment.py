#!/usr/bin/env python3
"""
ORF alignment and reference-based coverage calculation

This script aligns query ORFs to reference CDS sequences using
Biopython PairwiseAligner in global alignment mode and calculates
reference-based coverage.

Coverage definition:
Coverage = number of non-gap reference positions aligned to non-gap
query positions / reference length (non-gap positions).

Outputs:
1. Raw coverage matrix (species x genes)
2. Gene-normalized coverage matrix
"""

import os
import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner

# ======================
# Paths and input/output
# ======================

BASE_DIR = "CDS"
INPUT_FILE = os.path.join(BASE_DIR, "input_sequences.tsv")
OUTPUT_RAW = os.path.join(BASE_DIR, "coverage_raw.tsv")
OUTPUT_NORM = os.path.join(BASE_DIR, "coverage_normalized.tsv")

REF_COL = "REF"   # column name for reference CDS

# ======================
# Alignment parameters
# ======================

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

# ======================
# Coverage calculation
# ======================

def compute_coverage(ref_aln, qry_aln):
    """
    Calculate reference-based coverage from aligned sequences.

    Parameters
    ----------
    ref_aln : str
        Aligned reference sequence
    qry_aln : str
        Aligned query sequence

    Returns
    -------
    float
        Coverage value (0–1)
    """
    ref_len = sum(1 for r in ref_aln if r != "-")
    if ref_len == 0:
        return np.nan

    covered = sum(
        1 for r, q in zip(ref_aln, qry_aln)
        if r != "-" and q != "-"
    )

    return covered / ref_len

# ======================
# Load input table
# ======================

df = pd.read_csv(INPUT_FILE, sep="\t")
df.set_index(df.columns[0], inplace=True)  # gene names as index

genes = df.index.tolist()
species = [c for c in df.columns if c != REF_COL]

# ======================
# Main computation
# ======================

coverage = pd.DataFrame(index=species, columns=genes, dtype=float)

for gene in genes:

    ref_seq = str(df.loc[gene, REF_COL])

    if not isinstance(ref_seq, str) or ref_seq.strip() == "":
        continue

    ref_seq = ref_seq.replace("N", "").replace("n", "")

    for sp in species:

        qry_seq = df.loc[gene, sp]

        if not isinstance(qry_seq, str) or qry_seq.strip() == "":
            coverage.loc[sp, gene] = np.nan
            continue

        qry_seq = qry_seq.replace("N", "").replace("n", "")

        aln = aligner.align(ref_seq, qry_seq)[0]
        ref_aln, qry_aln = aln.sequences

        cov = compute_coverage(str(ref_aln), str(qry_aln))
        coverage.loc[sp, gene] = cov

# ======================
# Output raw coverage
# ======================

coverage.to_csv(OUTPUT_RAW, sep="\t", na_rep="")

# ======================
# Gene-wise normalization
# ======================

coverage_norm = coverage.copy()

for gene in coverage_norm.columns:
    max_val = coverage_norm[gene].max(skipna=True)
    if pd.notna(max_val) and max_val > 0:
        coverage_norm[gene] /= max_val
    else:
        coverage_norm[gene] = np.nan

coverage_norm.to_csv(OUTPUT_NORM, sep="\t", na_rep="")

print("Coverage calculation finished.")
print("Raw coverage:", OUTPUT_RAW)
print("Normalized coverage:", OUTPUT_NORM)
