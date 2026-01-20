#!/usr/bin/env python3
"""
Calculate strand-specific FPKM from gene counts
FPKM = (1e6 * gene_counts) / (total_mapped_reads * gene_length_in_kb)
For strand-specific RNA-seq, counts are split by forward/reverse strand,
but denominator is total mapped reads in the sample.
"""

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Calculate strand-specific FPKM")
parser.add_argument("--fwd", required=True, help="Forward strand counts file")
parser.add_argument("--rev", required=True, help="Reverse strand counts file")
parser.add_argument("--out_fwd", required=True, help="Output FPKM forward")
parser.add_argument("--out_rev", required=True, help="Output FPKM reverse")
args = parser.parse_args()

# -------------------
# Read counts
# -------------------
df_fwd = pd.read_csv(args.fwd, sep="\t")
df_rev = pd.read_csv(args.rev, sep="\t")

# Ensure 'Length' column exists (gene length in bp)
for df, strand in zip([df_fwd, df_rev], ["forward", "reverse"]):
    if 'Length' not in df.columns:
        raise ValueError(f"Counts file for {strand} strand must include 'Length' column")

# -------------------
# Total mapped reads (sum of both strands)
# -------------------
total_mapped = df_fwd['Counts'].sum() + df_rev['Counts'].sum()
if total_mapped == 0:
    raise ValueError("Total mapped reads is zero, cannot calculate FPKM.")

# -------------------
# FPKM calculation
# -------------------
df_fwd['FPKM'] = 1e6 * df_fwd['Counts'] / (total_mapped * (df_fwd['Length']/1000))
df_rev['FPKM'] = 1e6 * df_rev['Counts'] / (total_mapped * (df_rev['Length']/1000))

# -------------------
# Output
# -------------------
os.makedirs(os.path.dirname(args.out_fwd), exist_ok=True)
df_fwd[['GeneID','FPKM']].to_csv(args.out_fwd, sep="\t", index=False)
df_rev[['GeneID','FPKM']].to_csv(args.out_rev, sep="\t", index=False)

print(f"Forward FPKM saved to {args.out_fwd}")
print(f"Reverse FPKM saved to {args.out_rev}")
