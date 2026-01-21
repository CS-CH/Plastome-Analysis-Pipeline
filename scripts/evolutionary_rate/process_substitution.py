#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
process_substitution.py

A script to convert PAML CODEML pairwise output (2ML.dN, 2ML.dS, or 2base.t) 
into a long-format table suitable for plotting and statistical analysis.

Important:
    - Only one type of substitution rate can be processed at a time: dN, dS, or d.
    - Do not combine different rate types in a single input file.

Usage:
    python process_substitution.py <input_file> <output_file> <group_labels_file>

Arguments:
    input_file          PAML pairwise matrix (2ML.dN, 2ML.dS, or 2base.t)
    output_file         Output long-format table (tab-delimited)
    group_labels_file   Text file with two columns: species_name \t group_label (A/B etc.)

Output:
    Tab-delimited table with columns: groups, subgroups, value
    Example:
        groups    subgroups           value
        A         gene1               0.0123
        B         gene1               0.0245
        ...
"""

import numpy as np
import argparse

# ----------------------------
# Argument parsing
# ----------------------------
parser = argparse.ArgumentParser(description="Process PAML pairwise substitution matrix (dN, dS, or d)")
parser.add_argument("input_file", help="Path to PAML pairwise matrix (2ML.dN, 2ML.dS, or 2base.t)")
parser.add_argument("output_file", help="Path to output long-format table")
parser.add_argument("group_labels_file", help="Path to group labels file (species_name \\t group)")
args = parser.parse_args()

# ----------------------------
# Load group labels
# ----------------------------
group_dict = {}
with open(args.group_labels_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            species, group = parts
            group_dict[species] = group

# ----------------------------
# Read pairwise matrix
# ----------------------------
with open(args.input_file) as f:
    lines = [l.strip() for l in f if l.strip()]

species_order = lines[0].split()
matrix_lines = lines[1:]

# Convert to numeric matrix
matrix = []
for line in matrix_lines:
    values = [float(x) for x in line.split()]
    matrix.append(values)
matrix = np.array(matrix)

# ----------------------------
# Convert to long format
# ----------------------------
long_table = []

n = len(species_order)
for i in range(n):
    for j in range(n):
        if i != j:
            species_i = species_order[i]
            value = matrix[i, j]
            group = group_dict.get(species_i, "NA")
            long_table.append([group, species_i, value])

# ----------------------------
# Save output
# ----------------------------
with open(args.output_file, "w") as out:
    out.write("groups\tsubgroups\tvalue\n")
    for row in long_table:
        out.write(f"{row[0]}\t{row[1]}\t{row[2]:.4f}\n")

print(f"Processed {args.input_file} → {args.output_file}")
print("Reminder: Only one substitution rate type (dN, dS, or d) processed per run.")
