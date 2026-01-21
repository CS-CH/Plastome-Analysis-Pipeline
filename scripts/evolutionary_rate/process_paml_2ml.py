#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evolutionary Rate Matrix Processing for PAML Output

This script reads PAML 2ML.dN or 2ML.dS files and generates:
1. Standardized square matrices with species labels
2. Pairwise table for downstream analysis or plotting

Author: Your Name
Date: 2026-01-21
"""

import numpy as np

# ---------------------------
# User Settings
# ---------------------------

# Ordered list of species for matrix
species_order = [
    'Nicotiana', 'Pouteria', 'Hydrocera', 'Marcgravia', 'Barringtonia',
    'Fouquieria', 'Saurauia', 'Polemonium', 'Diospyros', 'Aegiceras',
    'Sladenia', 'Euryodendron', 'Apterosperma', 'Symplocos', 'Pterostyrax',
    'Clethra', 'Roridula', 'Heliamphora', 'Darlingtonia', 'Sarracenia'
]

# Mapping species to numeric index (optional for matrix reordering)
species_index = {name: str(i+1) for i, name in enumerate(species_order)}

# Input file (PAML 2ML.dN or 2ML.dS)
input_file = "substitution_rate/2ML.dS"

# Output files
output_matrix_file = "substitution_rate/CDS_dS_standard.txt"
output_pairwise_file = "substitution_rate/CDS_dS_pairwise.txt"

# ---------------------------
# Read PAML Output
# ---------------------------

with open(input_file, "r") as f:
    lines = [line.strip() for line in f if line.strip() != ""]

# Dictionary to store gene: matrix
gene_dict = {}

i = 0
while i < len(lines):
    gene_name = lines[i]
    length = int(lines[i+1])
    matrix_data = lines[i+2:i+2+length]
    gene_dict[gene_name] = matrix_data
    i += 2 + length

# ---------------------------
# Convert to numeric matrix and reorder
# ---------------------------

for gene, rows in gene_dict.items():
    matrix = []
    # First row: add zero padding for symmetry
    matrix.append([0]*(len(rows)+1))
    
    for row in rows:
        row_items = row.split()
        sp_name = row_items[0]
        # Convert species name to index if in species_index
        if sp_name in species_index:
            row_items[0] = species_index[sp_name]
        # Convert to float and pad with zeros if needed
        values = [float(x) for x in row_items if x != '']
        values.extend([0]*(len(rows)-len(values)+1))
        matrix.append(values)
    
    # Make symmetric: X + X^T - diag
    X = np.array(matrix)
    X = X + X.T - np.diag(np.diag(X))
    # Format to 4 decimals
    X = [[f"{x:.4f}" for x in row] for row in X]

    # Reorder according to species_order
    label_dict = {int(row[0]): idx for idx, row in enumerate(X)}
    new_order = [label_dict[int(species_index[sp])] for sp in species_order]

    transformed = [[X[i][j] for j in new_order] for i in new_order]

    # Replace first column with species names
    transformed[0] = ['0'] + species_order
    for idx, row in enumerate(transformed[1:], start=1):
        row[0] = species_order[idx-1]
    
    gene_dict[gene] = transformed

# ---------------------------
# Write Standardized Matrix
# ---------------------------

with open(output_matrix_file, "w") as f:
    for gene, matrix in gene_dict.items():
        f.write(f"\n{gene}\n")
        for row in matrix:
            f.write(" ".join(str(x) for x in row) + "\n")

# ---------------------------
# Write Pairwise Table
# ---------------------------

with open(output_pairwise_file, "w") as f:
    for gene, matrix in gene_dict.items():
        data_rows = matrix[1:]  # skip header
        for i in range(len(data_rows)):
            for j in range(len(data_rows[i])):
                if i != j:
                    f.write(f"{gene:<20}{data_rows[i][0]:<20}{data_rows[j+1][0]:<20}{data_rows[i][j+1]}\n")
