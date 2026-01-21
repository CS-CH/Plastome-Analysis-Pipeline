#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate non-redundant length and GC content of gene and intergenic regions.

Input:
  - BED file containing CDS, rRNA and tRNA coordinates
  - Plastome FASTA files named as <species>.fasta

Output:
  - gene_intergenic_gc_summary.txt

Author: Shengxin Chang
Project: Plastome-Analysis-Pipeline
"""

import os
import argparse


# -------------------------------
# Functions
# -------------------------------

def read_fasta(fasta_file):
    seq_dict = {}
    with open(fasta_file) as f:
        header = ''
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    seq_dict[header] = ''.join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            seq_dict[header] = ''.join(seq_lines)
    return seq_dict


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = []
    cur_start, cur_end = intervals[0]

    for s, e in intervals[1:]:
        if s <= cur_end + 1:
            cur_end = max(cur_end, e)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = s, e

    merged.append((cur_start, cur_end))
    return merged


def calc_len_gc(seq_list):
    if not seq_list:
        return 0, 0.0
    total_seq = ''.join(seq_list).upper()
    length = len(total_seq)
    gc = (total_seq.count('G') + total_seq.count('C')) / length * 100
    return length, round(gc, 2)


def get_intergenic_intervals(merged_genes, chrom_len):
    intergenic = []
    prev_end = 0
    for start, end in merged_genes:
        if start - 1 > prev_end:
            intergenic.append((prev_end + 1, start - 1))
        prev_end = end
    if prev_end < chrom_len:
        intergenic.append((prev_end + 1, chrom_len))
    return intergenic


# -------------------------------
# Main
# -------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Calculate gene and intergenic region GC content and length"
    )
    parser.add_argument("-b", "--bed", required=True, help="BED file with gene coordinates")
    parser.add_argument("-f", "--fasta_dir", default=".", help="Directory of FASTA files")
    parser.add_argument("-o", "--out", default="gene_intergenic_gc_summary.txt", help="Output file")

    args = parser.parse_args()

    bed_file = args.bed
    fasta_dir = args.fasta_dir
    output_file = args.out

    # -------------------------------
    # Read BED file
    # -------------------------------
    gene_intervals = {}
    with open(bed_file) as f:
        for line in f:
            if line.strip() == '':
                continue
            cols = line.strip().split('\t')
            species = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            gene_intervals.setdefault(species, []).append((start, end))

    # -------------------------------
    # Process each species
    # -------------------------------
    with open(output_file, 'w') as out_f:
        out_f.write("Species\tRegion_type\tNonredundant_length(bp)\tGC(%)\n")

        for species, intervals in gene_intervals.items():
            fasta_file = os.path.join(fasta_dir, f"{species}.fasta")
            if not os.path.exists(fasta_file):
                print(f"[Warning] FASTA not found for {species}, skipped.")
                continue

            seq_dict = read_fasta(fasta_file)

            # assume single plastome sequence
            seq_key = list(seq_dict.keys())[0]
            chrom_seq = seq_dict[seq_key]
            chrom_len = len(chrom_seq)

            # merge gene intervals
            merged_genes = merge_intervals(intervals)
            gene_seqs = [chrom_seq[s-1:e] for s, e in merged_genes]
            gene_len, gene_gc = calc_len_gc(gene_seqs)

            # intergenic regions
            intergenic_intervals = get_intergenic_intervals(merged_genes, chrom_len)
            intergenic_seqs = [chrom_seq[s-1:e] for s, e in intergenic_intervals]
            inter_len, inter_gc = calc_len_gc(intergenic_seqs)

            # write results
            out_f.write(f"{species}\tGene\t{gene_len}\t{gene_gc}\n")
            out_f.write(f"{species}\tIntergenic\t{inter_len}\t{inter_gc}\n")

    print("Analysis finished.")
    print(f"Results saved in: {output_file}")


if __name__ == "__main__":
    main()
