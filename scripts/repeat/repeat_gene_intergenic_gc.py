#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate GC content and non-redundant length of repeat regions
located in gene and intergenic regions.

Input:
  - FASTA genomes
  - BED file of gene coordinates
  - Non-redundant repeat interval files

Output:
  - repeat_gene_intergenic_summary.txt

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


def merge_intervals(segs):
    if not segs:
        return []
    segs.sort(key=lambda x: x[0])
    merged = []
    cur_start, cur_end = segs[0]

    for s, e in segs[1:]:
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


# -------------------------------
# Main
# -------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Classify repeat regions into gene and intergenic regions and calculate GC content"
    )
    parser.add_argument("-p", "--path", required=True, help="Working directory containing fasta and repeat results")
    parser.add_argument("-e", "--evalue_dir", required=True, help="Repeat e-value directory name (e.g. 6)")
    parser.add_argument("-b", "--bed", required=True, help="Gene BED file")
    parser.add_argument("-o", "--out", default="repeat_gene_intergenic_summary.txt", help="Output summary file")

    args = parser.parse_args()

    path = args.path
    evalue_dir = os.path.join(path, args.evalue_dir)
    gene_bed_file = args.bed
    output_file = args.out

    # -------------------------------
    # Load gene intervals
    # -------------------------------
    gene_intervals = {}
    with open(gene_bed_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            seqid = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            gene_intervals.setdefault(seqid, []).append((start, end))

    # -------------------------------
    # Process each genome
    # -------------------------------
    with open(output_file, 'w') as out_f:
        out_f.write("Genome\tRegion_type\tNonredundant_length(bp)\tGC(%)\n")

        for file in os.listdir(path):
            if not file.endswith(".fasta"):
                continue

            genome = file.replace(".fasta", "")
            fasta_file = os.path.join(path, file)
            repeat_file = os.path.join(evalue_dir, genome + "_nonredundant_intervals.txt")

            if not os.path.exists(repeat_file):
                print(f"[Warning] {genome}: repeat file not found, skipped.")
                continue

            print(f"\nProcessing genome: {genome}")

            seq_dict = read_fasta(fasta_file)

            repeat_intervals = {}
            with open(repeat_file) as f:
                for line in f:
                    seqid, s, e = line.strip().split('\t')
                    repeat_intervals.setdefault(seqid, []).append((int(s), int(e)))

            gene_total_seqs = []
            intergenic_total_seqs = []

            for seqid, repeats in repeat_intervals.items():
                genes = gene_intervals.get(seqid, [])
                merged_repeats = merge_intervals(repeats)

                for r_start, r_end in merged_repeats:
                    overlap = False
                    for g_start, g_end in genes:
                        if r_end >= g_start and r_start <= g_end:
                            overlap = True
                            break

                    seq = seq_dict[seqid][r_start-1:r_end]
                    if overlap:
                        gene_total_seqs.append(seq)
                    else:
                        intergenic_total_seqs.append(seq)

            gene_len, gene_gc = calc_len_gc(gene_total_seqs)
            inter_len, inter_gc = calc_len_gc(intergenic_total_seqs)

            out_f.write(f"{genome}\tGene\t{gene_len}\t{gene_gc}\n")
            out_f.write(f"{genome}\tIntergenic\t{inter_len}\t{inter_gc}\n")

            print(f"{genome}: Gene repeats {gene_len} bp, GC={gene_gc}% | "
                  f"Intergenic repeats {inter_len} bp, GC={inter_gc}%")

    print("\nAll genomes processed successfully.")


if __name__ == "__main__":
    main()
