#!/usr/bin/env python3
"""
Repeat analysis based on BLASTN self-alignment.

This script performs BLASTN self-alignment of plastome genomes under multiple
E-value thresholds to identify short repeats. Redundant and overlapping repeat
regions are merged to obtain nonredundant repeat segments. The total nonredundant
repeat length and GC content are calculated for each genome.

Author: YourName
"""

import os
import argparse

# -------------------------------
# Argument parser
# -------------------------------
parser = argparse.ArgumentParser(description="Plastome repeat analysis using BLASTN self-alignment")
parser.add_argument("-i", "--input", required=True, help="Directory containing plastome fasta files")
parser.add_argument("-o", "--output", required=True, help="Output directory")
args = parser.parse_args()

path = os.path.abspath(args.input)
out_root = os.path.abspath(args.output)

EVALUE_LIST = ['9', '6', '3', '1', '1e-1', '1e-3', '1e-6']
EVALUE_FLOAT = [9, 6, 3, 1, 1e-1, 1e-3, 1e-6]
MAX_REPEAT_LEN = 500

os.makedirs(out_root, exist_ok=True)

# create evalue folders
for e in EVALUE_LIST:
    os.makedirs(os.path.join(out_root, e), exist_ok=True)

summary_dir = os.path.join(out_root, "summary")
os.makedirs(summary_dir, exist_ok=True)

summary_file = os.path.join(summary_dir, "repeat_nonredundant_length.txt")
summary = open(summary_file, "w")
summary.write("Genome\tEvalue\tNonredundant_repeat_length(bp)\tGC(%)\n")

# -------------------------------
# FASTA reader
# -------------------------------
def read_fasta(fasta_file):
    seq_dict = {}
    with open(fasta_file) as f:
        header = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seq_dict[header] = "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            seq_dict[header] = "".join(seq)
    return seq_dict

# -------------------------------
# Main loop
# -------------------------------
for fasta in os.listdir(path):
    if not fasta.endswith(".fasta"):
        continue

    base = fasta.replace(".fasta", "")
    fasta_path = os.path.join(path, fasta)
    print(f"Processing genome: {base}")

    seq_dict = read_fasta(fasta_path)

    # genome length (assume single sequence plastome)
    seqid = list(seq_dict.keys())[0]
    G_len = len(seq_dict[seqid])

    for evalue_str, evalue_float in zip(EVALUE_LIST, EVALUE_FLOAT):

        print(f"  E-value: {evalue_str}")
        out_dir = os.path.join(out_root, evalue_str)
        blast_out = os.path.join(out_dir, f"{base}.blastn.out")

        # BLASTN self-alignment
        cmd = (
            f"cd {path} && "
            f"makeblastdb -in {fasta} -dbtype nucl -parse_seqids -out {base} && "
            f"blastn -query {fasta} -db {base} -outfmt 6 -word_size 7 "
            f"-evalue {evalue_float} > {blast_out}"
        )
        os.system(cmd)

        # read blast results
        intervals = {}

        with open(blast_out) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue

                qid = cols[0]
                sid = cols[1]
                aln_len = int(cols[3])
                qstart = int(cols[6])
                qend = int(cols[7])

                # remove full-length self hit
                if qid == sid and qstart == 1 and qend == G_len:
                    continue

                # remove long repeats
                if aln_len > MAX_REPEAT_LEN:
                    continue

                s = min(qstart, qend)
                e = max(qstart, qend)
                intervals.setdefault(qid, []).append((s, e))

        # merge intervals and export
        interval_file = os.path.join(out_dir, f"{base}_nonredundant_intervals.txt")
        interval_fasta = os.path.join(out_dir, f"{base}_nonredundant_intervals.fasta")

        nonredundant_len = 0
        total_sequence = []

        with open(interval_file, "w") as f_txt, open(interval_fasta, "w") as f_fa:

            for seqid, segs in intervals.items():
                if not segs:
                    continue

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

                for s, e in merged:
                    f_txt.write(f"{seqid}\t{s}\t{e}\n")
                    seq = seq_dict[seqid][s-1:e]

                    f_fa.write(f">{seqid}_{s}_{e}\n")
                    for i in range(0, len(seq), 60):
                        f_fa.write(seq[i:i+60] + "\n")

                    nonredundant_len += e - s + 1
                    total_sequence.append(seq)

        if total_sequence:
            total_seq = "".join(total_sequence).upper()
            gc = round((total_seq.count("G") + total_seq.count("C")) / len(total_seq) * 100, 2)
        else:
            gc = 0.0

        summary.write(f"{base}\t{evalue_str}\t{nonredundant_len}\t{gc}\n")

    print(f"{base} finished\n")

summary.close()
print("All repeat analysis finished.")
