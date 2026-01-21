import os, re

def calc_indel_gene(seq_dir="results", output_file="results/0Indel_gene_unit_len_GC_seq.txt"):
    """
    Calculate nucleotide composition of Indel sequences intersected with CDS regions.

    Inputs:
    - no_*.indel.seq files in seq_dir
    - corresponding *.bed.merge files for matching regions

    Output:
    - Tab-delimited file with sequence name, CDS, length, sequence, and A/T/C/G proportions
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    indel_seqs = [f for f in os.listdir(seq_dir) if f.endswith(".indel.seq")]

    with open(output_file, 'w') as fo:
        for f in indel_seqs:
            bed_file = f"{seq_dir}/{f.replace('.indel.seq','.bed.merge')}"
            insert_bed = open(bed_file).read().replace('\n','\t') if os.path.exists(bed_file) else ""
            
            with open(f"{seq_dir}/{f}") as ft:
                name, gene = "", ""
                for line in ft:
                    if line.startswith(">"):
                        # Sequence header
                        name = line.strip(">\n")
                        # Extract CDS name from BED intersection
                        gene_match = re.findall(re.sub('[:-]', '\t', name) + '(.*)', insert_bed)
                        gene = gene_match[0].split('\t')[1] if gene_match else "NA"
                    else:
                        seq = line.strip()
                        # Calculate nucleotide proportions
                        a = seq.count('A')/len(seq)
                        t = seq.count('T')/len(seq)
                        c = seq.count('C')/len(seq)
                        g = seq.count('G')/len(seq)
                        fo.write(f"{name}\t{gene}\t{len(seq)}\t{seq}\t{a:.3f}\t{t:.3f}\t{c:.3f}\t{g:.3f}\n")

if __name__ == "__main__":
    calc_indel_gene()
