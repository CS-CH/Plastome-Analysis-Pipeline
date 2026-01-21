import os

def extract_indels(fasta_dir="data/fasta", db_dir="data/genome_db", result_dir="results"):
    """
    Extract Indel sequences from input FASTA files.
    
    Steps:
    1. Build a no-self BLAST database for each FASTA.
    2. Run BLASTn for each FASTA against its no-self database.
    3. Merge matched regions using bedtools.
    4. Subtract matched regions to obtain unmatched Indel sequences.
    """
    os.makedirs(result_dir, exist_ok=True)
    fa_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fa')]

    for fa in fa_files:
        no_self_fa = f"{db_dir}/no_{fa}"
        db_name = no_self_fa.replace('.fa', '')

        # Build BLAST database
        os.system(f"makeblastdb -in {no_self_fa} -dbtype nucl -out {db_name}")
        print(f"Database {db_name} created.")

        # Run BLASTn
        query_file = f"{fasta_dir}/{fa}"
        out_file = f"{result_dir}/{fa}.out"
        os.system(f"blastn -query {query_file} -db {db_name} -outfmt 6 -word_size 7 > {out_file}")

        # Merge BLAST matched regions
        bed_file = f"{result_dir}/{fa}.bed.merge"
        os.system(f"awk '{{print($1\"\\t\"$7\"\\t\"$8)}}' {out_file} | sort -k1,1V -k2,2n -k3,3n | bedtools merge > {bed_file}")

        # Create genome file for bedtools subtract
        genome_file = f"{result_dir}/{fa}.genome"
        seq_len = sum(1 for l in open(query_file) if not l.startswith('>'))
        with open(genome_file, 'w') as gf:
            gf.write(f"{fa}\t1\t{seq_len}\n")

        # Subtract matched regions to get Indel sequences
        no_indel_bed = f"{result_dir}/no_{fa}.bed"
        no_indel_seq = f"{result_dir}/no_{fa}.indel.seq"
        os.system(f"bedtools subtract -a {genome_file} -b {bed_file} > {no_indel_bed}")
        os.system(f"bedtools getfasta -fi {query_file} -bed {no_indel_bed} -fo {no_indel_seq}")

        print(f"Indels for {fa} extracted to {no_indel_seq}")

if __name__ == "__main__":
    extract_indels()
