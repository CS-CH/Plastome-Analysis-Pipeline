#!/bin/bash
# BLAST comparison and merge matched regions
# Usage: ./0Indel_blastn_matches.sh <DB> <QUERY.fa>
set -ue

DB=$1       # BLAST database (no_self)
QUERY=$2    # Input FASTA file

# Run BLASTn to find matches
blastn -query $QUERY -db $DB -outfmt 6 -word_size 7 > "${QUERY}.out"

# Extract matched regions and merge overlapping intervals
awk '{print($1"\t"$7"\t"$8)}' "${QUERY}.out" | sort -k1,1V -k2,2n -k3,3n | bedtools merge > "${QUERY}.bed.merge"
