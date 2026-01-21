A fully reproducible pipeline for plastome annotation, transcriptome validation, repeat and mismatch region analysis, and evolutionary rate estimation.

1. Read Quality Control and Trimming

Genome reads: fastp v0.23.1

Parameters: -l 50, --qualified_quality_phred 20, --n_base_limit 5, --unqualified_percent_limit 40, --trim_poly_g, --trim_poly_x

Script: scripts/genome_assembly/qc_fastp.sh

RNA-seq reads: Trim Galore v0.6.8

Parameters: minimum Phred score = 20, minimum read length = 50, adapter stringency = 3

Script: scripts/transcriptome/qc_trimgalore.sh

2. Plastome Assembly

Tool: NOVOPlasty v4.3.1

Chloroplast-type configuration, seed sequence required

Script: scripts/genome_assembly/assembly_novoplasty.sh

3. Read Mapping and IR/SC Verification

Tool: Bowtie2 v2.3.5

Script: scripts/genome_assembly/mapping_bowtie2.sh

Manual verification of inverted repeat (IR) and single-copy (SC) region boundaries

4. Genome Annotation

Geneious Prime v2023.0.1 for initial annotation

Reference: model plants and closely related plastomes, similarity 60%

Manual curation of ORFs using model plant CDS as reference

Protein-coding ORFs: aligned to model plant CDS using Bio.Align.PairwiseAligner

Coverage >70% → intact genes

Coverage <70% → pseudogenes/missing

Special ndh rule: if any ndh gene pseudogenized or lost, all ndh genes considered functionally lost

tRNA annotation: tRNAscan-SE v2.0

Parameters: organelle mode (-O), intron detection (-I)

Script: scripts/annotation/trna_scan.sh

ORF alignment: scripts/annotation/orf_alignment.py

5. Transcriptome Analysis

RNA-seq reads mapped to corresponding plastomes using TOPHAT2

Parameters: --library-type fr-firststrand --read-mismatches 4 --read-gap-length 0 --read-edit-dist 4 --max-insertion-length 0 --max-deletion-length 0 --coverage-search

Read counts for designated regions obtained with Geneious Prime

FPKM calculation: scripts/transcriptome/fpkm_calculation.py

Differential gene expression: DESeq2 v3.22 (scripts/transcriptome/dge_deseq2.R)

RNA editing sites on protein-coding transcripts identified using Geneious Prime with thresholds:

minimum coverage = 60

minimum variant frequency = 10%

maximum variant P-value = 10^-6

strand-bias P-value threshold = 10^-5 for >65% bias

CDS alignment, cleaning, and phylogenetic analysis performed using PhyloSuite v1.2.3 (packaged workflow):

Alignment: MAFFT v7.511 (Codon mode)

Alignment optimization: MACSE v2.06

Cleaning: Gblocks v0.91b (min block length = 5, max contiguous nonconserved = 8, no gaps)

Phylogeny: IQ-TREE v2.0.7 (ModelFinder, BIC)

6. Indel Analysis

BLASTN searches and indel detection

Scripts:

scripts/indel/Indel_blastn_matches.sh

scripts/indel/calc_indel_gene.py

scripts/indel/extract_indels.py

7. Genome Repeat and Mismatch Region Analysis

BLASTN self-alignment of plastomes (word size = 7)

Repeat filtering and merging: Bedtools v2.30.0

Genome alignment visualization: Mauve v2.4.0

Scripts:

scripts/repeat/repeat_blastn_analysis.py

scripts/repeat/gene_intergenic_gc.py

scripts/repeat/repeat_gene_intergenic_gc.py

8. Evolutionary Rate and Selection Analysis

CDS alignment: MAFFT v7.511 (Codon mode) → MACSE v2.06 → Gblocks v0.91b

Gblocks parameters: min block length = 5, max contiguous nonconserved = 8, no gaps

Phylogeny: PhyloSuite v1.2.3 → IQ-TREE v2.0.7 (ModelFinder, BIC)

Branch length estimation: HyPhy v2.2.4 (MG94×GTR)

Selection pressure: CODEML (PAML v4.10.7)

Null model H0 vs alternative HA

Likelihood ratio test: LMAP v1.0.2

Base substitution rates:

CDS: CODEML (runmode=pairwise, codonFreq=F3×4)

Non-CDS: BASEML (tree + HKY85)
