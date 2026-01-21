Plastome Analysis Pipeline

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

Chloroplast-type configuration with seed sequence

Script: scripts/genome_assembly/assembly_novoplasty.sh

3. Read Mapping and IR/SC Verification

Tool: Bowtie2 v2.3.5

Script: scripts/genome_assembly/mapping_bowtie2.sh

Manual verification of inverted repeat (IR) and single-copy (SC) region boundaries

4. Genome Annotation

Initial annotation: Geneious Prime v2023.0.1

Reference: model plants and closely related plastomes (≥60% similarity)

ORFs manually curated using model plant CDS as reference:

Coverage ≥70% → predicted as intact genes

Coverage <70% → predicted as pseudogenes/missing

NDH genes special rule: if any NDH gene is pseudogenized or lost, all NDH genes are considered functionally lost

tRNA annotation: tRNAscan-SE v2.0 (organelle mode -O, intron detection -I)
Script: scripts/annotation/trna_scan.sh

ORF alignment: scripts/annotation/orf_alignment.py

5. Transcriptome Analysis

RNA-seq reads mapped to plastomes using TOPHAT2
Parameters: --library-type fr-firststrand, --read-mismatches 4, --read-gap-length 0, --read-edit-dist 4, --max-insertion-length 0, --max-deletion-length 0, --coverage-search

Read counts obtained with Geneious Prime

FPKM calculation: scripts/transcriptome/fpkm_calculation.py

Differential expression: DESeq2 v3.22 (scripts/transcriptome/dge_deseq2.R)

RNA editing sites detected in protein-coding transcripts using Geneious Prime
Thresholds: minimum coverage 60, minimum variant frequency 10%, maximum variant P-value 10^-6, strand-bias P-value threshold 10^-5 for >65% bias

CDS alignment and phylogeny (PhyloSuite v1.2.3 workflow):

Alignment: MAFFT v7.511 (codon mode)

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

CDS regions:

Alignment: MAFFT v7.511 (codon mode) → MACSE v2.06 → Gblocks v0.91b

Gblocks parameters: min block length = 5, max contiguous nonconserved = 8, no gaps

Phylogeny reconstruction: PhyloSuite v1.2.3 → IQ-TREE v2.0.7 (ModelFinder, BIC)

Branch length estimation: HyPhy v2.2.4 (MG94×GTR model)

Selection analysis: CODEML (PAML v4.10.7) with null model H0 vs alternative HA

Likelihood ratio test: LMAP v1.0.2

Substitution rate extraction and processing:

CDS: CODEML (pairwise mode, codonFreq = F3×4) → dN and dS

Non-CDS: BASEML (model = HKY85, using the phylogenetic tree) → base substitution rates

Pairwise matrices processed by Python script process_substitution.py → long-format table (groups, subgroups, value)

Only one substitution rate type is processed per run (dN, dS, or base substitution rate)

Visualization and statistical testing:

Long-format table plotted with R script plot_substitution_rate.R → boxplots, individual points, pairwise Wilcoxon tests per gene/subgroup
