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


Selection Analysis (Branch and Branch-site Models)

Selection analyses were conducted using codon-based maximum likelihood models implemented in CODEML (PAML v4.10.7) and interfaced through ETE3.

Model execution

Codon alignments were analyzed gene-by-gene using a fixed species tree in Newick format. Predefined foreground branches were specified a priori based on biological hypotheses.

Branch-site models were applied to detect episodic positive selection acting on specific sites along foreground branches, including:

bsA1 (null model, ω fixed at 1 for foreground sites)

bsA (alternative model allowing ω > 1 on foreground sites)

Scripts:

scripts/selection/ete3_paml_branchsite_batch.py
This script iterates over codon alignments, links each alignment to the species tree, marks predefined foreground branches, and runs branch-site models (bsA and bsA1) using ETE3. Model outputs are organized into gene-specific directories.

Likelihood ratio tests and summary statistics

Model fit was evaluated using likelihood ratio tests (LRTs) comparing nested models (e.g., bsA vs. bsA1). Test statistics were calculated as

2ΔlnL = 2(lnL_HA − lnL_H0)

and evaluated against a chi-square distribution with degrees of freedom determined by the difference in the number of free parameters.

Scripts:

scripts/selection/paml_branch_branchsite_LRT_summary.py
This script parses CODEML output files to extract log-likelihood values (lnL), number of parameters (np), and ω (dN/dS) estimates, performs LRTs, and summarizes results across genes.

Inference of selective regimes

For branch models, shifts in selective constraint were inferred based on comparisons between foreground and background ω values:

Relaxed constraint: significantly higher ω on the foreground branch

Strengthened constraint: significantly lower ω on the foreground branch

For branch-site models, genes with a significant bsA vs. bsA1 LRT (P < 0.05) were considered candidates for episodic positive selection.

Identification of positively selected sites

For genes showing significant evidence of branch-site positive selection, Bayes Empirical Bayes (BEB) analysis was used to identify sites under selection. Sites with posterior probability ≥ 0.95 were retained. The number and proportion of BEB sites were calculated relative to protein length.

Abnormally large or boundary ω estimates (e.g., ω = 0.001 or 999) were flagged for caution during downstream interpretation.


Substitution rate extraction and processing:

CDS: CODEML (pairwise mode, codonFreq = F3×4) → dN and dS

Non-CDS: BASEML (model = HKY85, using the phylogenetic tree) → base substitution rates

Pairwise matrices processed by Python script process_substitution.py → long-format table (groups, subgroups, value)

Only one substitution rate type is processed per run (dN, dS, or base substitution rate)

Visualization and statistical testing:

Long-format table plotted with R script plot_substitution_rate.R → boxplots, individual points, pairwise Wilcoxon tests per gene/subgroup
