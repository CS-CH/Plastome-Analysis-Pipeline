Plastome Analysis Pipeline

A fully reproducible pipeline for plastome annotation, transcriptome validation, repeat and mismatch region analysis, and evolutionary rate estimation.

Analysis Checklist
| Input               | Step                         | Output                             |
| ------------------- | ---------------------------- | ---------------------------------- |
| Plastome assemblies | Annotation and ORF curation  | Annotated plastome files           |
| Plastome assemblies | tRNA annotation              | tRNA annotation tables             |
| RNA-seq reads       | Quality control and mapping  | BAM files                          |
| BAM files           | Expression quantification    | FPKM matrices                      |
| BAM files           | RNA editing detection        | RNA editing site tables            |
| Plastome assemblies | Repeat detection             | Short repeat BED files             |
| Plastome assemblies | Mismatch region detection    | Mismatch BED files                 |
| CDS sequences       | Alignment and filtering      | Codon alignments                   |
| Alignments          | Evolutionary rate estimation | dN/dS and substitution rate tables |
| Statistical tables  | Statistical analysis         | Final figures and tables           |

Workflow Overview
1. Plastome annotation and ORF curation

Plastome assemblies were annotated using Geneious Prime v2023.0.1 with the Live Annotation and Predict module (similarity threshold = 60%).

ORFs were aligned to tobacco CDS using Biopython PairwiseAligner:
Global alignment mode
Match score = 1
Mismatch score = −1
Gap open penalty = −5
Gap extend penalty = −1

Coverage criteria:
Coverage > 70%: predicted as intact genes
Coverage < 70%: predicted as pseudogenes or missing

If any ndh gene was pseudogenized or lost, all ndh genes were annotated as functionally lost.

Script: scripts/annotation/orf_alignment.py

2. tRNA annotation

tRNAscan-SE v2.0 was used with default parameters.

Script: scripts/annotation/trna_scan.sh

3. RNA-seq processing

Quality control:

trim_galore *.fastq.gz


Read mapping:

tophat2 --library-type fr-firststrand \
        --read-mismatches 4 \
        --read-gap-length 0 \
        --read-edit-dist 4 \
        --max-insertion-length 0 \
        --max-deletion-length 0 \
        --coverage-search


Expression quantification was performed using custom Python scripts.

Script: scripts/transcriptome/fpkm_calculation.py

4. Differential expression analysis

DESeq2 v3.22 was used for differential expression analysis.

Genes with −log10(padj) > 1 were considered significant.

Script: scripts/transcriptome/deseq2_analysis.R

5. RNA editing detection

RNA editing sites were detected in Geneious Prime using the following parameters:

Parameter	Value
Minimum coverage	60
Minimum variant frequency	10%
Maximum variant P-value	1e−6
Strand-bias P-value	1e−5
6. Repeat and mismatch region analysis

Short repeats were detected using BLASTN (word size = 7) and filtered using bedtools.

Genome alignment visualization was performed using Mauve.

Mismatch regions were extracted using bedtools complement and searched against the NCBI nt database using BLASTN.

Scripts are provided in scripts/repeats/.

7. Evolutionary rate estimation

CDS alignments were generated using MAFFT, MACSE, and Gblocks.

Model selection was performed using IQ-TREE ModelFinder with BIC criterion.

Branch length optimization was conducted in HyPhy.

Selection pressure and substitution rate analyses were performed using PAML (CODEML and BASEML).

Scripts are provided in scripts/evolution/.

8. Statistical analysis

All statistical analyses were conducted in R.

Methods used include:

Mann–Whitney test

One-way ANOVA with LSD post hoc test

Shapiro–Wilk test

Levene’s test

Pearson and Spearman correlation analyses

Script: scripts/statistics/statistics.R

Reproducibility

All scripts and notebooks in this repository are sufficient to reproduce all figures and tables from raw input data.
Each script is documented and parameterized to match the methods described in the manuscript.

Citation

If you use this pipeline, please cite the corresponding publication.
