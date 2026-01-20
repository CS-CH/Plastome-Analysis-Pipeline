#!/usr/bin/env Rscript

# ================================
# Differential gene expression analysis using DESeq2 v3.22
# Input: count matrix (rows=genes, columns=samples), sample metadata
# Output: DESeq2 results table with log2 fold changes, p-values, and significance
# ================================

suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))

# -------------------------------
# Command line arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]   # Count matrix: row=gene, col=sample
meta_file   <- args[2]   # Sample metadata: SampleID, Condition
out_prefix  <- args[3]   # Output prefix

# -------------------------------
# Read input data
# -------------------------------
count_data <- read.delim(counts_file, row.names=1, check.names=FALSE)
col_data   <- read.delim(meta_file, row.names=1)

# Ensure that count columns match metadata samples
count_data <- count_data[, rownames(col_data)]

# -------------------------------
# Create DESeq2 dataset
# -------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = ~ Condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# -------------------------------
# Run DESeq2
# -------------------------------
dds <- DESeq(dds)

# Obtain results for specified contrast (Carnivorous vs Autotrophic)
res <- results(dds, contrast=c("Condition","Carnivorous","Autotrophic"))

# Shrink log2 fold changes for stability
resLFC <- lfcShrink(dds, coef="Condition_Carnivorous_vs_Autotrophic", type="apeglm")

# -------------------------------
# Export results
# -------------------------------
res_df <- as.data.frame(resLFC) %>%
  rownames_to_column("GeneID") %>%
  mutate(Significant = padj < 0.05)

write.table(res_df, file=paste0(out_prefix,"_DESeq2_results.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("DESeq2 differential expression analysis completed. Results saved to:", paste0(out_prefix,"_DESeq2_results.tsv"), "\n")
