#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ------------------------------------------------------------
# Plot substitution rates or base substitution rates for plastid genes
#
# Input:
# - A long-format table (processed from PAML 2ML.dN, 2ML.dS, or 2-base substitution output)
# - Required columns: groups, subgroups, value
#   - groups: categorical grouping, e.g., heterotrophic ("B") vs autotrophic ("A")
#   - subgroups: gene or intron name
#   - value: substitution rate (dN, dS) or base substitution rate (e.g., 2base.t)
# - Only one type of rate per input file

library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyr)

# ----------------------------
# User parameters
# ----------------------------
input_file <- "substitution_rate_long.txt"
output_svg <- "substitution_rate_plot.svg"
output_tiff <- "substitution_rate_plot.tiff"

# Optional: reorder subgroups if needed
# Example: subgroups_order <- c("gene1", "gene2", ...)
# If not specified, order will follow factor levels in the input
# subgroups_order <- c("gene1", "gene2", ...)

# ----------------------------
# Load data
# ----------------------------
data <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)

# Ensure subgroups are factors for plotting
if (exists("subgroups_order")) {
  data$subgroups <- factor(data$subgroups, levels = subgroups_order)
} else {
  data$subgroups <- factor(data$subgroups)
}

# ----------------------------
# Boxplot with points and jitter
# ----------------------------
p <- ggplot(data, aes(x = subgroups, y = value, color = groups)) +
  geom_boxplot(width = 0.8, fill = "transparent", outlier.shape = NA, position = position_dodge(width = 0.9)) +
  geom_jitter(aes(color = groups), width = 0.2, height = 0, size = 0.2, position = position_dodge(width = 0.9)) +
  stat_summary(fun = median, geom = "point", size = 0, aes(group = interaction(groups, subgroups)), position = position_dodge(width = 0.9)) +
  scale_color_manual(values = c("darkgreen", "purple")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "right") +
  xlab("Genes/Subgroups") +
  ylab("Substitution rate")

# ----------------------------
# Pairwise Wilcoxon test
# ----------------------------
# Only between groups, per gene/subgroup
stat.test <- data %>%
  group_by(subgroups) %>%
  pairwise_wilcox_test(value ~ groups) %>%
  add_xy_position(x = "subgroups", fun = "mean_sd", dodge = 0.8)

# Overlay stat test
p <- p + stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, hide.ns = TRUE)

# ----------------------------
# Save plots
# ----------------------------
ggsave(output_svg, p, width = 10, height = 6, dpi = 300)
ggsave(output_tiff, p, width = 10, height = 6)

cat("Plot saved as SVG and TIFF.\n")
