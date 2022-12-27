#!/usr/bin/env Rscript
'
Usage:
  DESeq2_DGEA.R --count_mat=<count_mat> --metadata=<metadata> [options]

Mandatory arguments:
  --count_mat=<count_mat>       Count matrix
  --metadata=<metadata>         Experiment metadata
  --output_file=<output_file>   Path to output file
' -> doc

# load required packages
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(conflicted)
library(readr)
library(dplyr)
library(tibble)
library(forcats)
library(DESeq2)
library(IHW)

conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

# Load parameters
count_mat <- read_csv(arguments$count_mat) |> column_to_rownames(var = "gene_id")
metadata <- read_csv(arguments$metadata)

output_file <- arguments$output_file

# Data wrangling
metadata <- metadata |>
  mutate_all(factor) |>
  mutate(group = fct_relevel(group, c("WT", "KO") ))

# Remove sample if paired sample missing because of too few cells in cell type pseudobulk
metadata_paired <- metadata |> group_by(batch) |> dplyr::filter(n() != 1) |> ungroup()
count_mat_paired <- count_mat[, colnames(count_mat) %in% metadata_paired$sample]

run_DESeq2 <- function(counts, meta) {

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta |> column_to_rownames(var = "sample"),
    design = ~batch + group)

  dds$group <- relevel(dds$group, ref = "WT")
  dds <- dds[rowSums(counts(dds)) >= 10, ]

# Likelihood ratio test
  dds_LRT <- DESeq(dds, test="LRT", full = ~batch + group, reduced=~batch)
  dds_LRT <- results(dds_LRT, contrast = c("group", "KO", "WT"), filterFun = ihw)
  dds_LRT <- as.data.frame(dds_LRT) |>
    rownames_to_column(var = "symbol") |>
    as_tibble() |>
    select(symbol, log2FoldChange, pvalue, padj) |>
    rename(lfc_LRT = log2FoldChange, pval_LRT = pvalue, padj_LRT = padj)

# LFC_shrink
  
  dds <- DESeq(dds)

  dds_lfc <- lfcShrink(dds, coef = "group_KO_vs_WT", type = "apeglm")
  dds_lfc <- as.data.frame(dds_lfc) |>
    rownames_to_column(var = "symbol") |>
    as_tibble() |>
    select(symbol, log2FoldChange, pvalue, padj) |>
    rename(lfc_Shrink = log2FoldChange, pval_Shrink = pvalue, padj_Shrink = padj)

# Wald test
  dds <- results(dds, filterFun = ihw, contrast = c("group", "KO", "WT"))
  dds <- as.data.frame(dds) |>
    rownames_to_column(var = "symbol") |>
    as_tibble() |>
    left_join(dds_lfc, by = "symbol") |>
    left_join(dds_LRT, by = "symbol") |>
    relocate(lfc_LRT, .after = log2FoldChange) |>
    relocate(padj_LRT, .after = padj) |>
    relocate(pval_LRT, .after = pvalue) |>
    relocate(lfc_Shrink, .after = lfc_LRT) |>
    relocate(padj_Shrink, .after = padj_LRT) |> 
    relocate(pval_Shrink, .after = pval_LRT) |>
    arrange(padj)

    return(dds)
}

DGEA <-  run_DESeq2(counts = count_mat_paired, meta = metadata_paired)

write_tsv(DGEA, output_file)