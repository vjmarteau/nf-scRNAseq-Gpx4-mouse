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
library(tidyverse)
library(DESeq2)
library(IHW)

conflict_prefer("filter", "dplyr")

# Load parameters
count_mat <- read_csv(arguments$count_mat)
metadata <- read_csv(arguments$metadata)
metadata <- metadata[, -1]

output_file <- arguments$output_file

# Data wrangling
count_mat <- count_mat |> column_to_rownames(var = "gene_id")
metadata <- metadata |>
  column_to_rownames(var = "sample") |>
  mutate_all(factor) |>
  mutate(group = fct_relevel(group, c("WT", "KO") ))

run_DESeq2 <- function(counts, meta) {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = meta,
                                  design = ~batch + group)
    dds$group <- relevel(dds$group, ref = "WT")
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)
    dds <- results(dds, filterFun = ihw, contrast = c("group", "KO", "WT"))
    dds <- as.data.frame(dds) |>
      rownames_to_column(var = "symbol") |>
      as_tibble() |>
      arrange(padj)
    
    return(dds)
}

DGEA <-  run_DESeq2(counts = count_mat, meta = metadata)

write_tsv(DGEA, output_file)