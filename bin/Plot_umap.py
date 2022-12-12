#!/usr/bin/env python3

"""
Usage:
  Plot_umap.py --adata=<adata> --marker_genes=<marker_genes> [options]

Mandatory arguments:
  --adata=<adata>       adata
  --marker_genes=<marker_genes>   marker genes as csv

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import scanpy_helpers as sh

# Only needed for processing and plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

args = docopt(__doc__)
adata = args["--adata"]
marker_genes = args["--marker_genes"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)

# Get marker genes
marker_genes = pd.read_csv(marker_genes)
ah = sh.annotation.AnnotationHelper(markers=marker_genes)

# Remove genes from being plotted
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Ccl6")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Spink4")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Derl3")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Igfbp4")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Ssr4")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "C1qa")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "C1qb")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "C1qc")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Ccl20")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Pecam1")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Pla2g7")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Cd74")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Derl3")]
marker_sub = marker_genes.loc[(marker_genes.gene_identifier != "Rgcc")]

for t in marker_sub.cell_type.unique():
        col = marker_sub.loc[(marker_sub.cell_type==t), 'gene_identifier'].values
        col_w = list(set(col).intersection(adata.var.index))
        #if 'Ccl6' in col_w:
        #    col_w.remove('Ccl6')
        
        with plt.rc_context({"figure.figsize": (6, 6)}):

            fig = sc.pl.umap(adata, color=col_w, cmap="inferno", return_fig=True, size=10)
            fig.savefig(f"{resDir}/umap_{t}_marker_genes.png", bbox_inches="tight")

# Plot covariates
for col in ["cell_type", "cell_type_rough", "sample", "group", "sex", "phase"]:
    with plt.rc_context({"figure.figsize": (6, 6)}):
        
        fig = sc.pl.umap(adata, color=col, return_fig=True, size=10)
        fig.savefig(f"{resDir}/umap_covar_{col}.png", bbox_inches="tight")

# Plot leiden
with plt.rc_context({"figure.figsize": (6, 6)}):
  fig = sc.pl.umap(adata, color=["leiden"], legend_loc="on data", legend_fontoutline=2, return_fig=True, size=10)
  fig.savefig(f"{resDir}/umap_leiden.png", bbox_inches="tight")

# Plot total counts
with plt.rc_context({"figure.figsize": (6, 6)}):
  fig = sc.pl.umap(adata, color="total_counts", vmax=20000, cmap="inferno", return_fig=True, size=10)
  fig.savefig(f"{resDir}/umap_qc_total_counts.png", bbox_inches="tight")

# Plot pct counts mito
with plt.rc_context({"figure.figsize": (6, 6)}):
  fig = sc.pl.umap(adata, color="pct_counts_mito", cmap="inferno", return_fig=True, size=10)
  fig.savefig(f"{resDir}/umap_qc_pct_counts_mito.png", bbox_inches="tight")

# Plot GOIs
with plt.rc_context({"figure.figsize": (6, 6)}):
  fig = sc.pl.umap(adata, color=['Gpx4', 'Gpx4-ps2', 'Mt1', 'Mt2'], cmap="inferno", return_fig=True, size=10)
  fig.savefig(f"{resDir}/umap_GOIs.png", bbox_inches="tight")