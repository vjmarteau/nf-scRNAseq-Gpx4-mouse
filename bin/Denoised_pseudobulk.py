#!/usr/bin/env python3

"""
Usage:
  Denoised_pseudobulk.py --adata=<adata> [options]

Mandatory arguments:
  --adata=<adata>                 adata

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import pandas as pd
import scanpy_helpers as sh
import decoupler as dc
from tqdm import tqdm

args = docopt(__doc__)
adata = args["--adata"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)

# Pseudobulk on raw count data not SCAR transformed data!
adata.X = adata.layers['denoised']

def save_pseudobulk(pb, samplesheet_filename, counts_filename):
    samplesheet = pb.obs.copy()
    samplesheet.drop(columns=["sample"], inplace=True)
    samplesheet.index.name = "sample"
    samplesheet.reset_index(inplace=True)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_filename, index=False)
    bulk_df.to_csv(counts_filename)

for ct, tmp_ad in tqdm(sh.util.split_anndata(adata, "cell_type")):
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="sample",
        groups_col="group",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=3
    )
    if pb.obs['group'].nunique() <= 1:
      print(f"Cell type {ct} does not have enough replicates per group")
    else:
      save_pseudobulk(pb, f"{resDir}/{ct}_denoised_samplesheet.csv", f"{resDir}/{ct}_denoised_counts.csv")


adata_rough = adata[adata.obs["cell_type_rough"].isin(["Enterocyte", "Progenitor", "T cell"]), :].copy()

for ct, tmp_ad in tqdm(sh.util.split_anndata(adata_rough, "cell_type_rough")):
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="sample",
        groups_col="group",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=3
    )
    if pb.obs['group'].nunique() <= 1:
      print(f"Cell type {ct} does not have enough replicates per group")
    else:
      save_pseudobulk(pb, f"{resDir}/{ct}_denoised_samplesheet.csv", f"{resDir}/{ct}_denoised_counts.csv")
