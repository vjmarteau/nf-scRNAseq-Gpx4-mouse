#!/usr/bin/env python3

"""
Usage:
  Annotate_cell_types.py --adata=<adata> --marker_genes=<marker_genes> [options]

Mandatory arguments:
  --adata=<adata>                 adata
  --marker_genes=<marker_genes>   marker genes as csv

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
marker_genes = args["--marker_genes"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)


# Define function to remove barcodes from adata
def rmBarcodes(adata, rm):
    barcodes = adata.obs.index.tolist()
    rmbarcodes = rm.obs.index.tolist()
    filteredBarcodes = [c for c in barcodes if c not in rmbarcodes]
    return(adata[adata.obs.index.isin(filteredBarcodes)].copy())


sc.tl.leiden(adata, key_added="leiden", resolution=0.4)

# Get marker genes
marker_genes = pd.read_csv(marker_genes)

ah = sh.annotation.AnnotationHelper(markers=marker_genes)

ah.annotate_cell_types(
    adata,
    {
        "B cell": [2],
        "B cell dividing": [11],
        "T cell": [5],
        "Stem_Progenitor": [0, 1],
        "Enterocyte": [3, 7],
        "Enteroendocrine": [9],
        "Paneth": [4],
        "Tuft": [10],
        "Goblet": [6, 8],
        "MO DC": [12],
    },
    column="leiden",
    key_added="cell_type",
)

## NK/T compartment
adata_t = adata[adata.obs["cell_type"] == "T cell", :].copy()
ah.reprocess_adata_subset_scvi(adata_t, leiden_res=0.5)

ah.annotate_cell_types(
    adata_t,
    {
        "T cells CD8": [2, 3, 5, 6, 7, 9],
        "T cells CD4": [0, 1, 4],
        "NK cells": [],
        "potentially empty droplets": [8],
    },
)

ah.integrate_back(adata, adata_t)

# Remove identified empty droplets
adata_rm = adata[adata.obs["cell_type"].isin(["potentially empty droplets"]), :]
adata = rmBarcodes(adata, adata_rm)


## Stem cell compartment
adata_s = adata[adata.obs["cell_type"] == "Stem_Progenitor", :].copy()
ah.reprocess_adata_subset_scvi(adata_s, leiden_res=0.5)

ah.annotate_cell_types(
    adata_s,
    {
        "Stem Progenitor": [0, 2, 4, 6, 7],
        "TA Progenitor": [1, 3, 5, 8],
    },
)

ah.integrate_back(adata, adata_s)

## Enterocyte compartment
adata_t = adata[adata.obs["cell_type"] == "Enterocyte", :].copy()
sc.tl.leiden(adata_t, key_added="leiden", resolution=0.5)

ah.annotate_cell_types(
    adata_t,
    {
        "Enterocyte dist": [0, 3, 6, 7],
        "Enterocyte prox": [1, 2, 4, 5],
    },
)
ah.integrate_back(adata, adata_t)

## MO DC empty droplets
adata_rm = adata[adata.obs["cell_type"].isin(["MO DC"]), :]
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 3000, :].copy()
adata = rmBarcodes(adata, adata_rm)


ah.annotate_cell_types(
    adata,
    {
        "B cell": ["B cell", "B cell dividing"],
        "T cell": ["T cells CD4", "T cells CD8"],
        "Progenitor": ["Stem Progenitor", "TA Progenitor"],
        "Enterocyte": ["Enterocyte dist", "Enterocyte prox"],
        "Enteroendocrine": ["Enteroendocrine"],
        "Paneth": ["Paneth"],
        "Tuft": ["Tuft"],
        "Goblet": ["Goblet"],
        "MO DC": ["MO DC"],
    },
    column="cell_type",
    key_added="cell_type_rough",
)


adata.write(f"{resDir}/annotated_adata.h5ad", compression="gzip")

# Pseudobulk on raw count data not SCAR transformed data!
adata.X = adata.layers['counts']

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
      save_pseudobulk(pb, f"{resDir}/{ct}_samplesheet.csv", f"{resDir}/{ct}_counts.csv")


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
      save_pseudobulk(pb, f"{resDir}/{ct}_samplesheet.csv", f"{resDir}/{ct}_counts.csv")
