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

# Get marker genes
marker_genes = pd.read_csv(marker_genes)

ah = sh.annotation.AnnotationHelper(markers=marker_genes)

sc.tl.leiden(adata, key_added="leiden", resolution=0.25)

ah.annotate_cell_types(
    adata,
    { "leiden_res0_25"
        "B cells": [0],
        "Cholangiocytes": [1],
        "Endothelial cells": [2],
        "Hepatocytes": [3],
        "Plasma cells": [4],
        "Progenitor cells": [5, 6],
        "Mast cells": [7],
        "Neutrophils": [8],
        "myeloid": [9, 10],
    },
)

adata.write(f"{resDir}/annotated_adata.h5ad", compression="gzip")

def save_pseudobulk(pb, samplesheet_filename, counts_filename):
    samplesheet = pb.obs.copy()
    samplesheet.drop(columns=["sample"], inplace=True)
    samplesheet.index.name = "sample"
    samplesheet.reset_index(inplace=True)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_filename)
    bulk_df.to_csv(counts_filename)

for ct, tmp_ad in tqdm(sh.util.split_anndata(adata, "leiden_res0_25")):
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="sample",
        groups_col="group",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=3
    )
    save_pseudobulk(pb, f"{resDir}/leiden_{ct}_samplesheet.csv", f"{resDir}/leiden_{ct}_counts.csv")