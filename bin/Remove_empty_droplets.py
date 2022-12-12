#!/usr/bin/env python3

"""
Usage:
  Remove_empty_droplets.py --adata=<adata> [options]

Mandatory arguments:
  --adata=<adata>       adata

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import scvi
import scanpy as sc
import pandas as pd

from threadpoolctl import threadpool_limits
import multiprocessing

def set_all_seeds(seed=0):
    import os
    import random
    import numpy as np
    import torch
    scvi.settings.seed = seed
    os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    np.random.seed(seed)  # Numpy random
    random.seed(seed)  # Python random
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # For multiGPU

args = docopt(__doc__)
adata = args["--adata"]
resDir = args["--resDir"]

threadpool_limits(16)
set_all_seeds()

# Define function to remove barcodes from adata
def rmBarcodes(adata, rm):
    barcodes = adata.obs.index.tolist()
    rmbarcodes = rm.obs.index.tolist()
    filteredBarcodes = [c for c in barcodes if c not in rmbarcodes]
    return(adata[adata.obs.index.isin(filteredBarcodes)].copy())

adata = sc.read_h5ad(adata)

# scVI first data integration

sc.pp.highly_variable_genes(adata, batch_key='sample', flavor="seurat_v3", n_top_genes=4000)

adata_raw = adata.copy()
sc.pp.normalize_total(adata_raw)
sc.pp.log1p(adata_raw)
sc.tl.pca(adata_raw)
sc.tl.pca(adata_raw, use_highly_variable=True)
adata.raw = adata_raw

adata_scvi = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="sample")

model = scvi.model.SCVI(adata_scvi)

model.train(early_stopping=True, use_gpu=True)

adata.obsm["X_scVI"] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="leiden", resolution=0.3)


# Remove empty droplets 1

# Clean small cell clusters of single sample origin
adata_rm = adata[adata.obs["leiden"] == "6", :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "1", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 5000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "7", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 1200, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "9", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 1500, :].copy()
adata = rmBarcodes(adata, adata_rm)


adata_rm = adata[adata.obs["leiden"].isin(["0", "2", "3", "4", "5", "8", "10", "11", "12"]), :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 3000, :].copy()
adata = rmBarcodes(adata, adata_rm)

# scVI

adata_scvi = adata
adata_scvi = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="sample")

model = scvi.model.SCVI(adata_scvi)

model.train(early_stopping=True, use_gpu=True)

adata.obsm["X_scVI"] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="leiden", resolution=0.3)


# Remove empty droplets 2
adata_rm = adata[adata.obs["leiden"].isin(["5", "12", "13"]), :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 8000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "10", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 6000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"].isin(["2", "11"]), :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 5000, :].copy()
adata = rmBarcodes(adata, adata_rm)


# scVI

adata_scvi = adata
adata_scvi = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="sample")

model = scvi.model.SCVI(adata_scvi)

model.train(early_stopping=True, use_gpu=True)

adata.obsm["X_scVI"] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="leiden", resolution=0.3)


adata_rm = adata[adata.obs["leiden"].isin(["7", "13"]), :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "12", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 5000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "0", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 10000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"].isin(["9"]), :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 7000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"].isin(["2", "11"]), :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 4000, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "10", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 3000, :].copy()
adata = rmBarcodes(adata, adata_rm)

# scVI no doublets

adata_scvi = adata
adata_scvi = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="sample")

model = scvi.model.SCVI(adata_scvi)

model.train(early_stopping=True, use_gpu=True)

adata.obsm["X_scVI"] = model.get_latent_representation()

sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="leiden", resolution=0.3)

adata.write(f"{resDir}/cleaned_adata.h5ad", compression="gzip")