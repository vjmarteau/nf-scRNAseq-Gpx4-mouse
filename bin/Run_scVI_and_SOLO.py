#!/usr/bin/env python3

"""
Usage:
  Run_scVI_and_SOLO.py --adata=<adata> [options]

Mandatory arguments:
  --adata=<adata>       adata

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sc.settings.set_figure_params(figsize=(5,5))
sns.set(font_scale=2)

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

adata = sc.read_h5ad(adata)

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

#model.save(f"{resDir}/scVI_model", save_anndata=True)

# Doublet detection (SOLO)
solo_models = [scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=batch) for batch in adata.obs["sample"].unique()]

def run_solo(solo_model):
    solo_model.train()
    return solo_model.predict(soft=False)

solo_res = []
for solo_model in solo_models:
    solo_res.append(run_solo(solo_model))

adata.obs["is_doublet"] = pd.concat(solo_res)

solo_series = pd.concat(solo_res)
solo_series.index = solo_series.index.str.replace("-0$", "", regex=True)
adata.obs["is_doublet"] = solo_series

fig, ax = plt.subplots(1,1,figsize=(5,5))
sc.pl.umap(adata, color="is_doublet", ax=ax)
fig.savefig(f"{resDir}/is_doublet.png", bbox_inches="tight")

# Filter doublets and reprocess
adata_nodoublet = adata[adata.obs["is_doublet"] == "singlet", :].copy()
adata_nodoublet.obsm["X_scVI"] = model.get_latent_representation(adata_scvi[adata_nodoublet.obs_names, :])

sc.pp.neighbors(adata_nodoublet, use_rep="X_scVI")

sc.tl.umap(adata_nodoublet)

sc.tl.leiden(adata_nodoublet, key_added="leiden", resolution=0.3)

# Save adata nodoublet
adata_nodoublet.write(f"{resDir}/adata_nodoublet.h5ad", compression="gzip")


# Define function to remove barcodes from adata
def rmBarcodes(adata, rm):
    barcodes = adata.obs.index.tolist()
    rmbarcodes = rm.obs.index.tolist()
    filteredBarcodes = [c for c in barcodes if c not in rmbarcodes]
    return(adata[adata.obs.index.isin(filteredBarcodes)].copy())

adata = adata_nodoublet.copy()

adata_rm = adata[adata.obs["leiden"].isin(["3", "7"]), :]
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 3500, :].copy()
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"] == "1", :].copy()
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 7000, :].copy()
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

sc.tl.leiden(adata, key_added="leiden", resolution=0.5)

adata_rm = adata[adata.obs["leiden"].isin(["12", "14"]), :]
adata = rmBarcodes(adata, adata_rm)

adata_rm = adata[adata.obs["leiden"].isin(["7"]), :]
adata_rm = adata_rm[adata_rm.obs["total_counts"] < 5000, :].copy()
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

model.save(f"{resDir}/scVI_model", save_anndata=False)
adata.write(f"{resDir}/integrated_adata.h5ad", compression="gzip")