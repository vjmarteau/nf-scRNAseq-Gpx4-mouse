#!/usr/bin/env python3

"""
Usage:
  PROGENy.py --adata=<adata> --progeny=<progeny> --dorothea=<dorothea> [options]

Mandatory arguments:
  --adata=<adata>       adata
  --progeny=<progeny>
  --dorothea=<dorothea>

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import decoupler as dc

# Only needed for processing and plotting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

args = docopt(__doc__)
adata = args["--adata"]
progeny = args["--progeny"]
dorothea = args["--dorothea"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)
progeny = pd.read_csv(progeny)
dorothea = pd.read_csv(dorothea)

# Get pseudo-bulk profile
padata = dc.get_pseudobulk(adata, sample_col='sample', groups_col='cell_type', layer='counts', min_prop=0.05, min_smpls=3, min_cells=10,)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e6)
sc.pp.log1p(padata)

logFCs, pvals = dc.get_contrast(
  padata,
  group_col='cell_type',
  condition_col='group',
  condition='KO',
  reference='WT',
  method='t-test'
  )

# Retrieve PROGENy model weights
#progeny = dc.get_progeny(organism="mouse", top=300)
#progeny.to_csv(f'../tables/progeny_mouse_2022-12-12.csv', index=False)


# Retrieve DoRothEA gene regulatory network
#dorothea = dc.get_dorothea(organism='mouse', levels=["A", "B"])


# Infer pathway activities with mlm
pathway_acts, pathway_pvals = dc.dense_run(dc.run_consensus, mat=logFCs, net=progeny, verbose=True)

# Infer DoRothEA pathway activities with mlm
tf_acts, tf_pvals = dc.dense_run(dc.run_mlm, mat=logFCs, net=dorothea, verbose=True)
tf_acts.fillna(0, inplace=True)


for ct in pathway_acts.index.tolist():
  fig = dc.plot_barplot(pathway_acts,ct, top=25, vertical=False, save=f"{resDir}/progeny_{ct}.png")
  fig = dc.plot_barplot(tf_acts, ct, top=25, vertical=True, save=f"{resDir}/dorothea_{ct}.png")

# Set nans to zero to be able to plot
pathway_acts.fillna(0, inplace=True)

sns.clustermap(pathway_acts, center=0, cmap='coolwarm')
plt.savefig(f"{resDir}/progeny_cell_types.png", bbox_inches="tight")



adata_rough = adata[adata.obs["cell_type_rough"].isin(["Enterocyte", "Progenitor", "T cell"]), :].copy()

# Get pseudo-bulk profile
padata = dc.get_pseudobulk(adata_rough, sample_col='sample', groups_col='cell_type', layer='counts', min_prop=0.05, min_smpls=3, min_cells=10,)

# Normalize
sc.pp.normalize_total(padata, target_sum=1e6) # 1e4? CPM
sc.pp.log1p(padata)

logFCs, pvals = dc.get_contrast(
  padata,
  group_col='cell_type_rough',
  condition_col='group',
  condition='KO',
  reference='WT',
  method='t-test'
  )

pathway_acts, pathway_pvals = dc.dense_run(dc.run_consensus, mat=logFCs, net=progeny, verbose=True)

# Infer DoRothEA pathway activities with mlm
tf_acts, tf_pvals = dc.dense_run(dc.run_mlm, mat=logFCs, net=dorothea, verbose=True)
tf_acts.fillna(0, inplace=True)


for ct in pathway_acts.index.tolist():
  fig = dc.plot_barplot(pathway_acts,ct, top=25, vertical=False, save=f"{resDir}/progeny_{ct}.png")
  fig = dc.plot_barplot(tf_acts, ct, top=25, vertical=True, save=f"{resDir}/dorothea_{ct}.png")

# Set nans to zero to be able to plot
pathway_acts.fillna(0, inplace=True)

sns.clustermap(pathway_acts, center=0, cmap='coolwarm')
plt.savefig(f"{resDir}/progeny_cell_types_rough.png", bbox_inches="tight")

