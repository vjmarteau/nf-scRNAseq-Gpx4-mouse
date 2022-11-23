#!/usr/bin/env python3

"""
Usage:
  Load_adata0.py --samplesheet=<samplesheet> [options]

Mandatory arguments:
  --samplesheet=<samplesheet>       Samplesheet from nfcore/rnaseq pipeline

Optional arguments:
  --resDir=<resDir>                 Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.io

args = docopt(__doc__)
samplesheet = args["--samplesheet"]
resDir = args["--resDir"]

# Get metadata from samplesheet
meta = pd.read_csv(samplesheet, index_col='internal_id')
meta.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# Save metadata as csv
meta.to_csv(f"{resDir}/metadata.csv")

# load cellranger .h5 feature matrices
adatas_list = []
cnvan_key_l = []
for ind, sample in zip(meta.index, meta.to_dict(orient="records")):
    path_h5ads = f'./40_nfcore_scrnaseq_v2-0-0_mm39/cellranger/sample-{ind}/outs/raw_feature_bc_matrix.h5'
    tmp_adata = sc.read_10x_h5(path_h5ads)
    
    # save gene conversion key and switch index to ensembl ids before making unique
    cnvan_key_l.append(tmp_adata.var.copy())
    tmp_adata.var = tmp_adata.var.drop(columns=['feature_types','genome'])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas_list.append(tmp_adata)

for k in cnvan_key_l[1:]:
    assert np.all(k==cnvan_key_l[0])

cnvan_key = cnvan_key_l[-1]
adata = ad.concat(adatas_list, index_unique="_")

# Use conversion key to re-assign symbols to ensembl ids
adata.var.loc[cnvan_key.index,'gene_ids'] = cnvan_key.gene_ids

adata.X = scipy.sparse.csr_matrix(adata.X)

adata.write(f"{resDir}/raw_adata.h5ad", compression="gzip")

# Seperately filter sample 6
adata.obs["sample"] = pd.Categorical(adata.obs["sample"])

sample_d = dict()
for s in adata.obs['sample'].values.unique():
    _sampli = adata[adata.obs["sample"] == s, :]
    _sampli.obs["value"] = 0
    sample_d[s] = _sampli

sc.pp.filter_cells(sample_d['FG-6'], min_counts=3000)
sc.pp.filter_cells(sample_d['FG-6'], min_genes=1100)

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata.obs['sample'].values.unique():
    adatas_new_l.append(sample_d[s])
adata = ad.concat(adatas_new_l)

# very basic cell/gene filtering for all samples jointly
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=600)
sc.pp.filter_cells(adata, max_counts=100000)
sc.pp.filter_genes(adata, min_cells=10)

# annotate the group of mitochondrial genes as 'mito'
adata.var['mito'] = adata.var_names.str.startswith('mt-')
# ribosomal genes as ribo
adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs['pct_counts_ribo'] > 5, :]

adata.write(f"{resDir}/adata.h5ad", compression="gzip")
