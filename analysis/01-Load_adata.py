# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python [conda env:CRCA-2022-crca-scanpy]
#     language: python
#     name: conda-env-CRCA-2022-crca-scanpy-py
# ---

# # Gpx4 deficient mice:

# ## 1. Load scRNAseq data
# Load cellranger output matrices and concatenate all samples to single adata

import anndata as ad
import scanpy as sc
import scipy as sp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sc.settings.set_figure_params(figsize=(5,5))
sns.set(font_scale=2)
fig_dir="../results/figures"

# Get metadata from samplesheet
meta = pd.read_csv('../tables/samplesheet_Grabherr.csv',
                   index_col='internal_id')
meta.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# +
# load cellranger .h5 feature matrices
path_dir='/data/projects/2021/Grabherr-scRNAseq-mouse/30_nfcore_scrnaseq_v2-0-0/cellranger'
adatas_list = []
cnvan_key_l = []
for ind, sample in zip(meta.index, meta.to_dict(orient="records")):

    path_h5ads = f'{path_dir}/sample-{ind}/outs/filtered_feature_bc_matrix.h5'
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

# + [markdown] tags=[]
# ### Alternative way to load adata
# -

meta = pd.read_csv('../tables/samplesheet_Grabherr.csv')
meta.drop(axis='columns', labels=['fastq_1', 'fastq_2'], inplace=True)

# +
p_dir="/data/projects/2021/Grabherr-scRNAseq-mouse/30_nfcore_scrnaseq_v2-0-0/cellranger"
adatas = dict()
key_save_l = []
for sample in meta.to_dict(orient="records"):
    tmp_adata = sc.read_10x_h5(
        f"{path_dir}/sample-{sample['sample']}/outs/filtered_feature_bc_matrix.h5"
    )
    # save gene conversion key and switch index to ensembl ids before making unique
    key_save_l.append(tmp_adata.var.copy())
    tmp_adata.var = tmp_adata.var.drop(columns=['feature_types','genome'])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas[sample['sample']] = tmp_adata # assign sample_id to barcodes
    
# when concatenating all, columns in .var are somehow dropped
# index_unique in .concat appends sample ids to barcodes
adata = ad.concat(adatas, index_unique="_")

for k in key_save_l[1:]:
    assert np.all(k==key_save_l[0])

key = key_save_l[-1]

# Use conversion key to re-assign symbols to ensembl ids
adata.var.loc[key.index,'gene_ids'] = key.gene_ids
# -

# ### Filter

# -> before any refined filtering only remove barcodes with less than 200 genes and genes found in less than 3 cells

# Basic filter thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# ### Append label as column

adata.obs['label'] = adata.obs['sample'].astype(str) + '\n' + adata.obs['group'].astype(str)

# Convert "internal_id" to factor - makes next step much faster!
adata.obs["sample"] = pd.Categorical(adata.obs["sample"])

# Append "sample_counts" as column
# First initialize with NAN and loop to fill values
adata.obs['sample_counts'] = np.NaN
for s in adata.obs["sample"].cat.categories:
    index_w_sel_sample = adata.obs.where(adata.obs['sample']==s).dropna(how='all').index
    adata.obs.loc[index_w_sel_sample,"sample_counts"] = adata.obs['sample'].value_counts().loc[s]
adata.obs['sample_counts'] = adata.obs['sample_counts'].astype(int)

adata.obs['label'] = adata.obs['label'].astype(str) + '\nn=' + adata.obs['sample_counts'].astype(str)

# ### Summary stats raw adata

# Dimensions of adata - barcodes X Genes
adata.shape

# Order highest to lowest
print(adata.obs['sample'].value_counts())
print('')
print(adata.obs['group'].value_counts())
print('')
print(adata.obs['sex'].value_counts())
print('')
print(adata.obs['batch'].value_counts())

sc.pl.highest_expr_genes(adata, n_top=20)

# Calculate QC metrics for all samples
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# ### Quality control - plot QC metrics

# Sample quality plots
fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata, 'total_counts', groupby='label', size=2, log=True, cut=0, ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_total_counts_raw_sub.pdf", bbox_inches="tight")

fig, ax = plt.subplots(1,1,figsize=(30,10))
sc.pl.violin(adata, 'pct_counts_mt', groupby='label', ax=ax)
fig.savefig(f"{fig_dir}/01-Violin_pct_counts_mt_raw_sub.pdf", bbox_inches="tight")

# ### Write AnnData object to disk

# Save h5ad
adata.write('../results/raw.h5ad', compression='gzip')
# !h5ls '../results/raw.h5ad'
