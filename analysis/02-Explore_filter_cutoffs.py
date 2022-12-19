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

# +
# Import packages
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import seaborn as sns

sc.settings.set_figure_params(figsize=(5, 5))
sns.set(font_scale=2)
fig_dir = "../results/figures"
# -

# ## Load raw data
# The barcode cutoff (kneeplot) from cellranger was not looking right for all samples. Therefore, load raw data and make custom filter thresholds per sample

# Get metadata from samplesheet
meta = pd.read_csv("../tables/samplesheet.csv", index_col="internal_id")
meta.drop(axis="columns", labels=["fastq_1", "fastq_2"], inplace=True)

# +
# load cellranger .h5 feature matrices
path_dir = (
    "/data/projects/2021/Grabherr-scRNAseq-mouse/30_nfcore_scrnaseq_v2-0-0/cellranger"
)
adatas_list = []
cnvan_key_l = []
for ind, sample in zip(meta.index, meta.to_dict(orient="records")):

    path_h5ads = f"{path_dir}/sample-{ind}/outs/raw_feature_bc_matrix.h5"
    tmp_adata = sc.read_10x_h5(path_h5ads)

    # save gene conversion key and switch index to ensembl ids before making unique
    cnvan_key_l.append(tmp_adata.var.copy())
    tmp_adata.var = tmp_adata.var.drop(columns=["feature_types", "genome"])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas_list.append(tmp_adata)

for k in cnvan_key_l[1:]:
    assert np.all(k == cnvan_key_l[0])

cnvan_key = cnvan_key_l[-1]
adata_raw = ad.concat(adatas_list, index_unique="_")

# Use conversion key to re-assign symbols to ensembl ids
adata_raw.var.loc[cnvan_key.index, "gene_ids"] = cnvan_key.gene_ids
# -

adata_raw.obs["sample"] = pd.Categorical(adata_raw.obs["sample"])
adata_raw.obs["group"] = pd.Categorical(adata_raw.obs["group"])
adata_raw.obs["sex"] = pd.Categorical(adata_raw.obs["sex"])
adata_raw.obs["batch"] = pd.Categorical(adata_raw.obs["batch"])

# Requires "adata_raw.obs['sample']" to be categorical
sample_d = dict()
for s in adata_raw.obs["sample"].values.unique():
    _sampli = adata_raw[adata_raw.obs["sample"] == s, :]
    _sampli.obs["value"] = 0
    sample_d[s] = _sampli

# ### Basic filtering + first look at per Sample metrics
# Apply most basic filter thresholds per sample. Mostly removes empty barcodes without counts or genes.

# +
for s in adata_raw.obs["sample"].values.unique():
    print(s + ":")
    print(sample_d[s].shape)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sc.pl.highest_expr_genes(sample_d[s], n_top=20, ax=ax)

    sc.pp.filter_cells(sample_d[s], min_genes=30)
    sc.pp.filter_cells(sample_d[s], min_counts=30)
    # Calculate QC metrics for all samples
    sample_d[s].var["mito"] = sample_d[s].var_names.str.startswith(
        "mt-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        sample_d[s], qc_vars=["mito"], percent_top=None, log1p=False, inplace=True
    )

    # sample_d[s] = sample_d[s][sample_d[s].obs.pct_counts_mito < 40, :]

    print("After min filter:")
    print(sample_d[s].shape)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sc.pl.highest_expr_genes(sample_d[s], n_top=20, ax=ax)

    # Plot QC metrics
    sc.pl.violin(
        sample_d[s],
        ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
        jitter=0.4,
        multi_panel=True,
    )

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata_raw.obs["sample"].values.unique():
    adatas_new_l.append(sample_d[s])
adata_sub = ad.concat(adatas_new_l)

adata = adata_raw[adata_sub.obs.index, adata_sub.var.index]


# -

# ### Remove outliers
# Within each sample look for outliers (mito, ribo, hemoglobin, and housekeeping genes). Last one should hopefully get rid of doublets.
# For mito, additional global cutoff > 25%, to get a better colour scale in the scatter plots. I guess everything above 25 is not viable ...
# Mito cutoff will be lower for downstream filtering.

# +
# Filter outliers
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * M.mad()) | (
        np.median(M) + nmads * M.mad() < M
    )
    return outlier


sample_d = dict()
for s in adata.obs["sample"].values.unique():
    print(s)
    adata_f = adata[adata.obs["sample"] == s, :].copy()
    adata_f.obs["value"] = 0

    # mitochondrial genes
    adata_f.var["mito"] = adata_f.var_names.str.startswith("mt-")
    # ribosomal genes
    adata_f.var["ribo"] = adata_f.var_names.str.startswith(("Rps", "Rpl", "Gm42418"))
    # hemoglobin genes.
    adata_f.var["hb"] = adata_f.var_names.str.contains(("Hb"))
    # housekeeping genes
    adata_f.var["house"] = adata_f.var_names.str.endswith(
        (
            "Actb",
            "Atp5f1",
            "B2m",
            "Gapdh",
            "Hprt",
            "Pgk1",
            "Rer1",
            "Rpl13a",
            "Rpl27",
            "Sdha",
            "Tbp",
            "Ubc",
        )
    )

    sc.pp.calculate_qc_metrics(
        adata_f,
        qc_vars=["mito", "ribo", "hb", "house"],
        inplace=True,
        percent_top=[20],
        log1p=True,
    )

    adata_f.obs["outlier"] = (
        is_outlier(adata_f, "log1p_total_counts", 5)
        | is_outlier(adata_f, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata_f, "pct_counts_in_top_20_genes", 5)
    )
    print("outlier")
    print(adata_f.obs.outlier.value_counts())

    adata_f.obs["mito_outlier"] = is_outlier(adata_f, "pct_counts_mito", 4) | (
        adata_f.obs["pct_counts_mito"] > 25
    )
    print("mito_outlier & pct_counts_mito > 25%")
    print(adata_f.obs.mito_outlier.value_counts())

    adata_f = adata_f[(~adata_f.obs.outlier) & (~adata_f.obs.mito_outlier)].copy()
    sample_d[s] = adata_f

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata.obs["sample"].values.unique():
    adatas_new_l.append(sample_d[s])
adata_sub = ad.concat(adatas_new_l)

adata_outlier = adata[adata_sub.obs.index, adata_sub.var.index]
# -

sample_d = dict()
for s in adata_outlier.obs["sample"].values.unique():
    _sampli = adata_outlier[adata_outlier.obs["sample"] == s, :]
    _sampli.obs["value"] = 0
    sample_d[s] = _sampli

# ### Second look at per Sample metrics after outlier filtering
# mt genes should be gone from top 20 highest_expr_genes

# +
for s in adata_outlier.obs["sample"].values.unique():
    print(s + ":")
    print(sample_d[s].shape)

    sc.pp.filter_genes(sample_d[s], min_cells=3)

    # Calculate QC metrics for all samples
    sample_d[s].var["mito"] = sample_d[s].var_names.str.startswith(
        "mt-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        sample_d[s], qc_vars=["mito"], percent_top=None, log1p=False, inplace=True
    )

    print("After min filter:")
    print(sample_d[s].shape)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    sc.pl.highest_expr_genes(sample_d[s], n_top=20, ax=ax)

    # Plot QC metrics
    sc.pl.violin(
        sample_d[s],
        ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
        jitter=0.4,
        multi_panel=True,
    )

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata_outlier.obs["sample"].values.unique():
    adatas_new_l.append(sample_d[s])
adata_sub = ad.concat(adatas_new_l)

adata_filtered = adata_outlier[adata_sub.obs.index, adata_sub.var.index]
# -

# ### Define per sample thresholds
# Use last code cell to visualize individual sample metrics and decide on cutoff. Still highly subjective cutoffs!

# +
# FG-1
dat = sample_d["FG-1"]
print("Sample: FG-1")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=300)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
# Use solo doublet detection instead!
# sc.pp.filter_cells(dat, max_counts = 40000)
# print('Number of cells after max count filter: {:d}'.format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=200)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-1"] = dat

# FG-2
dat = sample_d["FG-2"]
print("Sample: FG-2")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=150)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 10, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=100)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-2"] = dat

# FG-3
dat = sample_d["FG-3"]
print("Sample: FG-3")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=300)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=200)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-3"] = dat

# FG-4
dat = sample_d["FG-4"]
print("Sample: FG-4")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=300)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=200)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-4"] = dat

# FG-5
dat = sample_d["FG-5"]
print("Sample: FG-5")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=400)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=250)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-5"] = dat

# FG-6
dat = sample_d["FG-6"]
print("Sample: FG-6")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=300)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 10, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=250)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-6"] = dat

# FG-7
dat = sample_d["FG-7"]
print("Sample: FG-7")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=300)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=200)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
print(" ")
sample_d["FG-7"] = dat

# FG-8
dat = sample_d["FG-8"]
print("Sample: FG-8")
# Filter cells according to identified QC thresholds:
print("Total number of cells: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_counts=250)
print("Number of cells after min count filter: {:d}".format(dat.n_obs))
dat = dat[dat.obs.pct_counts_mito < 5, :]
print("Number of cells after MT filter: {:d}".format(dat.n_obs))
sc.pp.filter_cells(dat, min_genes=200)
print("Number of cells after gene filter: {:d}".format(dat.n_obs))
sample_d["FG-8"] = dat

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata_filtered.obs["sample"].values.unique():
    adatas_new_l.append(sample_d[s])
adata_sub = ad.concat(adatas_new_l)

adata_filtered_per_sample = adata_filtered[adata_sub.obs.index, adata_sub.var.index]
# -

# ### Plot violin per Sample quality

# Calculate QC metrics for all samples
adata_filtered_per_sample.var[
    "mito"
] = adata_filtered_per_sample.var_names.str.startswith(
    "mt-"
)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(
    adata_filtered_per_sample,
    qc_vars=["mito"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

# +
# Sample quality plots
fig, ax = plt.subplots(1, 1, figsize=(30, 10))
sc.pl.violin(
    adata_filtered_per_sample,
    "total_counts",
    groupby="sample",
    size=2,
    log=True,
    cut=0,
    ax=ax,
)

fig, ax = plt.subplots(1, 1, figsize=(30, 10))
sc.pl.violin(adata_filtered_per_sample, "pct_counts_mito", groupby="sample", ax=ax)
# -

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sc.pl.highest_expr_genes(adata_filtered_per_sample, n_top=20, ax=ax)

# Usefull plots to decide on thresholds:

# +
# Decide on Filter thresholds
dat = adata_filtered_per_sample

p1 = sc.pl.scatter(dat, "total_counts", "n_genes_by_counts", color="pct_counts_mito")
p2 = sc.pl.scatter(
    dat[dat.obs["total_counts"] < 5000],
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_mito",
)

sc.pl.violin(
    dat,
    ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
    jitter=0.4,
    multi_panel=True,
)

# Thresholding decision: counts
p3 = sns.displot(dat.obs["total_counts"], kde=False, bins=60)
plt.show()

p4 = sns.displot(
    dat.obs["total_counts"][dat.obs["total_counts"] < 1000], kde=False, bins=100
)
plt.show()

p5 = sns.displot(
    dat.obs["total_counts"][dat.obs["total_counts"] > 1000], kde=False, bins=100
)
plt.show()

p6 = sns.displot(dat.obs["pct_counts_mito"], kde=False, bins=60)
plt.show()

p7 = sns.displot(
    dat.obs["pct_counts_mito"][dat.obs["pct_counts_mito"] < 20], kde=False, bins=100
)
plt.show()

# Thresholding decision: genes
p8 = sns.displot(dat.obs["n_genes_by_counts"], kde=False, bins=60)
plt.show()

p9 = sns.displot(
    dat.obs["n_genes_by_counts"][dat.obs["n_genes_by_counts"] < 500], kde=False, bins=60
)
plt.show()
