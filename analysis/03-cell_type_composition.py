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
import altair as alt
import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy_helpers as sh
import seaborn as sns
from nxfvars import nxfvars

sc.settings.set_figure_params(figsize=(4, 4), frameon=False)

# +
# adata = ad.read_h5ad("../results/artifacts/ANNOTATE_CELL_TYPES/annotated_adata.h5ad")
# artifact_dir = "../results/artifacts/Cell_type_composition"
# -

adata_path = nxfvars.get(
    "adata_path", "../results/artifacts/ANNOTATE_CELL_TYPES/annotated_adata.h5ad"
)
artifact_dir = nxfvars.get("artifact_dir", "../results/artifacts/")
adata = ad.read_h5ad(adata_path)

# # Get cell type fractions per Sample
#
# -> For both rough and fine cell type annotations

per_sample = (
    adata.obs.groupby(["sample", "group", "cell_type"]).size().reset_index(name="n")
)

for col in ["cell_type", "cell_type_rough"]:
    per_patient = adata.obs.groupby(["sample", col]).size().reset_index(name="n")
    ch = (
        alt.Chart(per_patient)
        .mark_bar()
        .encode(
            x=alt.X("n", stack="normalize"),
            y="sample",
            color=alt.Color(col, scale=sh.colors.altair_scale(col)),
        )
    )
    per_patient.to_csv(f"{artifact_dir}/{col}_counts_per_sample.csv")
    ch.save(f"{artifact_dir}/{col}_barchart_per_patient.svg")
    ch.display()

# # Fractions per sample average for condition - stacked bar chart

pb = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col="group",
    min_prop=0.05,
    min_cells=10,
    min_counts=1000,
    min_smpls=3,
)
sc.pp.normalize_total(pb, target_sum=1e6)
sc.pp.log1p(pb)

for col in ["cell_type", "cell_type_rough"]:
    cell_type_fractions = (
        adata.obs.groupby(["sample", "group", "batch"])
        .apply(lambda x: x[col].value_counts(normalize=True).rename_axis(col))
        .reset_index(name="frac")
    )
    tmp_df = cell_type_fractions.groupby(["group", col]).agg("mean").reset_index()
    ch = (
        alt.Chart(tmp_df)
        .encode(
            x="group",
            y="frac",
            color=alt.Color(col, scale=sh.colors.altair_scale(col)),
        )
        .mark_bar()
    )
    cell_type_fractions.frac = round(cell_type_fractions.frac, 4)
    cell_type_fractions.to_csv(
        f"{artifact_dir}/{col}_fractions_sample_average_condition.csv"
    )
    ch.save(f"{artifact_dir}/{col}_fractions_sample_average_stacked_bar_chart.svg")
    ch.display()

# # transcript counts per cell type

mean_counts = (
    adata.obs.groupby(["sample", "cell_type_rough"])
    .agg(total_counts=pd.NamedAgg("total_counts", "mean"))
    .reset_index()
)

order = (
    mean_counts.groupby("cell_type_rough")
    .agg("median")
    .sort_values("total_counts")
    .index.tolist()
)

ch = (
    alt.Chart(mean_counts)
    .mark_boxplot()
    .encode(
        y=alt.Y("cell_type_rough", sort=order[::-1]),
        x=alt.X("total_counts", title="UMI counts"),
        color=alt.Color(
            "cell_type_rough",
            scale=sh.colors.altair_scale("cell_type_rough"),
            legend=None,
        ),
    )
).properties(width=300)
ch.save(f"{artifact_dir}/transcript_counts_per_cell_type.svg")
ch.display()

pb_cell_type_coarse = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col="cell_type_rough",
    min_prop=0.05,
    min_cells=10,
    min_counts=800,
    min_smpls=3,
)
sc.pp.normalize_total(pb_cell_type_coarse, target_sum=1e6)
sc.pp.log1p(pb_cell_type_coarse)

tmp_df = pb_cell_type_coarse.obs
tmp_df["Gpx4"] = np.array(pb_cell_type_coarse[:, "Gpx4"].X[:, 0])

order = (
    tmp_df.groupby("cell_type_rough")
    .agg("median")
    .sort_values("Gpx4", ascending=False)
    .index.values
)

# +
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "darkgrey"},
    "medianprops": {"color": "darkgrey"},
    "whiskerprops": {"color": "darkgrey"},
    "capprops": {"color": "darkgrey"},
}

fig, ax = plt.subplots(1, 1, figsize=(7, 4))
sns.stripplot(
    x="cell_type_rough",
    y="Gpx4",
    hue="group",
    data=tmp_df,
    ax=ax,
    order=order,
    palette=sh.colors.COLORS.group,
    size=7,
    linewidth=1,
)
sns.boxplot(
    x="cell_type_rough",
    y="Gpx4",
    ax=ax,
    data=tmp_df,
    order=order,
    color="white",
    **PROPS,
    showfliers=False,
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
_ = plt.xticks(rotation=90)
fig.savefig(f"{artifact_dir}/Gpx4_fractions.png", bbox_inches="tight")
# -

# # QC plots per cell type

for obs, figsize in zip(["cell_type_rough", "sample"], [(6, 5), (7, 5)]):
    with plt.rc_context({"figure.figsize": figsize}):
        for var, label in zip(
            ["n_genes_by_counts", "pct_counts_mito"],
            ["n_genes", "% mitochondrial reads"],
        ):
            fig, ax = plt.subplots()
            sc.pl.violin(adata, var, groupby=obs, rotation=90, ylabel=label, ax=ax)
            fig.savefig(
                f"{artifact_dir}/qc_{var}_violin_by_{obs}.png", bbox_inches="tight"
            )
