import altair as alt
import pandas as pd
from natsort import natsorted


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable, *, data=None, data_col=None, **kwargs):
    """
    Discrete color scale for altair based on our color definitions.

    Parameters:
    -----------
    variable
        name of the color scale
    data
        Data frame used for the chart. If specified, will only show values that actually occur in the data.
    data_col
        If specified, check this column in `data` instead of `variable`

    Returns
    -------
    Altair color scale
    """
    tmp_colors = getattr(COLORS, variable)
    if data is not None:
        data_col = variable if data_col is None else data_col
        tmp_colors = {k: tmp_colors[k] for k in sorted(data[data_col].unique())}

    return alt.Scale(
        domain=list(tmp_colors.keys()),
        range=list(tmp_colors.values()),
        **kwargs,
    )


def altair_scale_mpl(scheme, **kwargs):
    """
    Use a continuous color scheme from mpl with altair
    """
    from matplotlib import cm
    from matplotlib.colors import to_hex

    return alt.Scale(
        range=[to_hex(x) for x in cm.get_cmap(scheme, 1000)(range(1000))], **kwargs
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")

class COLORS:

    group = {
        "WT": "#4daf4a",
        "KO": "#ff7f00",
    }

    sample = {
        "FG-1": "#000000",
        "FG-2": "#E69F00",
        "FG-3": "#56B4E9",
        "FG-4": "#009E73",
        "FG-5": "#F0E442",
        "FG-6": "#0072B2",
        "FG-7": "#D55E00",
        "FG-8": "#CC79A7",
    }

    batch = {
        "1": "#1b9e77",
        "2": "#d95f02",
        "3": "#7570b3",
        "4": "#0072B2",
    }

    phase = {
        "G1": "#1b9e77",
        "G2M": "#d95f02",
        "S": "#7570b3",
    }

    cell_type = dict(
        natsorted(
            {
                "B cell": "#b5bd61",
                "B cell dividing": "#f7b6d2",
                "Enterocyte prox": "#8c564b",
                "Enterocyte dist": "#17becf",
                "Enteroendocrine": "#ffbb78",
                "TA Progenitor": "#aa40fc",
                "Goblet": "#1f77b4",
                "Paneth": "#e377c2",
                "Stem Progenitor": "#555555",
                "T cells CD4": "#279e68",
                "T cells CD8": "#006d2c",
                "Tuft": "#98df8a",
                "MO DC": "#aec7e8",
            }.items()
        )
    )

    cell_type_rough = dict(
        natsorted(
            {
                "B cell": "#b5bd61",
                "Progenitor": "#f7b6d2",
                "Enterocyte": "#8c564b",
                "Goblet": "#17becf",
                "Enteroendocrine": "#ffbb78",
                "Paneth": "#aa40fc",
                "Tuft": "#1f77b4",
                "MO DC": "#e377c2",
                "T cell": "#279e68",
            }.items()
        )
    )
    leiden = {
        "0": "#4C78A8",
        "1": "#F58518",
        "2": "#E45756",
    }
    neutro_clusters = {
        "N0": "#4C78A8",
        "N1": "#F58518",
        "N2": "#E45756",
        "N3": "#72B7B2",
    }
    LT = {
        "No": "#4daf4a",
        "Yes": "#ff7f00",
    }
    ECD = {
        "No": "#7570b3",
        "Yes": "#d95f02",
    }
