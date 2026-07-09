if "snakemake" in locals():
    coords_path = snakemake.input["coords"]
    plot_path = snakemake.output["plot"]
else:
    coords_path = "results/glue/coembedding_umap.csv"
    plot_path = "plots/figures/Figure1/glue_coembedding_rna_atac_annotations.pdf"

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


UNKNOWN_LABEL = "Unannotated"
UNKNOWN_COLOR = "#D3D3D3"

RNA_COLORS = {
    "PT": "#C05A3A",
    "aPT": "#D66F45",
    "dPT": "#8F4B2F",
    "TAL": "#C2B35A",
    "DCT": "#3F8D3A",
    "CNT": "#6BD868",
    "PC": "#4FAE78",
    "IC-A": "#3D8F93",
    "IC-B": "#8AD5CC",
    "FIB": "#7C55A5",
    "VSM-P": "#5E3E88",
    "MES": "#B38AD5",
    "EC": "#5E8FBA",
    "Immune": "#C83FA8",
    "POD": "#D8908C",
    "PEC": "#E8B1AA",
    UNKNOWN_LABEL: UNKNOWN_COLOR,
}

ATAC_COLORS = {
    "PT": RNA_COLORS["PT"],
    "aPT": RNA_COLORS["aPT"],
    "dPT": RNA_COLORS["dPT"],
    "TAL": RNA_COLORS["TAL"],
    "DCT": RNA_COLORS["DCT"],
    "EC": RNA_COLORS["EC"],
    "Immune": RNA_COLORS["Immune"],
    "PEC": RNA_COLORS["PEC"],
    "Pod": RNA_COLORS["POD"],
    "PT_aPT": RNA_COLORS["aPT"],
    "CNT_PC": RNA_COLORS["CNT"],
    "IC": RNA_COLORS["IC-A"],
    "Fib_VSMC": RNA_COLORS["FIB"],
    UNKNOWN_LABEL: UNKNOWN_COLOR,
}

REQUIRED_COLUMNS = {
    "modality",
    "UMAP1",
    "UMAP2",
    "rna_annotation",
    "atac_cluster_label",
}


def clean_label(value):
    if pd.isna(value):
        return UNKNOWN_LABEL
    value = str(value).strip()
    if value == "" or value == "NA" or value.lower() == "nan":
        return UNKNOWN_LABEL
    return value


def validate_columns(data):
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(
            "Coembedding UMAP table is missing required columns: "
            f"{', '.join(sorted(missing))}."
        )


def validate_colors(data, label_column, colors, label):
    observed = set(data[label_column].map(clean_label))
    missing = sorted(observed.difference(colors))
    if missing:
        raise ValueError(
            f"Missing {label} colors for labels: {', '.join(missing)}."
        )


def ordered_labels(data, label_column, colors):
    observed = set(data[label_column])
    labels = [label for label in colors if label in observed and label != UNKNOWN_LABEL]
    if UNKNOWN_LABEL in observed:
        labels.append(UNKNOWN_LABEL)
    return labels


def plot_panel(ax, data, label_column, colors, title, point_size):
    for label in ordered_labels(data, label_column, colors):
        label_data = data.loc[data[label_column] == label]
        ax.scatter(
            label_data["UMAP1"],
            label_data["UMAP2"],
            s=point_size,
            alpha=0.65,
            c=[colors[label]],
            label=f"{label} (n={len(label_data):,})",
            linewidths=0,
            rasterized=True,
        )

    ax.set_title(title)
    ax.set_xlabel("GLUE UMAP 1")
    ax.set_ylabel("GLUE UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
        markerscale=8,
        fontsize=6,
    )


coords = pd.read_csv(coords_path, low_memory=False)
validate_columns(coords)

coords["rna_annotation"] = coords["rna_annotation"].map(clean_label)
coords["atac_cluster_label"] = coords["atac_cluster_label"].map(clean_label)

rna = coords.loc[coords["modality"] == "RNA"].copy()
atac = coords.loc[coords["modality"] == "ATAC"].copy()
if rna.empty:
    raise ValueError("No RNA cells were found in the coembedding UMAP table.")
if atac.empty:
    raise ValueError("No ATAC cells were found in the coembedding UMAP table.")

validate_colors(rna, "rna_annotation", RNA_COLORS, "RNA annotation")
validate_colors(atac, "atac_cluster_label", ATAC_COLORS, "ATAC annotation")

x_limits = (coords["UMAP1"].min(), coords["UMAP1"].max())
y_limits = (coords["UMAP2"].min(), coords["UMAP2"].max())

if "snakemake" not in locals():
    Path(plot_path).parent.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
plot_panel(
    axes[0],
    rna,
    "rna_annotation",
    RNA_COLORS,
    "RNA coembedding by annotation",
    point_size=0.6,
)
plot_panel(
    axes[1],
    atac,
    "atac_cluster_label",
    ATAC_COLORS,
    "ATAC coembedding by annotation",
    point_size=2.0,
)

for ax in axes:
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)

fig.savefig(plot_path, dpi=300)
plt.close(fig)
