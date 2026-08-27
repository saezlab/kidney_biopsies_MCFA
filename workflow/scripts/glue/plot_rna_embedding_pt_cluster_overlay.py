if "snakemake" in locals():
    coords_path = snakemake.input["coords"]
    plot_path = snakemake.output["plot"]
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--coords", default="results/glue/coembedding_umap.csv")
    parser.add_argument(
        "--plot",
        default="plots/glue/rna_embedding_pt_C3_pt_C6_overlay.pdf",
    )
    args = parser.parse_args()

    coords_path = args.coords
    plot_path = args.plot

from pathlib import Path

import matplotlib.lines as mlines
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

OVERLAY_COLORS = {
    "pt_C3": "#111111",
    "pt_C6": "#00A6D6",
}

REQUIRED_COLUMNS = {
    "modality",
    "UMAP1",
    "UMAP2",
    "rna_annotation",
    "atac_cluster",
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


def ordered_rna_labels(data):
    observed = set(data["rna_annotation"])
    labels = [label for label in RNA_COLORS if label in observed and label != UNKNOWN_LABEL]
    if UNKNOWN_LABEL in observed:
        labels.append(UNKNOWN_LABEL)
    return labels


def validate_rna_colors(data):
    observed = set(data["rna_annotation"].map(clean_label))
    missing = sorted(observed.difference(RNA_COLORS))
    if missing:
        raise ValueError(f"Missing RNA colors for labels: {', '.join(missing)}.")


coords = pd.read_csv(coords_path, low_memory=False)
validate_columns(coords)

coords["rna_annotation"] = coords["rna_annotation"].map(clean_label)
coords["atac_cluster"] = coords["atac_cluster"].map(clean_label)

rna = coords.loc[coords["modality"] == "RNA"].copy()
overlay = coords.loc[
    (coords["modality"] == "ATAC") & (coords["atac_cluster"].isin(OVERLAY_COLORS))
].copy()

if rna.empty:
    raise ValueError("No RNA cells were found in the coembedding UMAP table.")
if overlay.empty:
    raise ValueError("No ATAC pt_C3 or pt_C6 cells were found in the coembedding UMAP table.")

validate_rna_colors(rna)

expected_overlay_clusters = set(OVERLAY_COLORS)
observed_overlay_clusters = set(overlay["atac_cluster"])
missing_overlay_clusters = sorted(expected_overlay_clusters.difference(observed_overlay_clusters))
if missing_overlay_clusters:
    raise ValueError(
        "Missing expected ATAC overlay cluster(s): "
        f"{', '.join(missing_overlay_clusters)}."
    )

plot_coords = pd.concat(
    [rna[["UMAP1", "UMAP2"]], overlay[["UMAP1", "UMAP2"]]],
    ignore_index=True,
)
x_limits = (plot_coords["UMAP1"].min(), plot_coords["UMAP1"].max())
y_limits = (plot_coords["UMAP2"].min(), plot_coords["UMAP2"].max())

if "snakemake" not in locals():
    Path(plot_path).parent.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(7.2, 7.0), constrained_layout=True)

for label in ordered_rna_labels(rna):
    label_data = rna.loc[rna["rna_annotation"] == label]
    ax.scatter(
        label_data["UMAP1"],
        label_data["UMAP2"],
        s=0.35,
        alpha=0.32,
        c=[RNA_COLORS[label]],
        linewidths=0,
        rasterized=True,
    )

for cluster, color in OVERLAY_COLORS.items():
    cluster_data = overlay.loc[overlay["atac_cluster"] == cluster]
    ax.scatter(
        cluster_data["UMAP1"],
        cluster_data["UMAP2"],
        s=8.0,
        alpha=0.65,
        c="white",
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        cluster_data["UMAP1"],
        cluster_data["UMAP2"],
        s=10.0,
        facecolors="none",
        edgecolors=color,
        linewidths=0.45,
        alpha=0.95,
        rasterized=True,
    )

legend_handles = [
    mlines.Line2D(
        [],
        [],
        color=color,
        marker="o",
        linestyle="None",
        markersize=5,
        markerfacecolor="white",
        markeredgewidth=1.0,
        label=f"{cluster} ATAC overlay (n={len(overlay.loc[overlay['atac_cluster'] == cluster]):,})",
    )
    for cluster, color in OVERLAY_COLORS.items()
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    frameon=True,
    framealpha=0.92,
    fontsize=7,
)
ax.set_title("GLUE RNA embedding with ATAC pt_C3 and pt_C6 overlay")
ax.set_xlabel("GLUE UMAP 1")
ax.set_ylabel("GLUE UMAP 2")
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(x_limits)
ax.set_ylim(y_limits)
ax.set_aspect("equal", adjustable="box")

fig.savefig(plot_path, dpi=300)
plt.close(fig)

print(f"RNA cells plotted: {len(rna):,}")
for cluster in OVERLAY_COLORS:
    print(f"{cluster} ATAC overlay cells: {len(overlay.loc[overlay['atac_cluster'] == cluster]):,}")
print(f"Saved overlay plot to: {plot_path}")
