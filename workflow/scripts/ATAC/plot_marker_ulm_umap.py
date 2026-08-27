if "snakemake" in locals():
    scores_csv = snakemake.input["scores"]
    pvals_csv = snakemake.input["pvals"]
    coords_csv = snakemake.input["coords"]
    plot_path = snakemake.output["plot"]
    signature_label = snakemake.params.get("signature_label", "Marker signature")
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--scores-csv", default="results/ATAC/level1_marker_ulm_scores.csv")
    parser.add_argument("--pvals-csv", default="results/ATAC/level1_marker_ulm_pvals.csv")
    parser.add_argument("--coords-csv", default="data/ATAC/coordinates.csv")
    parser.add_argument("--plot-path", default="plots/ATAC/level1_marker_ulm_umap.pdf")
    parser.add_argument("--signature-label", default="Level1 marker signature")
    args = parser.parse_args()

    scores_csv = args.scores_csv
    pvals_csv = args.pvals_csv
    coords_csv = args.coords_csv
    plot_path = args.plot_path
    signature_label = args.signature_label

from math import ceil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REQUIRED_COORD_COLUMNS = {"cell", "X", "Y"}
CELL_COLUMN = "cell"
POINT_SIZE = 0.18
N_COLUMNS = 5


def validate_columns(data, required, label):
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}.")


coords = pd.read_csv(coords_csv)
coords = coords.loc[:, ~coords.columns.str.startswith("Unnamed:")]
validate_columns(coords, REQUIRED_COORD_COLUMNS, "ATAC UMAP coordinates")

scores = pd.read_csv(scores_csv)
pvals = pd.read_csv(pvals_csv)
validate_columns(scores, {CELL_COLUMN}, "ULM score table")
validate_columns(pvals, {CELL_COLUMN}, "ULM p-value table")

score_columns = [column for column in scores.columns if column != CELL_COLUMN]
pval_columns = [column for column in pvals.columns if column != CELL_COLUMN]
if not score_columns:
    raise ValueError(
        f"ULM score table does not contain any {signature_label} score columns."
    )
if score_columns != pval_columns:
    raise ValueError(
        "ULM score and p-value tables have different signature columns or ordering."
    )
if scores[CELL_COLUMN].duplicated().any():
    raise ValueError("ULM score table contains duplicated cell IDs.")
if coords[CELL_COLUMN].duplicated().any():
    raise ValueError("ATAC UMAP coordinates contain duplicated cell IDs.")

plot_data = coords[[CELL_COLUMN, "X", "Y"]].merge(scores, on=CELL_COLUMN, how="inner")
if plot_data.empty:
    raise ValueError("No overlapping cells were found between UMAP coordinates and ULM scores.")

missing_score_cells = set(scores[CELL_COLUMN]).difference(set(coords[CELL_COLUMN]))
if missing_score_cells:
    raise ValueError(
        "ATAC UMAP coordinates are missing "
        f"{len(missing_score_cells)} cells from the ULM score table."
    )

score_values = plot_data[score_columns].to_numpy(dtype=float)
finite_scores = score_values[np.isfinite(score_values)]
if finite_scores.size == 0:
    raise ValueError("ULM score table does not contain any finite score values to plot.")

vmin = float(np.min(finite_scores))
vmax = float(np.percentile(finite_scores, 99))
if vmin == vmax:
    vmin -= 1.0
    vmax += 1.0

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
cmap = plt.get_cmap("viridis")
n_panels = len(score_columns)
n_rows = ceil(n_panels / N_COLUMNS)

fig, axes = plt.subplots(
    n_rows,
    N_COLUMNS,
    figsize=(2.25 * N_COLUMNS, 2.35 * n_rows),
    squeeze=False,
    constrained_layout=True,
)

x_limits = (plot_data["X"].min(), plot_data["X"].max())
y_limits = (plot_data["Y"].min(), plot_data["Y"].max())

for index, column in enumerate(score_columns):
    ax = axes[index // N_COLUMNS][index % N_COLUMNS]
    panel_data = plot_data.sort_values(column, ascending=True)
    ax.scatter(
        panel_data["X"],
        panel_data["Y"],
        c=panel_data[column],
        cmap=cmap,
        norm=norm,
        s=POINT_SIZE,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(column, fontsize=9)
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)

for index in range(n_panels, n_rows * N_COLUMNS):
    axes[index // N_COLUMNS][index % N_COLUMNS].set_visible(False)

scalar_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
scalar_mappable.set_array([])
fig.colorbar(
    scalar_mappable,
    ax=axes,
    shrink=0.82,
    pad=0.01,
    label=f"{signature_label} ULM enrichment score",
)

if "snakemake" not in locals():
    Path(plot_path).parent.mkdir(parents=True, exist_ok=True)

fig.savefig(plot_path, dpi=300)
plt.close(fig)
