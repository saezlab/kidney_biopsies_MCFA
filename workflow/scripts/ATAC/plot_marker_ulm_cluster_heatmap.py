if "snakemake" in locals():
    cluster_mean_scores_csv = snakemake.input["cluster_mean_scores"]
    heatmap_pdf = snakemake.output["heatmap"]
    signature_label = snakemake.params.get("signature_label", "Marker signature")
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cluster-mean-scores-csv",
        default="results/ATAC/level1_marker_ulm_cluster_mean_scores.csv",
    )
    parser.add_argument(
        "--heatmap-pdf",
        default="plots/ATAC/level1_marker_ulm_cluster_mean_scores_heatmap.pdf",
    )
    parser.add_argument("--signature-label", default="Level1 marker signature")
    args = parser.parse_args()

    cluster_mean_scores_csv = args.cluster_mean_scores_csv
    heatmap_pdf = args.heatmap_pdf
    signature_label = args.signature_label

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


CLUSTER_COLUMN = "cluster"


def validate_summary(data, label):
    if CLUSTER_COLUMN not in data.columns:
        raise ValueError(f"{label} is missing required column: {CLUSTER_COLUMN}.")
    value_columns = [column for column in data.columns if column != CLUSTER_COLUMN]
    if not value_columns:
        raise ValueError(f"{label} does not contain any signature value columns.")
    if data[CLUSTER_COLUMN].duplicated().any():
        raise ValueError(f"{label} contains duplicated clusters.")
    return value_columns


def matrix_from_summary(data, value_columns):
    matrix = data.set_index(CLUSTER_COLUMN)[value_columns]
    return matrix.apply(pd.to_numeric, errors="coerce")


def row_zscore(matrix):
    row_means = matrix.mean(axis=1)
    row_sds = matrix.std(axis=1, ddof=1)
    zscores = matrix.sub(row_means, axis=0).div(row_sds.replace(0, np.nan), axis=0)
    return zscores.fillna(0.0)


def plot_heatmap(matrix, plot_path, cmap, norm, colorbar_label):
    figure_width = max(7.0, 0.55 * matrix.shape[1] + 2.4)
    figure_height = max(5.0, 0.32 * matrix.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(figure_width, figure_height), constrained_layout=True)
    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        norm=norm,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": colorbar_label},
    )
    ax.set_xlabel(signature_label)
    ax.set_ylabel("ATAC cluster")
    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)


cluster_mean_scores = pd.read_csv(cluster_mean_scores_csv)

score_columns = validate_summary(cluster_mean_scores, "Mean score summary")
score_matrix = matrix_from_summary(cluster_mean_scores, score_columns)
score_plot_matrix = row_zscore(score_matrix)

score_values = score_plot_matrix.to_numpy(dtype=float)
finite_scores = score_values[np.isfinite(score_values)]
if finite_scores.size == 0:
    raise ValueError("Mean score summary does not contain any finite values to plot.")
score_vmin = float(np.min(finite_scores))
score_vmax = float(np.max(finite_scores))
if score_vmin == score_vmax:
    score_vmin -= 1.0
    score_vmax += 1.0
score_norm = mpl.colors.Normalize(vmin=score_vmin, vmax=score_vmax)

if "snakemake" not in locals():
    Path(heatmap_pdf).parent.mkdir(parents=True, exist_ok=True)

plot_heatmap(
    score_plot_matrix,
    heatmap_pdf,
    cmap="RdBu_r",
    norm=score_norm,
    colorbar_label="Row Z-score of mean ULM enrichment score",
)
