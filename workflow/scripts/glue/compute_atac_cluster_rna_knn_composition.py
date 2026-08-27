if "snakemake" in locals():
    coembedding_csv = snakemake.input["coembedding"]
    annotations_csv = snakemake.input["annotations"]
    composition_csv = snakemake.output["composition"]
    enrichment_csv = snakemake.output["enrichment"]
    summary_csv = snakemake.output["summary"]
    heatmaps_pdf = snakemake.output["heatmaps"]
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--coembedding", default="results/glue/coembedding.csv")
    parser.add_argument("--annotations", default="results/glue/coembedding_umap.csv")
    parser.add_argument(
        "--composition",
        default="results/glue/atac_cluster_rna_level2_knn_composition.csv",
    )
    parser.add_argument(
        "--enrichment",
        default="results/glue/atac_cluster_rna_level2_knn_enrichment.csv",
    )
    parser.add_argument(
        "--summary",
        default="results/glue/atac_cluster_rna_level2_knn_summary.csv",
    )
    parser.add_argument(
        "--heatmaps",
        default="plots/glue/atac_cluster_rna_level2_knn_heatmaps.pdf",
    )
    args = parser.parse_args()

    coembedding_csv = args.coembedding
    annotations_csv = args.annotations
    composition_csv = args.composition
    enrichment_csv = args.enrichment
    summary_csv = args.summary
    heatmaps_pdf = args.heatmaps

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import NearestNeighbors


K_VALUES = [25, 50, 100, 200]
PRIMARY_K = 100
QUERY_CHUNK_SIZE = 5000
UNKNOWN_LABEL = "Unannotated"
RNA_LEVEL2_ORDER = [
    "PT",
    "aPT",
    "dPT",
    "TAL",
    "DCT",
    "CNT",
    "PC",
    "IC-A",
    "IC-B",
    "FIB",
    "VSM-P",
    "MES",
    "EC",
    "Immune",
    "POD",
    "PEC",
]
ATAC_CLUSTER_ORDER = [
    "atlas_C1",
    "atlas_C12",
    "atlas_C14",
    "atlas_C3",
    "atlas_C4",
    "pc_cnt_dct_C1",
    "pc_cnt_dct_C2",
    "pt_C1",
    "pt_C2",
    "pt_C3",
    "pt_C4",
    "pt_C5",
    "pt_C6",
    "pt_C8",
    "vsm_p_fib_C1",
    "vsm_p_fib_C2",
    "vsm_p_fib_C4",
    "vsm_p_fib_C5",
]
REQUIRED_ANNOTATION_COLUMNS = {
    "coembedding_id",
    "modality",
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


def ordered_observed(values, preferred_order):
    observed = pd.Index(pd.Series(values).dropna().astype(str).unique())
    ordered = [value for value in preferred_order if value in set(observed)]
    ordered.extend(sorted(set(observed).difference(ordered)))
    return ordered


def read_inputs():
    embedding = pd.read_csv(coembedding_csv, index_col=0)
    embedding.index = embedding.index.astype(str)
    embedding.index.name = "coembedding_id"
    embedding = embedding.apply(pd.to_numeric, errors="raise")

    annotations = pd.read_csv(annotations_csv, low_memory=False)
    missing = REQUIRED_ANNOTATION_COLUMNS.difference(annotations.columns)
    if missing:
        raise ValueError(
            "Annotation table is missing required columns: "
            f"{', '.join(sorted(missing))}."
        )
    if annotations["coembedding_id"].duplicated().any():
        raise ValueError("Annotation table contains duplicated coembedding_id values.")

    annotations = annotations.set_index("coembedding_id")
    annotations.index = annotations.index.astype(str)
    annotations["rna_annotation"] = annotations["rna_annotation"].map(clean_label)
    annotations["atac_cluster"] = annotations["atac_cluster"].map(clean_label)

    joined = embedding.join(
        annotations[["modality", "rna_annotation", "atac_cluster"]],
        how="inner",
    )
    missing_annotation_count = len(embedding) - len(joined)
    if missing_annotation_count:
        print(f"Embedding rows missing annotations: {missing_annotation_count:,}")
    if joined.empty:
        raise ValueError("No embedding rows matched the annotation table.")
    return joined, embedding.shape[1]


def build_knn_counts(rna_matrix, rna_label_codes, atac_matrix, atac_clusters, n_labels):
    max_k = max(K_VALUES)
    model = NearestNeighbors(
        n_neighbors=max_k,
        algorithm="ball_tree",
        metric="euclidean",
        n_jobs=-1,
    )
    model.fit(rna_matrix)

    cluster_levels = ordered_observed(atac_clusters, ATAC_CLUSTER_ORDER)
    cluster_to_position = {cluster: i for i, cluster in enumerate(cluster_levels)}
    counts = {
        k: np.zeros((len(cluster_levels), n_labels), dtype=np.int64)
        for k in K_VALUES
    }
    cluster_cell_counts = np.zeros(len(cluster_levels), dtype=np.int64)

    for start in range(0, atac_matrix.shape[0], QUERY_CHUNK_SIZE):
        stop = min(start + QUERY_CHUNK_SIZE, atac_matrix.shape[0])
        chunk_neighbors = model.kneighbors(
            atac_matrix[start:stop],
            return_distance=False,
        )
        chunk_labels = rna_label_codes[chunk_neighbors]
        chunk_clusters = atac_clusters[start:stop]

        for cluster in np.unique(chunk_clusters):
            cluster_position = cluster_to_position[cluster]
            query_mask = chunk_clusters == cluster
            cluster_cell_counts[cluster_position] += int(query_mask.sum())
            label_matrix = chunk_labels[query_mask]
            for k in K_VALUES:
                counts[k][cluster_position, :] += np.bincount(
                    label_matrix[:, :k].reshape(-1),
                    minlength=n_labels,
                )

        print(f"Processed ATAC query cells {stop:,}/{atac_matrix.shape[0]:,}")

    return cluster_levels, cluster_cell_counts, counts


def make_tables(cluster_levels, cluster_cell_counts, counts, rna_label_levels, rna_background_counts):
    background_total = int(rna_background_counts.sum())
    background_fraction = rna_background_counts / background_total
    composition_rows = []
    enrichment_rows = []

    for k in K_VALUES:
        for cluster_index, cluster in enumerate(cluster_levels):
            n_atac_cells = int(cluster_cell_counts[cluster_index])
            denominator = n_atac_cells * k
            for label_index, label in enumerate(rna_label_levels):
                neighbor_count = int(counts[k][cluster_index, label_index])
                neighbor_fraction = neighbor_count / denominator if denominator else np.nan
                bg_fraction = float(background_fraction[label_index])
                enrichment = (
                    neighbor_fraction / bg_fraction
                    if bg_fraction > 0 and np.isfinite(neighbor_fraction)
                    else np.nan
                )
                log2_enrichment = np.log2(enrichment) if enrichment > 0 else np.nan

                composition_rows.append(
                    {
                        "k": k,
                        "atac_cluster": cluster,
                        "rna_level2": label,
                        "n_atac_cells": n_atac_cells,
                        "neighbor_count": neighbor_count,
                        "neighbor_fraction": neighbor_fraction,
                    }
                )
                enrichment_rows.append(
                    {
                        "k": k,
                        "atac_cluster": cluster,
                        "rna_level2": label,
                        "n_atac_cells": n_atac_cells,
                        "neighbor_count": neighbor_count,
                        "neighbor_fraction": neighbor_fraction,
                        "rna_background_count": int(rna_background_counts[label_index]),
                        "rna_background_fraction": bg_fraction,
                        "enrichment": enrichment,
                        "log2_enrichment": log2_enrichment,
                    }
                )

    composition = pd.DataFrame(composition_rows)
    enrichment = pd.DataFrame(enrichment_rows)

    summary_rows = []
    for (k, cluster), group in enrichment.groupby(["k", "atac_cluster"], sort=False):
        top_fraction = group.sort_values(
            ["neighbor_fraction", "log2_enrichment"],
            ascending=[False, False],
        ).iloc[0]
        top_enrichment = group.sort_values(
            ["log2_enrichment", "neighbor_fraction"],
            ascending=[False, False],
        ).iloc[0]
        fractions = group["neighbor_fraction"].to_numpy(dtype=float)
        positive_fractions = fractions[fractions > 0]
        entropy = -np.sum(positive_fractions * np.log2(positive_fractions))
        summary_rows.append(
            {
                "k": k,
                "atac_cluster": cluster,
                "n_atac_cells": int(top_fraction["n_atac_cells"]),
                "top_fraction_rna_level2": top_fraction["rna_level2"],
                "top_neighbor_fraction": top_fraction["neighbor_fraction"],
                "top_fraction_log2_enrichment": top_fraction["log2_enrichment"],
                "top_enriched_rna_level2": top_enrichment["rna_level2"],
                "top_log2_enrichment": top_enrichment["log2_enrichment"],
                "top_enriched_neighbor_fraction": top_enrichment["neighbor_fraction"],
                "neighbor_fraction_entropy": entropy,
            }
        )

    summary = pd.DataFrame(summary_rows)
    return composition, enrichment, summary


def pivot_metric(enrichment, k, metric, cluster_levels, rna_label_levels):
    data = enrichment.loc[enrichment["k"] == k]
    matrix = data.pivot(index="atac_cluster", columns="rna_level2", values=metric)
    return matrix.reindex(index=cluster_levels, columns=rna_label_levels)


def plot_heatmap_page(pdf, matrix, title, cmap, center=None, vmin=None, vmax=None, cbar_label=None):
    fig_width = max(9.0, 0.45 * matrix.shape[1] + 3.0)
    fig_height = max(7.0, 0.32 * matrix.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    sns.heatmap(
        matrix,
        ax=ax,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": cbar_label} if cbar_label else None,
    )
    ax.set_title(title)
    ax.set_xlabel("snRNA-seq level2 label")
    ax.set_ylabel("ATAC cluster")
    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)
    pdf.savefig(fig)
    plt.close(fig)


def write_heatmaps(enrichment, cluster_levels, rna_label_levels):
    sns.set_theme(style="white", context="notebook")
    with PdfPages(heatmaps_pdf) as pdf:
        fraction_primary = pivot_metric(
            enrichment,
            PRIMARY_K,
            "neighbor_fraction",
            cluster_levels,
            rna_label_levels,
        )
        log2_primary = pivot_metric(
            enrichment,
            PRIMARY_K,
            "log2_enrichment",
            cluster_levels,
            rna_label_levels,
        )
        finite_log2 = log2_primary.to_numpy(dtype=float)
        finite_log2 = finite_log2[np.isfinite(finite_log2)]
        max_abs_log2 = float(np.max(np.abs(finite_log2))) if finite_log2.size else 1.0
        max_abs_log2 = max(max_abs_log2, 1.0)

        plot_heatmap_page(
            pdf,
            fraction_primary,
            f"RNA level2 kNN fractions around ATAC clusters (k={PRIMARY_K})",
            cmap="viridis",
            vmin=0,
            vmax=max(0.25, float(np.nanmax(fraction_primary.to_numpy(dtype=float)))),
            cbar_label="Neighbor fraction",
        )
        plot_heatmap_page(
            pdf,
            log2_primary,
            f"RNA level2 kNN log2 enrichment around ATAC clusters (k={PRIMARY_K})",
            cmap="RdBu_r",
            center=0,
            vmin=-max_abs_log2,
            vmax=max_abs_log2,
            cbar_label="log2 enrichment vs RNA background",
        )

        for k in [value for value in K_VALUES if value != PRIMARY_K]:
            fraction_matrix = pivot_metric(
                enrichment,
                k,
                "neighbor_fraction",
                cluster_levels,
                rna_label_levels,
            )
            plot_heatmap_page(
                pdf,
                fraction_matrix,
                f"Sensitivity: RNA level2 kNN fractions around ATAC clusters (k={k})",
                cmap="viridis",
                vmin=0,
                vmax=max(0.25, float(np.nanmax(fraction_matrix.to_numpy(dtype=float)))),
                cbar_label="Neighbor fraction",
            )


joined, n_dims = read_inputs()
dimension_columns = [column for column in joined.columns if column not in {"modality", "rna_annotation", "atac_cluster"}]
if len(dimension_columns) != n_dims:
    raise ValueError("Unexpected joined coordinate column count.")

rna = joined.loc[joined["modality"] == "RNA"].copy()
atac = joined.loc[
    (joined["modality"] == "ATAC") & (joined["atac_cluster"] != UNKNOWN_LABEL)
].copy()

if rna.empty:
    raise ValueError("No RNA cells found in joined GLUE embedding.")
if atac.empty:
    raise ValueError("No named ATAC cells found in joined GLUE embedding.")

rna_label_levels = ordered_observed(rna["rna_annotation"], RNA_LEVEL2_ORDER)
rna_label_to_code = {label: i for i, label in enumerate(rna_label_levels)}
rna_label_codes = rna["rna_annotation"].map(rna_label_to_code).to_numpy(dtype=np.int64)
rna_background_counts = np.bincount(rna_label_codes, minlength=len(rna_label_levels))

rna_matrix = rna[dimension_columns].to_numpy(dtype=np.float32)
atac_matrix = atac[dimension_columns].to_numpy(dtype=np.float32)
atac_clusters = atac["atac_cluster"].to_numpy(dtype=str)

print(f"Embedding dimensions: {len(dimension_columns)}")
print(f"RNA reference cells: {len(rna):,}")
print(f"Named ATAC query cells: {len(atac):,}")
print(f"Excluded Unannotated ATAC cells: {int(((joined['modality'] == 'ATAC') & (joined['atac_cluster'] == UNKNOWN_LABEL)).sum()):,}")
print(f"k values: {', '.join(map(str, K_VALUES))}")

cluster_levels, cluster_cell_counts, counts = build_knn_counts(
    rna_matrix,
    rna_label_codes,
    atac_matrix,
    atac_clusters,
    len(rna_label_levels),
)
composition, enrichment, summary = make_tables(
    cluster_levels,
    cluster_cell_counts,
    counts,
    rna_label_levels,
    rna_background_counts,
)

for path in [composition_csv, enrichment_csv, summary_csv, heatmaps_pdf]:
    Path(path).parent.mkdir(parents=True, exist_ok=True)

composition.to_csv(composition_csv, index=False)
enrichment.to_csv(enrichment_csv, index=False)
summary.to_csv(summary_csv, index=False)
write_heatmaps(enrichment, cluster_levels, rna_label_levels)

print(f"ATAC clusters analyzed: {len(cluster_levels)}")
print(f"RNA level2 labels analyzed: {len(rna_label_levels)}")
print(f"Saved composition table to: {composition_csv}")
print(f"Saved enrichment table to: {enrichment_csv}")
print(f"Saved summary table to: {summary_csv}")
print(f"Saved heatmaps to: {heatmaps_pdf}")
