if "snakemake" in locals():
    scores_csv = snakemake.input["scores"]
    annotations_csv = snakemake.input["annotations"]
    cell_zscores_csv = snakemake.output["cell_zscores"]
    cluster_mean_scores_csv = snakemake.output["cluster_mean_scores"]
    signature_label = snakemake.params.get("signature_label", "marker")
    cluster_subset = snakemake.params.get("cluster_subset", [])
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--scores-csv", default="results/ATAC/level1_marker_ulm_scores.csv")
    parser.add_argument(
        "--annotations-csv",
        default="data/ATAC/label_transfer_RNA_ATAC_V3_clean.csv",
    )
    parser.add_argument(
        "--cell-zscores-csv",
        default="results/ATAC/level1_marker_ulm_cell_zscores.csv",
    )
    parser.add_argument(
        "--cluster-mean-scores-csv",
        default="results/ATAC/level1_marker_ulm_cluster_mean_scores.csv",
    )
    parser.add_argument("--signature-label", default="level1")
    parser.add_argument("--cluster-subset", nargs="*", default=[])
    args = parser.parse_args()

    scores_csv = args.scores_csv
    annotations_csv = args.annotations_csv
    cell_zscores_csv = args.cell_zscores_csv
    cluster_mean_scores_csv = args.cluster_mean_scores_csv
    signature_label = args.signature_label
    cluster_subset = args.cluster_subset

from pathlib import Path

import numpy as np
import pandas as pd


CELL_COLUMN = "cell"
CLUSTER_COLUMN = "cluster"
REQUIRED_ANNOTATION_COLUMNS = {CELL_COLUMN, CLUSTER_COLUMN}


def validate_columns(data, required, label):
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}.")


def read_annotations(path):
    annotations = pd.read_csv(path)
    unnamed_columns = [col for col in annotations.columns if col.startswith("Unnamed:")]
    if CELL_COLUMN not in annotations.columns:
        if len(unnamed_columns) != 1:
            raise ValueError(
                "Annotation table must contain a cell column or exactly one unnamed "
                "cell ID column."
            )
        annotations = annotations.rename(columns={unnamed_columns[0]: CELL_COLUMN})
    else:
        annotations = annotations.drop(columns=unnamed_columns, errors="ignore")

    validate_columns(annotations, REQUIRED_ANNOTATION_COLUMNS, "Annotation table")
    annotations = annotations[[CELL_COLUMN, CLUSTER_COLUMN]].copy()
    annotations[CELL_COLUMN] = annotations[CELL_COLUMN].astype(str)
    annotations[CLUSTER_COLUMN] = annotations[CLUSTER_COLUMN].astype(str)

    if annotations[CELL_COLUMN].duplicated().any():
        raise ValueError("Annotation table contains duplicated cell IDs.")
    return annotations


scores = pd.read_csv(scores_csv)
validate_columns(scores, {CELL_COLUMN}, "ULM score table")
if scores[CELL_COLUMN].duplicated().any():
    raise ValueError("ULM score table contains duplicated cell IDs.")

score_columns = [column for column in scores.columns if column != CELL_COLUMN]
if not score_columns:
    raise ValueError(
        f"ULM score table does not contain any {signature_label} score columns."
    )

for column in score_columns:
    scores[column] = pd.to_numeric(scores[column], errors="coerce")

annotations = read_annotations(annotations_csv)
cluster_order = annotations[CLUSTER_COLUMN].drop_duplicates().tolist()
cluster_subset = [str(cluster) for cluster in cluster_subset]
if cluster_subset:
    missing_clusters = sorted(set(cluster_subset).difference(cluster_order))
    if missing_clusters:
        raise ValueError(
            "Requested ATAC cluster(s) were not found in annotations: "
            f"{', '.join(missing_clusters)}."
        )
    annotations = annotations.loc[annotations[CLUSTER_COLUMN].isin(cluster_subset)].copy()
    cluster_order = [cluster for cluster in cluster_subset if cluster in cluster_order]

joined = annotations.merge(scores, on=CELL_COLUMN, how="inner")
if joined.empty:
    raise ValueError("No overlapping cells were found between annotations and ULM scores.")

score_matrix = joined[score_columns].to_numpy(dtype=float)
row_means = np.nanmean(score_matrix, axis=1, keepdims=True)
row_sds = np.nanstd(score_matrix, axis=1, ddof=1, keepdims=True)
zscores = np.divide(
    score_matrix - row_means,
    row_sds,
    out=np.zeros_like(score_matrix, dtype=float),
    where=np.isfinite(row_sds) & (row_sds != 0),
)

cell_zscores = pd.concat(
    [
        joined[[CELL_COLUMN]].reset_index(drop=True),
        pd.DataFrame(zscores, columns=score_columns),
    ],
    axis=1,
)

cluster_categories = pd.CategoricalDtype(categories=cluster_order, ordered=True)
joined[CLUSTER_COLUMN] = joined[CLUSTER_COLUMN].astype(cluster_categories)

cluster_mean_scores = (
    joined.groupby(CLUSTER_COLUMN, observed=True)[score_columns]
    .mean()
    .reset_index()
)
if "snakemake" not in locals():
    Path(cell_zscores_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(cluster_mean_scores_csv).parent.mkdir(parents=True, exist_ok=True)

cell_zscores.to_csv(cell_zscores_csv, index=False)
cluster_mean_scores.to_csv(cluster_mean_scores_csv, index=False)
