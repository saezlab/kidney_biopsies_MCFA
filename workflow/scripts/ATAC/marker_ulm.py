if "snakemake" in locals():
    h5ad_path = snakemake.input["h5ad"]
    markers_csv = snakemake.input["markers_csv"]
    scores_csv = snakemake.output["scores"]
    pvals_csv = snakemake.output["pvals"]
    signature_label = snakemake.params.get("signature_label", "marker")
else:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", default="data/ATAC/gene_score_matrix.h5ad")
    parser.add_argument(
        "--markers-csv",
        default="results/seurat/subset_high_quality/markers/findallmarkers_level1.csv",
    )
    parser.add_argument(
        "--scores-csv",
        default="results/ATAC/level1_marker_ulm_scores.csv",
    )
    parser.add_argument(
        "--pvals-csv",
        default="results/ATAC/level1_marker_ulm_pvals.csv",
    )
    parser.add_argument("--signature-label", default="level1")
    args = parser.parse_args()

    h5ad_path = args.h5ad
    markers_csv = args.markers_csv
    scores_csv = args.scores_csv
    pvals_csv = args.pvals_csv
    signature_label = args.signature_label

from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd


MARKER_P_ADJ_CUTOFF = 0.05
MIN_TARGETS_PER_SOURCE = 5
REQUIRED_MARKER_COLUMNS = {"p_val_adj", "avg_log2FC", "cluster", "gene"}


def validate_columns(data, required, label):
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}.")


def build_marker_network(markers, genes):
    validate_columns(markers, REQUIRED_MARKER_COLUMNS, "Marker table")

    marker_net = markers.loc[
        markers["p_val_adj"] < MARKER_P_ADJ_CUTOFF,
        ["cluster", "gene", "avg_log2FC"],
    ].copy()
    marker_net = marker_net.rename(
        columns={"cluster": "source", "gene": "target", "avg_log2FC": "weight"}
    )
    marker_net["source"] = marker_net["source"].astype(str).str.strip()
    marker_net["target"] = marker_net["target"].astype(str).str.strip()
    marker_net["weight"] = pd.to_numeric(marker_net["weight"], errors="coerce")
    marker_net = marker_net.replace({"": np.nan})
    marker_net = marker_net.replace([np.inf, -np.inf], np.nan)
    marker_net = marker_net.dropna(subset=["source", "target", "weight"])

    if marker_net.empty:
        raise ValueError(
            "No marker rows remained after filtering to "
            f"p_val_adj < {MARKER_P_ADJ_CUTOFF}."
        )

    marker_net = (
        marker_net.groupby(["source", "target"], as_index=False)["weight"]
        .mean()
        .sort_values(["source", "target"])
    )
    marker_net = marker_net.loc[marker_net["target"].isin(genes)].copy()

    if marker_net.empty:
        raise ValueError(
            "No significant marker genes overlap the ATAC gene score matrix features."
        )

    target_counts = marker_net.groupby("source")["target"].nunique()
    low_target_sources = target_counts[target_counts < MIN_TARGETS_PER_SOURCE]
    if not low_target_sources.empty:
        formatted = ", ".join(
            f"{source} ({count})" for source, count in low_target_sources.items()
        )
        raise ValueError(
            f"Some {signature_label} marker sets have fewer than "
            f"{MIN_TARGETS_PER_SOURCE} overlapping genes: {formatted}."
        )

    return marker_net


adata = ad.read_h5ad(h5ad_path)
markers = pd.read_csv(markers_csv)
net = build_marker_network(markers, set(adata.var_names))

ulm_result = dc.run_ulm(
    mat=adata,
    net=net,
    source="source",
    target="target",
    weight="weight",
    min_n=MIN_TARGETS_PER_SOURCE,
    use_raw=False,
)
if ulm_result is None:
    scores = adata.obsm["ulm_estimate"]
    pvals = adata.obsm["ulm_pvals"]
else:
    scores, pvals = ulm_result

scores = scores.rename_axis("cell").reset_index()
pvals = pvals.rename_axis("cell").reset_index()

if "snakemake" not in locals():
    Path(scores_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(pvals_csv).parent.mkdir(parents=True, exist_ok=True)

scores.to_csv(scores_csv, index=False)
pvals.to_csv(pvals_csv, index=False)
