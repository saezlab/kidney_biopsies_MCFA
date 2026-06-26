if "snakemake" in locals():
    mtx_path = snakemake.input["mtx"]
    barcodes_path = snakemake.input["barcodes"]
    features_path = snakemake.input["features"]
    variable_features_path = snakemake.input["variable_features"]
    metadata_csv_path = snakemake.input["metadata"]
    harmony_csv_path = snakemake.input["harmony"]
    h5ad_path = snakemake.output["h5ad"]
else:
    mtx_path = "results/seurat/subset_high_quality/raw_counts.mtx"
    barcodes_path = "results/seurat/subset_high_quality/raw_counts_barcodes.tsv.gz"
    features_path = "results/seurat/subset_high_quality/raw_counts_features.tsv.gz"
    variable_features_path = "results/seurat/subset_high_quality/sct_variable_features.csv"
    metadata_csv_path = "results/seurat/subset_high_quality/metadata_annotations.csv"
    harmony_csv_path = "results/seurat/subset_high_quality/harmony_embeddings.csv"
    h5ad_path = "results/seurat/subset_high_quality/adata_20260717.h5ad"

import os

import anndata as ad
import pandas as pd
from scipy import io


def require_columns(frame, required, label):
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{label} is missing required columns: {', '.join(missing)}."
        )


def require_existing_file(path, label):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} does not exist: {path}")


def main():
    for path, label in (
        (mtx_path, "RNA matrix"),
        (barcodes_path, "RNA barcode table"),
        (features_path, "RNA feature table"),
        (variable_features_path, "SCT variable feature table"),
        (metadata_csv_path, "annotated metadata CSV"),
        (harmony_csv_path, "Harmony embedding CSV"),
    ):
        require_existing_file(path, label)

    print("Reading RNA counts matrix...")
    counts = io.mmread(mtx_path).tocsc()

    print("Reading RNA annotations and embeddings...")
    barcodes = pd.read_csv(barcodes_path, sep="\t", compression="gzip")
    features = pd.read_csv(features_path, sep="\t", compression="gzip")
    variable_features = pd.read_csv(variable_features_path)
    metadata = pd.read_csv(metadata_csv_path)
    harmony = pd.read_csv(harmony_csv_path)

    require_columns(barcodes, ["cellid"], "RNA barcode table")
    require_columns(features, ["feature"], "RNA feature table")
    require_columns(variable_features, ["x"], "SCT variable feature table")
    require_columns(metadata, ["cellid", "level1", "level2"], "annotated metadata CSV")
    require_columns(harmony, ["cellid"], "Harmony embedding CSV")

    barcode_ids = pd.Index(barcodes["cellid"].astype(str), name="cellid")
    feature_ids = pd.Index(features["feature"].astype(str), name="feature")
    variable_feature_ids = pd.Index(variable_features["x"].astype(str))

    if counts.shape[0] != len(feature_ids):
        raise ValueError(
            "RNA counts matrix row count does not match feature count: "
            f"{counts.shape[0]} vs {len(feature_ids)}."
        )
    if counts.shape[1] != len(barcode_ids):
        raise ValueError(
            "RNA counts matrix column count does not match barcode count: "
            f"{counts.shape[1]} vs {len(barcode_ids)}."
        )
    if barcode_ids.has_duplicates:
        raise ValueError("RNA barcode table contains duplicated cell IDs.")
    if feature_ids.has_duplicates:
        raise ValueError("RNA feature table contains duplicated feature names.")
    if variable_feature_ids.has_duplicates:
        raise ValueError("SCT variable feature table contains duplicated genes.")
    if variable_feature_ids.empty:
        raise ValueError("SCT variable feature table is empty.")

    missing_variable_features = variable_feature_ids.difference(feature_ids)
    if len(missing_variable_features) > 0:
        raise ValueError(
            "SCT variable features are missing from the RNA count matrix, including "
            f"{missing_variable_features[:5].tolist()}."
        )

    if metadata["cellid"].duplicated().any():
        raise ValueError("annotated metadata CSV contains duplicated cell IDs.")
    if harmony["cellid"].duplicated().any():
        raise ValueError("Harmony embedding CSV contains duplicated cell IDs.")

    metadata_indexed = metadata.set_index(metadata["cellid"].astype(str))
    metadata_indexed.index.name = "cellid"
    harmony_indexed = harmony.set_index(harmony["cellid"].astype(str))
    harmony_indexed.index.name = "cellid"

    annotated_ids = pd.Index(metadata_indexed.index.astype(str), name="cellid")
    missing_metadata_barcodes = annotated_ids.difference(barcode_ids)
    if len(missing_metadata_barcodes) > 0:
        raise ValueError(
            "Annotated metadata contains cells missing from exported barcodes, including "
            f"{missing_metadata_barcodes[:5].tolist()}."
        )

    missing_harmony = annotated_ids.difference(harmony_indexed.index)
    if len(missing_harmony) > 0:
        raise ValueError(
            "Harmony embeddings are missing annotated cells, including "
            f"{missing_harmony[:5].tolist()}."
        )

    barcode_positions = pd.Series(range(len(barcode_ids)), index=barcode_ids)
    annotated_positions = barcode_positions.loc[annotated_ids].to_numpy()

    print(f"Subsetting to {len(annotated_ids)} annotated cells.")
    counts_annotated = counts[:, annotated_positions].transpose().tocsr()

    obs = metadata_indexed.loc[annotated_ids].copy()
    obs = obs.drop(columns=["cellid"], errors="ignore")

    var = pd.DataFrame(index=feature_ids)
    var["highly_variable"] = var.index.isin(variable_feature_ids)

    rna = ad.AnnData(X=counts_annotated, obs=obs, var=var)

    harmony_aligned = harmony_indexed.loc[annotated_ids].drop(columns=["cellid"], errors="ignore")
    rna.obsm["X_harmony"] = harmony_aligned.to_numpy()

    if "level1" not in rna.obs or "level2" not in rna.obs:
        raise ValueError("AnnData obs is missing level1/level2 annotations.")
    if rna.X is None:
        raise ValueError("AnnData X is missing the RNA counts matrix.")
    expected_shape = (len(annotated_ids), len(feature_ids))
    if rna.X.shape != expected_shape:
        raise ValueError(
            "AnnData X shape does not match annotated cells and features: "
            f"{rna.X.shape} vs {expected_shape}."
        )
    if "X_harmony" not in rna.obsm:
        raise ValueError("AnnData is missing obsm['X_harmony'].")

    print(f"Writing AnnData to {h5ad_path}...")
    rna.write_h5ad(h5ad_path)


if __name__ == "__main__":
    main()
