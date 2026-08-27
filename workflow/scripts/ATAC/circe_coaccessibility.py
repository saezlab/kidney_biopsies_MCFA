if "snakemake" in locals():
    h5ad_path = snakemake.input["h5ad"]
    lsi_csv_path = snakemake.input["lsi"]
    coaccessibility_csv_path = snakemake.output["coaccessibility"]
    peaks_csv_path = snakemake.output["peaks"]
    cluster_subset = list(snakemake.params["clusters"])
    min_cells_per_peak = int(snakemake.params["min_cells_per_peak"])
    min_peaks_per_cell = int(snakemake.params["min_peaks_per_cell"])
    njobs = int(snakemake.threads)
else:
    h5ad_path = "data/ATAC/peak_mat2.h5ad"
    lsi_csv_path = "data/ATAC/harmony_lsi.csv"
    coaccessibility_csv_path = "results/ATAC/circe_coaccessibility.csv"
    peaks_csv_path = "results/ATAC/circe_peaks.csv"
    cluster_subset = ["pt_C3", "pt_C6"]
    min_cells_per_peak = 500
    min_peaks_per_cell = 20
    njobs = 8

import os

import circe as ci
import numpy as np
import pandas as pd
import scanpy as sc


CLUSTER_COLUMN = "cluster"
LSI_KEY = "X_lsi"


if __name__ == "__main__":
    # Validate the fixed pipeline inputs before loading the large peak matrix.
    for path, label in (
        (h5ad_path, "ATAC AnnData"),
        (lsi_csv_path, "ATAC LSI CSV"),
    ):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} does not exist: {path}")

    # Load ATAC peaks and use the raw count layer that is present in peak_mat2.h5ad.
    print("Loading ATAC peak matrix.", flush=True)
    atac = sc.read_h5ad(h5ad_path)
    if "counts" not in atac.layers:
        raise ValueError("ATAC AnnData is missing required layer 'counts'.")
    atac.X = atac.layers["counts"]

    atac.obs_names = atac.obs_names.astype(str)
    atac.var_names = atac.var_names.astype(str)
    if atac.obs_names.has_duplicates:
        examples = atac.obs_names[atac.obs_names.duplicated()][:5].tolist()
        raise ValueError(
            f"ATAC AnnData contains duplicated cell IDs, including {examples}."
        )
    if atac.var_names.has_duplicates:
        examples = atac.var_names[atac.var_names.duplicated()][:5].tolist()
        raise ValueError(
            f"ATAC AnnData contains duplicated peak IDs, including {examples}."
        )
    if atac.n_obs == 0 or atac.n_vars == 0:
        raise ValueError("ATAC AnnData is empty.")
    if CLUSTER_COLUMN not in atac.obs:
        raise ValueError(f"ATAC AnnData is missing obs['{CLUSTER_COLUMN}'].")

    # Load the Harmony LSI table, align it to the H5AD cells, and ignore extra rows.
    print(f"Loading LSI embedding from {lsi_csv_path}.", flush=True)
    lsi = pd.read_csv(lsi_csv_path, index_col=0)
    lsi.index = lsi.index.astype(str)
    if lsi.index.duplicated().any():
        examples = lsi.index[lsi.index.duplicated()].unique()[:5].tolist()
        raise ValueError(
            f"ATAC LSI table contains duplicated cell IDs, including {examples}."
        )

    atac_obs_names = pd.Index(atac.obs_names)
    missing_lsi = atac_obs_names.difference(lsi.index)
    if len(missing_lsi) > 0:
        raise ValueError(
            "ATAC LSI embeddings are missing ATAC barcodes, including "
            f"{missing_lsi[:5].tolist()}."
        )

    extra_lsi = lsi.index.difference(atac_obs_names)
    if len(extra_lsi) > 0:
        print(
            f"Ignoring {len(extra_lsi)} ATAC LSI rows that are not present "
            "in the ATAC peak matrix.",
            flush=True,
        )

    atac.obsm[LSI_KEY] = lsi.reindex(atac_obs_names).to_numpy()

    # Keep the PT/aPT high-confidence clusters for CIRCE linking.
    clusters = [str(cluster) for cluster in cluster_subset]
    keep_cells = atac.obs[CLUSTER_COLUMN].astype(str).isin(clusters)
    if int(keep_cells.sum()) == 0:
        raise ValueError(
            "No ATAC cells matched requested clusters: "
            f"{', '.join(clusters)}."
        )
    print(
        f"Retaining {int(keep_cells.sum())} ATAC cells in {', '.join(clusters)}.",
        flush=True,
    )
    atac = atac[keep_cells, :].copy()

    # CIRCE expects peak names in chrom_start_end format.
    split_peaks = atac.var_names.to_series().str.split("_")
    invalid_peaks = split_peaks.map(len) != 3
    if invalid_peaks.any():
        examples = atac.var_names[invalid_peaks][:5].tolist()
        raise ValueError(
            "ATAC peak names must follow 'chrom_start_end', but invalid values "
            f"include: {examples}."
        )

    # Follow CIRCE's filtering guidance after cluster subsetting to avoid
    # fitting graphical-lasso models over hundreds of thousands of rare peaks.
    print(
        f"Filtering peaks with min_cells={min_cells_per_peak} and "
        f"cells with min_genes={min_peaks_per_cell}.",
        flush=True,
    )
    sc.pp.filter_genes(atac, min_cells=min_cells_per_peak)
    sc.pp.filter_cells(atac, min_genes=min_peaks_per_cell)
    sc.pp.filter_genes(atac, min_counts=1)
    if atac.n_obs == 0 or atac.n_vars == 0:
        raise ValueError("No ATAC data retained after CIRCE preprocessing filters.")
    print(
        f"Retained {atac.n_obs} cells and {atac.n_vars} peaks after filtering.",
        flush=True,
    )

    ci.add_region_infos(atac, inplace=True)

    # Build the co-accessibility network in CIRCE's metacell representation.
    print("Computing CIRCE metacells.", flush=True)
    metacells = ci.metacells.compute_metacells(atac, dim_reduction=LSI_KEY)
    if metacells.n_obs == 0 or metacells.n_vars == 0:
        raise ValueError("CIRCE metacell AnnData is empty.")
    print(
        f"Computed {metacells.n_obs} metacells across {metacells.n_vars} peaks.",
        flush=True,
    )

    pd.DataFrame({"peak_id": metacells.var_names.astype(str)}).to_csv(
        peaks_csv_path,
        index=False,
    )
    print(f"Wrote retained CIRCE peaks to {peaks_csv_path}.", flush=True)

    print("Computing CIRCE peak co-accessibility network.", flush=True)
    ci.compute_atac_network(
        metacells,
        njobs=njobs,
        verbose=2,
        organism="human",
    )

    # Extract links using CIRCE's Peak1, Peak2, score output columns.
    print("Extracting CIRCE co-accessibility links.", flush=True)
    coaccessibility = ci.extract_atac_links(metacells)
    required_coaccessibility_columns = ["Peak1", "Peak2", "score"]
    missing_columns = [
        column
        for column in required_coaccessibility_columns
        if column not in coaccessibility.columns
    ]
    if missing_columns:
        raise ValueError(
            "CIRCE co-accessibility links are missing required columns: "
            f"{', '.join(missing_columns)}."
        )
    coaccessibility = coaccessibility.loc[:, required_coaccessibility_columns].copy()
    if coaccessibility.empty:
        raise ValueError("CIRCE produced an empty co-accessibility network.")

    coaccessibility.to_csv(coaccessibility_csv_path, index=False)
    print(f"Wrote co-accessibility network to {coaccessibility_csv_path}.", flush=True)
