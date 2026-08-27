if "snakemake" in locals():
    coaccessibility_csv_path = snakemake.input["coaccessibility"]
    peaks_csv_path = snakemake.input["peaks"]
    peak_to_gene_csv_path = snakemake.output["peak_to_gene"]
else:
    coaccessibility_csv_path = "results/ATAC/circe_coaccessibility.csv"
    peaks_csv_path = "results/ATAC/circe_peaks.csv"
    peak_to_gene_csv_path = "results/ATAC/circe_peak_to_gene.csv"

import os

import pandas as pd
from celloracle import motif_analysis as ma


# Validate the CIRCE outputs before running CellOracle TSS annotation.
for path, label in (
    (coaccessibility_csv_path, "CIRCE co-accessibility CSV"),
    (peaks_csv_path, "CIRCE peaks CSV"),
):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} does not exist: {path}")


# Load CIRCE links and retained peak names from the upstream CIRCE container step.
coaccessibility = pd.read_csv(coaccessibility_csv_path)
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
    raise ValueError("CIRCE co-accessibility network is empty.")

peaks = pd.read_csv(peaks_csv_path)
if "peak_id" not in peaks.columns:
    raise ValueError("CIRCE peaks CSV is missing required column: peak_id.")
peaks = peaks["peak_id"].dropna().astype(str)
if peaks.empty:
    raise ValueError("CIRCE peaks CSV is empty.")
if peaks.duplicated().any():
    examples = peaks[peaks.duplicated()].unique()[:5].tolist()
    raise ValueError(f"CIRCE peaks CSV contains duplicated peak IDs, including {examples}.")


# Use CellOracle to connect co-accessible peaks to gene TSSs.
print("Annotating TSS peaks with CellOracle ref_genome=hg38.")
tss_annotated = ma.get_tss_info(
    peak_str_list=peaks.to_numpy(),
    ref_genome="hg38",
)
if tss_annotated.empty:
    raise ValueError("CellOracle did not identify any TSS-overlapping peaks.")

print("Integrating TSS peaks with CIRCE co-accessibility links.")
celloracle_connections = coaccessibility.rename(columns={"score": "coaccess"}).loc[
    :, ["Peak1", "Peak2", "coaccess"]
]
peak_to_gene = ma.integrate_tss_peak_with_cicero(
    tss_peak=tss_annotated,
    cicero_connections=celloracle_connections,
)
required_peak_to_gene_columns = ["peak_id", "gene_short_name", "coaccess"]
missing_columns = [
    column
    for column in required_peak_to_gene_columns
    if column not in peak_to_gene.columns
]
if missing_columns:
    raise ValueError(
        "CellOracle peak-to-gene network is missing required columns: "
        f"{', '.join(missing_columns)}."
    )
if peak_to_gene.empty:
    raise ValueError("CellOracle produced an empty peak-to-gene network.")

peak_to_gene.to_csv(peak_to_gene_csv_path, index=False)
print(f"Wrote peak-to-gene network to {peak_to_gene_csv_path}.")
