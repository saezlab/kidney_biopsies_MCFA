if "snakemake" in locals():
    mtx_path = snakemake.input["mtx"]
    barcodes_path = snakemake.input["barcodes"]
    features_path = snakemake.input["features"]
    variable_features_path = snakemake.input["variable_features"]
    harmony_csv_path = snakemake.input["harmony_csv"]
    metadata_csv_path = snakemake.input["metadata_csv"]
    atac_h5ad_path = snakemake.input["atac_h5ad"]
    atac_lsi_csv_path = snakemake.input["atac_lsi_csv"]
    gtf_path = snakemake.input["gtf"]
    out_csv_path = snakemake.output["csv"]
    model_dir = snakemake.output["model_dir"]
    cuda_device = int(snakemake.params["cuda_device"])
else:
    mtx_path = "results/seurat/raw_counts.mtx"
    barcodes_path = "results/seurat/raw_counts_barcodes.tsv.gz"
    features_path = "results/seurat/raw_counts_features.tsv.gz"
    variable_features_path = "results/seurat/sct_variable_features.csv"
    harmony_csv_path = "results/seurat/harmony_embeddings.csv"
    metadata_csv_path = "results/seurat/metadata.csv"
    atac_h5ad_path = "data/ATAC/peak_mat2.h5ad"
    atac_lsi_csv_path = "data/ATAC/harmony_lsi.csv"
    gtf_path = "data/reference/gencode.v38.primary_assembly.annotation.gtf.gz"
    out_csv_path = "results/glue/coembedding.csv"
    model_dir = "results/glue/model"
    cuda_device = 0

import os
import sys

os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_device)


def expose_nixos_nvidia_driver():
    driver_dir = "/run/opengl-driver/lib"
    reexec_marker = "GLUE_NVIDIA_DRIVER_REEXEC"
    if not os.path.exists(os.path.join(driver_dir, "libcuda.so.1")):
        return
    if os.environ.get(reexec_marker) == "1":
        return

    library_paths = [
        path for path in os.environ.get("LD_LIBRARY_PATH", "").split(":")
        if path
    ]
    if driver_dir in library_paths:
        return

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = ":".join([driver_dir, *library_paths])
    env[reexec_marker] = "1"
    os.execve(sys.executable, [sys.executable, *sys.argv], env)


expose_nixos_nvidia_driver()

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scglue
import torch
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


def stochastic_round_sparse_counts(counts, seed=0, chunk_size=1_000_000):
    values = counts.data
    noninteger_count = 0
    for start in range(0, values.size, chunk_size):
        stop = min(start + chunk_size, values.size)
        chunk = values[start:stop]
        if not np.isfinite(chunk).all():
            raise ValueError("RNA counts matrix contains non-finite values.")
        if np.any(chunk < 0):
            raise ValueError("RNA counts matrix contains negative values.")
        noninteger_count += int(np.count_nonzero(chunk != np.floor(chunk)))

    if noninteger_count == 0:
        return counts

    print(
        f"Stochastically rounding {noninteger_count} fractional SoupX values "
        f"to integers with seed {seed}."
    )
    rng = np.random.default_rng(seed)
    for start in range(0, values.size, chunk_size):
        stop = min(start + chunk_size, values.size)
        chunk = values[start:stop]
        floors = np.floor(chunk)
        chunk[:] = floors + rng.binomial(1, chunk - floors)

    counts.eliminate_zeros()
    for start in range(0, counts.data.size, chunk_size):
        stop = min(start + chunk_size, counts.data.size)
        chunk = counts.data[start:stop]
        if np.any(chunk != np.floor(chunk)):
            raise ValueError(
                "RNA stochastic rounding did not produce integer values."
            )
    return counts


def check_cuda():
    print(f"Requested physical GPU: {cuda_device}")
    print(f"CUDA_VISIBLE_DEVICES: {os.environ['CUDA_VISIBLE_DEVICES']}")
    print(f"LD_LIBRARY_PATH: {os.environ.get('LD_LIBRARY_PATH', '')}")
    print(f"PyTorch version: {torch.__version__}")
    print(f"PyTorch CUDA runtime: {torch.version.cuda}")

    if torch.version.cuda is None:
        raise RuntimeError(
            "The GLUE conda environment contains a CPU-only PyTorch build. "
            "Recreate the environment from workflow/envs/glue.yaml."
        )
    if not torch.cuda.is_available():
        raise RuntimeError(
            "GLUE in this workflow currently requires a CUDA-capable GPU, "
            "but torch.cuda.is_available() returned False. Verify that the "
            "NVIDIA driver is loaded, CUDA_VISIBLE_DEVICES selects a visible GPU, "
            "and the conda environment was rebuilt from workflow/envs/glue.yaml."
        )
    print(f"CUDA device count: {torch.cuda.device_count()}")
    print(f"Logical CUDA device 0: {torch.cuda.get_device_name(0)}")

    try:
        probe = torch.zeros(1, device="cuda")
        del probe
    except Exception as error:
        raise RuntimeError(
            "PyTorch detected CUDA but could not allocate a tensor on the selected "
            f"GPU {cuda_device}."
        ) from error
    print("CUDA tensor allocation succeeded.")


def load_rna():
    counts = io.mmread(mtx_path).tocsr()
    counts = stochastic_round_sparse_counts(counts)

    barcodes = pd.read_csv(barcodes_path, sep="\t", compression="gzip")
    features = pd.read_csv(features_path, sep="\t", compression="gzip")
    variable_features = pd.read_csv(variable_features_path)
    harmony = pd.read_csv(harmony_csv_path)
    metadata = pd.read_csv(metadata_csv_path)

    require_columns(barcodes, ["cellid"], "RNA barcode table")
    require_columns(features, ["feature"], "RNA feature table")
    require_columns(variable_features, ["x"], "SCT variable feature table")
    require_columns(harmony, ["cellid"], "RNA harmony table")
    require_columns(metadata, ["cellid"], "RNA metadata table")

    barcode_ids = pd.Index(barcodes["cellid"])
    feature_ids = pd.Index(features["feature"])
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
    if harmony["cellid"].duplicated().any():
        raise ValueError("RNA harmony table contains duplicated cell IDs.")
    if metadata["cellid"].duplicated().any():
        raise ValueError("RNA metadata table contains duplicated cell IDs.")

    rna = ad.AnnData(
        X=counts.transpose().tocsr(),
        obs=pd.DataFrame(index=barcode_ids),
        var=pd.DataFrame(index=feature_ids),
    )
    rna.layers["counts"] = rna.X.copy()
    rna.var["highly_variable"] = rna.var_names.isin(variable_feature_ids)

    harmony_indexed = harmony.set_index("cellid")
    missing_harmony = barcode_ids.difference(harmony_indexed.index)
    extra_harmony = harmony_indexed.index.difference(barcode_ids)
    if len(missing_harmony) > 0:
        raise ValueError(
            "RNA harmony embeddings are missing exported barcodes, including "
            f"{missing_harmony[:5].tolist()}."
        )
    if len(extra_harmony) > 0:
        raise ValueError(
            "RNA harmony embeddings contain unexpected barcodes, including "
            f"{extra_harmony[:5].tolist()}."
        )
    harmony_aligned = harmony_indexed.reindex(rna.obs_names)
    rna.obsm["X_harmony"] = harmony_aligned.to_numpy()

    metadata_indexed = metadata.set_index("cellid")
    missing_metadata = barcode_ids.difference(metadata_indexed.index)
    extra_metadata = metadata_indexed.index.difference(barcode_ids)
    if len(missing_metadata) > 0:
        raise ValueError(
            "RNA metadata is missing exported barcodes, including "
            f"{missing_metadata[:5].tolist()}."
        )
    if len(extra_metadata) > 0:
        raise ValueError(
            "RNA metadata contains unexpected barcodes, including "
            f"{extra_metadata[:5].tolist()}."
        )
    metadata_aligned = metadata_indexed.reindex(rna.obs_names)
    if metadata_aligned.isnull().all(axis=1).any():
        raise ValueError("RNA metadata alignment introduced missing rows.")
    rna.obs = metadata_aligned

    return rna


def add_rna_gene_annotations(rna):
    scglue.data.get_gene_annotation(rna, gtf=gtf_path, gtf_by="gene_name")

    required_columns = ["chrom", "chromStart", "chromEnd"]
    missing_coords = rna.var.loc[:, required_columns].isnull().any(axis=1)
    rna = rna[:, ~missing_coords].copy()
    if rna.n_vars == 0:
        raise ValueError("No RNA genes retained after adding genomic coordinates.")
    retained_variable_features = int(rna.var["highly_variable"].sum())
    if retained_variable_features == 0:
        raise ValueError(
            "No SCT variable features retained after adding genomic coordinates."
        )
    print(
        f"Retained {retained_variable_features} SCT variable features after "
        "GENCODE coordinate filtering."
    )
    return rna


def load_atac():
    atac = sc.read_h5ad(atac_h5ad_path)
    if "counts" in atac.layers:
        atac.X = atac.layers["counts"]

    lsi = pd.read_csv(atac_lsi_csv_path, index_col=0)
    lsi.index = lsi.index.astype(str)
    atac_obs_names = pd.Index(atac.obs_names.astype(str))

    if lsi.index.duplicated().any():
        raise ValueError("ATAC LSI table contains duplicated cell IDs.")

    missing_lsi = atac_obs_names.difference(lsi.index)
    extra_lsi = lsi.index.difference(atac_obs_names)
    if len(missing_lsi) > 0:
        raise ValueError(
            "ATAC LSI embeddings are missing ATAC barcodes, including "
            f"{missing_lsi[:5].tolist()}."
        )
    if len(extra_lsi) > 0:
        print(
            f"Ignoring {len(extra_lsi)} ATAC LSI rows that are not present "
            "in the ATAC peak matrix."
        )

    lsi_aligned = lsi.reindex(atac_obs_names)
    atac.obsm["X_lsi"] = lsi_aligned.to_numpy()

    split = atac.var_names.to_series().str.split("_")
    invalid = split.map(len) != 3
    if invalid.any():
        examples = atac.var_names[invalid][:5].tolist()
        raise ValueError(
            "ATAC peak names must follow 'chrom_start_end', but invalid values "
            f"include: {examples}."
        )

    atac.var["chrom"] = split.map(lambda x: x[0]).to_numpy()
    atac.var["chromStart"] = split.map(lambda x: int(x[1])).to_numpy()
    atac.var["chromEnd"] = split.map(lambda x: int(x[2])).to_numpy()
    return atac


def main():
    for path, label in (
        (mtx_path, "RNA matrix"),
        (barcodes_path, "RNA barcode table"),
        (features_path, "RNA feature table"),
        (variable_features_path, "SCT variable feature table"),
        (harmony_csv_path, "RNA harmony CSV"),
        (metadata_csv_path, "RNA metadata CSV"),
        (atac_h5ad_path, "ATAC AnnData"),
        (atac_lsi_csv_path, "ATAC harmony LSI CSV"),
        (gtf_path, "GENCODE annotation"),
    ):
        require_existing_file(path, label)

    check_cuda()

    print("Loading RNA data...")
    rna = load_rna()
    print("Adding RNA gene annotations...")
    rna = add_rna_gene_annotations(rna)

    print("Loading ATAC data...")
    atac = load_atac()

    print("Constructing guidance graph...")
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])

    print("Configuring datasets...")
    scglue.models.configure_dataset(
        rna, "NB",
        use_layer="counts",
        use_rep="X_harmony"
    )
    scglue.models.configure_dataset(
        atac, "NB",
        use_rep="X_lsi"
    )

    print("Fitting SCGLUE model...")
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac},
        guidance,
        fit_kws={"directory": model_dir}
    )
    glue.save(os.path.join(model_dir, "glue.dill"))

    print("Encoding coembedding...")
    x_rna = pd.DataFrame(
        glue.encode_data("rna", rna),
        index=["rna_" + obs_name for obs_name in rna.obs_names]
    )
    x_atac = pd.DataFrame(
        glue.encode_data("atac", atac),
        index=["atac_" + obs_name for obs_name in atac.obs_names]
    )
    coembedding = pd.concat([x_rna, x_atac], axis=0)
    coembedding.to_csv(out_csv_path)
    print(f"Wrote coembedding to {out_csv_path}")


if __name__ == "__main__":
    main()
