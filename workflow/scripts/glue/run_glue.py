if "snakemake" in locals():
    rna_h5ad_path = snakemake.input["rna_h5ad"]
    atac_h5ad_path = snakemake.input["atac_h5ad"]
    atac_lsi_csv_path = snakemake.input["atac_lsi_csv"]
    gtf_path = snakemake.input["gtf"]
    out_csv_path = snakemake.output["csv"]
    model_dir = snakemake.output["model_dir"]
    cuda_device = int(snakemake.params["cuda_device"])
else:
    rna_h5ad_path = "results/seurat/subset_high_quality/adata_20260717.h5ad"
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

import numpy as np
import pandas as pd
import scanpy as sc
import scglue
import torch


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
    rna = sc.read_h5ad(rna_h5ad_path)

    if rna.X is None:
        raise ValueError("RNA AnnData X is missing the counts matrix.")
    if "X_harmony" not in rna.obsm:
        raise ValueError("RNA AnnData is missing obsm['X_harmony'].")
    require_columns(rna.obs, ["level1", "level2"], "RNA AnnData obs")
    require_columns(rna.var, ["highly_variable"], "RNA AnnData var")

    if rna.obs_names.has_duplicates:
        raise ValueError("RNA AnnData contains duplicated cell IDs.")
    if rna.var_names.has_duplicates:
        raise ValueError("RNA AnnData contains duplicated feature names.")
    if rna.n_obs == 0 or rna.n_vars == 0:
        raise ValueError("RNA AnnData is empty.")
    if int(rna.var["highly_variable"].sum()) == 0:
        raise ValueError("RNA AnnData has no highly variable genes.")
    if rna.obsm["X_harmony"].shape[0] != rna.n_obs:
        raise ValueError(
            "RNA Harmony embedding row count does not match AnnData observations."
        )

    rna.X = stochastic_round_sparse_counts(rna.X.tocsr())

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
        (rna_h5ad_path, "RNA AnnData"),
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
        fit_kws={"directory": model_dir},
        init_kws={"latent_dim": 30}
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
