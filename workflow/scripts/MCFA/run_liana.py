# 26.04.2023 - Charlotte Boys
# Updated 09.07.2026 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import liana as li


# %% Parse snakemake object
if "snakemake" in locals():
    h5ad_input = snakemake.input["h5ad_input"]
    h5ad_output = snakemake.output["h5ad_output"]
    sample_key = snakemake.params["sample_key"]
    groupby = snakemake.params["groupby"]
else:
    h5ad_input = "results/MCFA/adata.h5ad"
    h5ad_output = "results/MCFA/adata_LIANA.h5ad"
    sample_key = "sample"
    groupby = "level1"

# %% Read in AnnData object
print("INFO: Loading data")
adata = sc.read_h5ad(h5ad_input)

# %% Normalize and log-transform for LIANA scoring
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# %% Run LIANA
print("INFO: Running LIANA")
li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    sample_key=sample_key,
    expr_prop=0.1,
    use_raw=False,
    n_perms=100,
    return_all_lrs=False,
    verbose=True,
)

print(adata.uns["liana_res"].sort_values("magnitude_rank").head())

if "snakemake" in locals():
    print(f"INFO: Writing AnnData to {h5ad_output}")
    adata.write(h5ad_output)
