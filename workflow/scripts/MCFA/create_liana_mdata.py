# 26.04.2023 - Charlotte Boys
# Updated 09.07.2026 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import liana as li


# %% Parse snakemake object
if "snakemake" in locals():
    h5ad_input = snakemake.input["h5ad_input"]
    h5mu_output = snakemake.output["h5mu_output"]
    sample_key = snakemake.params["sample_key"]
else:
    h5ad_input = "results/MCFA/adata_LIANA.h5ad"
    h5mu_output = "results/MCFA/adata_LIANA.h5mu"
    sample_key = "sample"

# %% Read in AnnData object
print("INFO: Loading data")
adata = sc.read_h5ad(h5ad_input)

# %%
print("INFO: Creating LIANA multi-view structure")
mdata = li.multi.lrs_to_views(
    adata,
    sample_key=sample_key,
    score_key="magnitude_rank",
    obs_keys=[
        "disease",
        "sex",
        "IFTA",
        "corticosteroids",
        "classification",
        "jm_slope",
        "immediate_slope",
    ],
    lr_prop=0.1,
    lrs_per_sample=20,
    lrs_per_view=20,
    samples_per_view=10,
    min_variance=0,
    lr_fill=0,
    verbose=True,
)

print(mdata)

if "snakemake" in locals():
    print(f"INFO: Writing MuData to {h5mu_output}")
    mdata.write(h5mu_output)
