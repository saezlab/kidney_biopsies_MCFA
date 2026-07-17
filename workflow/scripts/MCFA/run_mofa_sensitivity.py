# 09.07.2026 Charlotte Boys

# %% Import libraries
import pandas as pd
import liana as li
import muon as mu


# %% Parse snakemake object
if "snakemake" in locals():
    h5mu_input = snakemake.input["h5mu_input"]
    model_output = snakemake.output["model_output"]
    scores_output = snakemake.output["scores_output"]
    scores_wide_output = snakemake.output["scores_wide_output"]
    loadings_output = snakemake.output["loadings_output"]
    nfactors = int(snakemake.params["nfactors"])
    hvg = snakemake.params["hvg"]
    exclude_samples = list(snakemake.params["exclude_samples"])
else:
    h5mu_input = "results/MCFA/adata.h5mu"
    model_output = "results/MCFA/model~all~5~exclude-TIN.h5ad"
    scores_output = "results/MCFA/scores~all~5~exclude-TIN.csv"
    scores_wide_output = "results/MCFA/widescores~all~5~exclude-TIN.csv"
    loadings_output = "results/MCFA/loadings~all~5~exclude-TIN.csv"
    nfactors = 5
    hvg = "all"
    exclude_samples = ["K12", "K14", "K16", "K17", "K19", "K25", "K26"]

# %% Read mdata
print("INFO: Loading assay")
mdata = mu.read_h5mu(h5mu_input)

# %% Exclude samples
exclude_samples = [str(sample) for sample in exclude_samples]
exclude_set = set(exclude_samples)
present_samples = set(map(str, mdata.obs_names))
samples_to_drop = sorted(exclude_set & present_samples)
missing_samples = sorted(exclude_set - present_samples)

print(f"INFO: Requested exclusion of {len(exclude_samples)} samples using {sample_key}")
if missing_samples:
    print(f"INFO: {len(missing_samples)} excluded samples were not present: {missing_samples}")
print(f"INFO: Dropping {len(samples_to_drop)} samples present in the MuData object: {samples_to_drop}")

for view in list(mdata.mod.keys()):
    adata = mdata.mod[view]
    # mdata built such that mdata.mod[view].obs_names is the "sample" column
    keep_mask = ~adata.obs_names.astype(str).isin(exclude_set)
    print(f"INFO: {view} keeping {keep_mask.sum()}/{adata.n_obs} samples after exclusion")
    mdata.mod[view] = adata[keep_mask, :].copy()

mdata.update()
if mdata.n_obs == 0:
    raise ValueError("No samples remain after applying the MOFA sensitivity exclusion list")

# %% Set up to run with HVGs
if hvg == "all":
    use_var = None
elif hvg == "hvg":
    use_var = None
    for view in list(mdata.mod.keys()):
        adata = mdata.mod[view]
        if "highly_variable" not in adata.var:
            raise ValueError(f"Missing highly_variable column for view {view}")

        hvg_mask = adata.var["highly_variable"].astype(bool).to_numpy()
        if hvg_mask.sum() == 0:
            raise ValueError(f"No highly variable genes retained for view {view}")

        print(f"INFO: {view} keeping {hvg_mask.sum()}/{adata.n_vars} HVGs")
        mdata.mod[view] = adata[:, hvg_mask].copy()
    mdata.update()
else:
    raise ValueError(f"Unsupported hvg value: {hvg}")

# %% Run mofa
print("INFO: Running MOFA")
mu.tl.mofa(
    mdata,
    groups_label=None,
    use_obs="union",
    convergence_mode="slow",
    n_factors=nfactors,
    seed=1337,
    outfile=model_output,
    use_var=use_var,
)
print(mdata)
print(mdata.varm)

# %% Extract factor scores and loadings
print("INFO: Extracting MOFA factor scores and loadings")
factor_scores_obs = li.utils.get_factor_scores(mdata, obsm_key="X_mofa")

print(factor_scores_obs.head())
factor_scores_wide = factor_scores_obs.set_index("index").filter(like="Factor", axis=1)
factor_scores_long = (
    factor_scores_wide.rename_axis("Sample")
    .reset_index()
    .melt(id_vars=["Sample"], var_name="Factor", value_name="Score")
)

variable_loadings = li.utils.get_variable_loadings(mdata, varm_key="LFs", view_sep=":")
print(variable_loadings.head())
variable_loadings_long = pd.melt(
    variable_loadings,
    id_vars=["view", "variable"],
    var_name="Factor",
    value_name="Loading",
)

# %% Save factor scores and loadings
print("INFO: Saving MOFA factor scores and loadings")
factor_scores_long.to_csv(scores_output)
factor_scores_wide.to_csv(scores_wide_output)
variable_loadings_long.to_csv(loadings_output)
