# 26.04.2024 Charlotte Boys
# Updated 09.07.2026 - Charlotte Boys

# %% Import libraries
import muon as mu


def strip_view_prefix(feature_name):
    feature_name = str(feature_name)
    if ":" in feature_name:
        return feature_name.split(":", 1)[1]
    return feature_name


def complex_genes(complex_name):
    return set(filter(None, strip_view_prefix(complex_name).split("_")))


def interaction_genes(interaction_name):
    interaction_name = strip_view_prefix(interaction_name)
    if "^" not in interaction_name:
        return set()

    ligand, receptor = interaction_name.split("^", 1)
    return complex_genes(ligand) | complex_genes(receptor)


# %% Parse snakemake object
if "snakemake" in locals():
    h5mu_input = snakemake.input["h5mu_input"]
    h5mu_liana = snakemake.input["h5mu_liana"]
    h5mu_output = snakemake.output["h5mu_output"]
else:
    h5mu_input = "results/MCFA/adata.h5mu"
    h5mu_liana = "results/MCFA/adata_LIANA.h5mu"
    h5mu_output = "results/MCFA/adata_liana_hvg.h5mu"

# %% Read mdata
print("INFO: Loading MuData objects")
mdata = mu.read_h5mu(h5mu_input)
mdata_liana = mu.read_h5mu(h5mu_liana)

# %% Align samples
print("INFO: Subsetting MuData objects to shared samples")
shared_samples = mdata.obs_names.intersection(mdata_liana.obs_names)
if len(shared_samples) == 0:
    raise ValueError("No shared samples found between transcriptome and LIANA MuData objects")

for view in list(mdata.mod.keys()):
    adata = mdata.mod[view]
    keep_mask = adata.obs_names.isin(shared_samples)
    mdata.mod[view] = adata[keep_mask, :].copy()
mdata.update()

for view in list(mdata_liana.mod.keys()):
    adata = mdata_liana.mod[view]
    keep_mask = adata.obs_names.isin(shared_samples)
    mdata_liana.mod[view] = adata[keep_mask, :].copy()
mdata_liana.update()

# %% Collect HVGs and filter transcriptome views
print("INFO: Filtering transcriptome views to HVGs")
hvg_genes_by_view = {}
for view in list(mdata.mod.keys()):
    adata = mdata.mod[view]
    if "highly_variable" not in adata.var:
        raise ValueError(f"Missing highly_variable column for transcriptome view {view}")

    hvg_mask = adata.var["highly_variable"].astype(bool).to_numpy()
    if hvg_mask.sum() == 0:
        raise ValueError(f"No HVGs retained for transcriptome view {view}")

    hvg_genes = {strip_view_prefix(feature) for feature in adata.var_names[hvg_mask]}
    hvg_genes_by_view[view] = hvg_genes
    print(f"INFO: {view} keeping {hvg_mask.sum()}/{adata.n_vars} HVGs")
    mdata.mod[view] = adata[:, hvg_mask].copy()
mdata.update()

# %% Filter LIANA interactions to those containing at least one sender/receiver HVG
print("INFO: Filtering LIANA interactions to LR pairs containing at least one HVG")
filtered_liana_views = {}
for view in list(mdata_liana.mod.keys()):
    if "&" not in view:
        print(f"INFO: Skipping LIANA view without sender/receiver separator: {view}")
        continue

    sender, receiver = view.split("&", 1)
    if sender not in hvg_genes_by_view or receiver not in hvg_genes_by_view:
        print(f"INFO: Skipping {view}; sender or receiver transcriptome view is absent")
        continue

    adata = mdata_liana.mod[view]
    sender_hvgs = hvg_genes_by_view[sender]
    receiver_hvgs = hvg_genes_by_view[receiver]

    keep_mask = []
    for feature in adata.var_names:
        feature_name = strip_view_prefix(feature)
        if "^" not in feature_name:
            keep_mask.append(False)
            continue

        ligand, receptor = feature_name.split("^", 1)
        ligand_genes = complex_genes(ligand)
        receptor_genes = complex_genes(receptor)
        keep_mask.append(
            bool(ligand_genes & sender_hvgs) or bool(receptor_genes & receiver_hvgs)
        )

    kept = sum(keep_mask)
    print(f"INFO: {view} keeping {kept}/{adata.n_vars} LR interactions")
    if kept == 0 or adata.n_obs == 0:
        continue

    filtered_liana_views[view] = adata[:, keep_mask].copy()

if not filtered_liana_views:
    raise ValueError("No LIANA views retained after HVG interaction filtering")

# %% Combine MuData
print("INFO: Combining transcriptome and LIANA MuData views")
combined_views = {}
for key in filtered_liana_views:
    combined_views[key] = filtered_liana_views[key]
for key in mdata.mod.keys():
    combined_views[key] = mdata.mod[key]

mdata_combined = mu.MuData(combined_views)
mdata_combined.obs = mdata.obs.loc[mdata_combined.obs_names].copy()

print(mdata_combined)

if "snakemake" in locals():
    print(f"INFO: Writing combined MuData to {h5mu_output}")
    mdata_combined.write(h5mu_output)
