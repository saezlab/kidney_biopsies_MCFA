# 31.10.2023 Charlotte Boys
# Updated 26.06.2026 - Charlotte Boys

# %% Import libraries
import platform as platform
import os as os
import yaml as yaml
import pandas as pd
import liana as li
import muon as mu
import mofax as mofa


# %% Parse snakemake object
if 'snakemake' in locals():
    h5mu_input = snakemake.input['h5mu_input']
    model_output = snakemake.output['model_output']
    scores_output = snakemake.output['scores_output']
    scores_wide_output = snakemake.output['scores_wide_output']
    loadings_output = snakemake.output['loadings_output']
    nfactors = snakemake.params['nfactors']
    hvg = snakemake.params['hvg']
    sample_key = snakemake.params['sample_key']
else:
    h5mu_input = 'results/MCFA/adata.h5mu'
    model_output = 'results/MCFA/model~all~5.h5ad'
    scores_output = 'results/MCFA/scores~all~5.csv'
    scores_wide_output = 'results/MCFA/widescores~all~5.csv'
    loadings_output = 'results/MCFA/loadings~all~5.csv'
    nfactors = 5
    hvg = 'all'
    sample_key = 'sample'

# %% Read mdata
print('INFO: Loading assay')
mdata = mu.read_h5mu(h5mu_input)

# %% Set up to run with HVGs
if hvg == 'all':
    use_var = None
elif hvg == 'hvg':
    use_var = None
    # Physically filter to HVGs. Direct AnnData slicing keeps layers aligned
    # with X/var, unlike mu.pp.filter_var() in this container stack.
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
    
# %% Run mofa
print('INFO: Running MOFA')
mu.tl.mofa(mdata,
           groups_label = None,
           use_obs='union',
           convergence_mode='slow',
           n_factors=nfactors,
           seed=1337,
           outfile=model_output,
           use_var=use_var
           )
print(mdata)
print(mdata.varm)
# %% Extract factor scores and loadings
print('INFO: Extracting MOFA factor scores and loadings')
factor_scores_obs = li.utils.get_factor_scores(mdata, obsm_key='X_mofa')

print(factor_scores_obs.head())
factor_scores_wide = factor_scores_obs.set_index('index').filter(like='Factor', axis=1)
factor_scores_long = factor_scores_wide.rename_axis('Sample').reset_index().melt(id_vars=['Sample'],
                                                           var_name='Factor',
                                                           value_name='Score')

variable_loadings =  li.utils.get_variable_loadings(mdata, varm_key = 'LFs', view_sep=':')
print(variable_loadings.head())
variable_loadings_long = pd.melt(variable_loadings, 
                                 id_vars=['view', 'variable'],
                                 var_name='Factor',
                                 value_name='Loading')
# %% Save factor scores and loadings
print('INFO: Saving MOFA factor scores and loadings')
factor_scores_long.to_csv(scores_output)
factor_scores_wide.to_csv(scores_wide_output)
variable_loadings_long.to_csv(loadings_output)
# %%
