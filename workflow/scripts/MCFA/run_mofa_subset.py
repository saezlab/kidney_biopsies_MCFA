# 31.10.2023 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import yaml as yaml
import pandas as pd
import liana as li
import muon as mu
import mofax as mofa


# %% Parse snakemake object TODO: update for new repo
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
    h5mu_input = 'results/MCFA/adata_combined.h5mu'
    model_output = 'results/MCFA/model~all~5.h5ad'
    scores_output = 'results/MCFA/scores~all~5.csv'
    scores_wide_output = 'results/MCFA/widescores~all~5.csv'
    loadings_output = 'results/MCFA/loadings~all~5.csv'
    nfactors = 5
    hvg = 'all'
    sample_key = 'Sample'

# %% Read mdata
print('INFO: Loading assay')
mdata = mu.read_h5mu(h5mu_input)

# %% Set up to run with HVGs
# First of all, if we're running with LIANA, the hvg selection has already been done or not depending
# on hvg parameter in rule combine_liana_mdata
if any("&" in key for key in mdata.mod.keys()):
    use_var = None
elif hvg == 'all':
    use_var = None
elif hvg == 'hvg':
    use_var = 'highly_variable' # I think there's a bug in muon that this doesn't actually make a difference
    # So we filter to only keep HVGs:
    for view in mdata.mod.keys():
        mu.pp.filter_var(mdata.mod[view], var = 'highly_variable')
    
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
