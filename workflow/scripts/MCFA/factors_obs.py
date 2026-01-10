# %% Import libraries
import liana as li
import decoupler as dc
import pandas as pd
import platform
import os
import scanpy as sc

# %% Parse snakemake object
if 'snakemake' in locals():
    scores_wide = snakemake.input['scores_wide']
    metadata_csv = snakemake.input['metadata_csv']
    clr = snakemake.input['clr']
    factors_obs_csv = snakemake.output['factors_obs_csv']
else:
    import os
    scores_wide = 'results/MCFA/widescores~all~5.csv'
    metadata_csv = 'data/metadata/patient_metadata_EDA.csv'
    clr = 'results/composition/clr~main_class.csv'
    factors_obs_csv = 'results/MCFA/factors_obs~all~5.csv'
    
# %% Read in factor scores
factor_scores = pd.read_csv(scores_wide, index_col=0).rename_axis('Sample').filter(like='Factor', axis=1)

# %% Read in proportions
clr_props = pd.read_csv(clr, index_col=0).reindex(factor_scores.index)
clr_props.columns = [col + '_clr' for col in clr_props.columns]

# %% Read in metadata
metadata = pd.read_csv(metadata_csv, index_col=0).reindex(factor_scores.index)

# %% Save data for plotting
obs = pd.concat([metadata, clr_props], axis='columns')
to_plot = pd.concat([factor_scores, obs], axis = 1)
to_plot.to_csv(factors_obs_csv)