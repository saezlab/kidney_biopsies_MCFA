# 26.04.2023 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml
import liana as li

from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object 
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    h5mu_output = snakemake.output['h5mu_output']
    sample_key = snakemake.params['sample_key']
    groupby = snakemake.params['groupby']
else:
    sample_key = 'Sample'
    groupby = 'main_class'
    h5ad_input = 'results/mofaxliana/adata_LIANA.h5ad'

# %% Read in AnnData object
print('INFO: Loading data')
adata = sc.read_h5ad(h5ad_input)

# %% 
print('INFO: Creating multi-view structure')
mdata = li.multi.lrs_to_views(adata,
                              sample_key = sample_key,
                              score_key='magnitude_rank',
                              obs_keys=['disease', 'sex', 'IFTA',
                                  'corticosteroids',
                                  'classification', 'jm_slope', 'immediate_slope'], # add those to mdata.obs
                              lr_prop = 0.1, # minimum required proportion of samples to keep an LR
                              lrs_per_sample = 20, # minimum number of interactions to keep a sample in a specific view
                              lrs_per_view = 20, # minimum number of interactions to keep a view
                              samples_per_view = 10, # minimum number of samples to keep a view
                              min_variance = 0, # minimum variance to keep an interaction
                              lr_fill = 0, # fill missing LR values across samples with this
                              verbose=True
                              )


mdata

if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5mu_output))
    mdata.write(h5mu_output)