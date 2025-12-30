# 26.04.2023 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml
import numpy as np
import pandas as pd
import muon as mu
import mofax as mofa
import liana as li
import decoupler as dc

# %% Parse snakemake object 
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    h5ad_output = snakemake.output['h5ad_output']
    sample_key = snakemake.params['sample_key']
    groupby = snakemake.params['groupby']
else:
    sample_key = 'Sample'
    groupby = 'main_class'
    h5ad_input = 'results/MCFA/adata.h5ad'

# %% Read in AnnData object
print('INFO: Loading data')
adata = sc.read_h5ad(h5ad_input)

# %% Filtering cells and genes has already been done
# log1p normalise the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# %% Run LIANA
print('INFO: Running LIANA')
li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    sample_key=sample_key, # sample key by which we which to loop
    expr_prop = 0.1,
    use_raw=False,
    n_perms=100, # reduce permutations for speed
    return_all_lrs=False, # we don't return all LR values to utilize MOFA's flexible views
    verbose=True, # use 'full' to show all information
    )

adata.uns["liana_res"].sort_values("magnitude_rank").head()

if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5ad_output))
    adata.write(h5ad_output)