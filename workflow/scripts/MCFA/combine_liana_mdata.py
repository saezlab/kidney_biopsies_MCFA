# 26.04.2024 Charlotte Boys
# %% Import libraries
import platform as platform
import os as os
import yaml as yaml
import muon as mu


# %% Parse snakemake object
if 'snakemake' in locals():
    h5mu_input = snakemake.input['h5mu_input']
    h5mu_liana = snakemake.input['h5mu_liana']
    h5mu_output = snakemake.output['h5mu_output']
    hvg = snakemake.params['hvg']
else:
    h5mu_input = 'results/MCFA/adata.h5mu'
    h5mu_liana = 'results/MCFA/adata_LIANA.h5mu'
    hvg = 'hvg' #'all'

# %% Read mdata
print('INFO: Loading assay')
mdata = mu.read_h5mu(h5mu_input)
mdata_liana = mu.read_h5mu(h5mu_liana)

# %%
print('INFO: Subsetting mudata objects')
# Make sure we're looking at the same samples
mu.pp.filter_obs(mdata, var = mdata_liana.obs_names)
mu.pp.filter_obs(mdata_liana, var = mdata.obs_names)

# %% Set up to run with HVGs
if hvg == 'hvg':
    for view in mdata.mod.keys():
        mu.pp.filter_var(mdata.mod[view], var = 'highly_variable')

# %% Filter interactions for those where one or both of the ligand and receptor are HVGs in their respective views
if hvg == 'hvg':
    for view in mdata_liana.mod.keys():
        sender = view.split('&')[0]
        receiver = view.split('&')[1]
        hvg_ligands = [i[1] for i in mdata.mod[sender].var_names.str.split(':')]
        hvg_receptors = [i[1] for i in mdata.mod[receiver].var_names.str.split(':')]
        interactions = [i[1] for i in mdata_liana.mod[view].var_names.str.split(':')]
        mask = [
                    interaction.split('^')[0] in hvg_ligands or interaction.split('^')[1] in hvg_receptors
                    for interaction in interactions
                ]
        mdata_liana.mod[view].var['highly_variable'] = mask
        mu.pp.filter_var(mdata_liana.mod[view], var = 'highly_variable')
# %%
print('INFO: Combining mudata')
# Initialize an empty dictionary to store AnnData objects
adict = {}

# Iterate through the keys of mdata_liana.mod
for key in mdata_liana.mod.keys():
    # Retrieve the AnnData object corresponding to the key
    adict[key] = mdata_liana.mod[key]
for key in mdata.mod.keys():
    adict[key] = mdata.mod[key]

mdata_combined = mu.MuData(adict)
mdata_combined.obs = mdata.obs

mdata_combined
# %%
if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5mu_output))
    mdata_combined.write(h5mu_output)