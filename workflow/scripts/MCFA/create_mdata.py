# 31.10.2023 - Charlotte Boys
# Updated 26.06.2026 - Charlotte Boys

# %% Import libraries
import scanpy as sc
import platform as platform
import os as os
import yaml as yaml
import numpy as np
import pandas as pd
import muon as mu
import mofax as mofa
from helper_functions import adata_to_views, preprocess_data

# %% Parse snakemake object
if 'snakemake' in locals():
    h5ad_input = snakemake.input['h5ad_input']
    coldata_path = snakemake.input['coldata_path']
    h5ad_output = snakemake.output['h5ad_output']
    h5mu_output = snakemake.output['h5mu_output']
    sample_key = snakemake.params['sample_key']
    condition_key = snakemake.params['condition_key']
    groupby = snakemake.params['groupby']
    min_cells = snakemake.params['min_cells']
    sample_subset = snakemake.params['sample_subset']
else:
    h5ad_input = 'results/seurat/subset_high_quality/adata_20260717.h5ad'
    sample_key = 'sample'
    condition_key = 'disease'
    groupby = 'level1'
    coldata_path = 'data/metadata/patient_metadata_EDA.csv'
    min_cells = 25
    sample_subset = ['K1', 'K100', 'K101',
       'K11', 'K12', 'K14', 'K16', 'K17', 'K19', 'K2',
       'K25', 'K26', 'K29', 'K3', 'K30', 'K31', 'K32', 'K33', 'K34', 'K35',
       'K36', 'K37', 'K38', 'K40', 'K41', 'K42', 'K43', 'K46', 'K48',
       'K49', 'K5', 'K50', 'K52', 'K53', 'K54', 'K55', 'K56', 'K58', 'K59',
       'K6', 'K61', 'K62', 'K67', 'K68', 'K70', 'K71',
       'K73', 'K74', 'K75', 'K76', 'K77', 'K78', 'K79',
       'K82', 'K84', 'K89', 'K92', 'K93', 'K94', 'K95', 'K97',
       'K98', 'K99']
# %% Read in AnnData object
adata = sc.read_h5ad(h5ad_input)

# Remove un-annotated cells
adata = adata[~adata.obs[groupby].astype(str).isin(['nan'])]

# %%
# Keep only samples in subset
if sample_subset != 'all':
    adata = adata[adata.obs[sample_key].isin(sample_subset)]

# %% Add coldata to adata.obs
coldata = pd.read_csv(coldata_path)
coldata = coldata.rename(columns= {"sample_ID": sample_key})
coldata_rna = coldata[coldata[sample_key].isin(adata.obs[sample_key])]

# %%
print("INFO: Merging metadata by sample key")
print(sample_key)
obs_id_col = "__cellid__"
obs = adata.obs.copy()
obs[obs_id_col] = obs.index
new_obs = obs.merge(coldata_rna, how = "left", on = [sample_key])
# %%
new_obs = new_obs.drop_duplicates(subset=obs_id_col).set_index(obs_id_col)
new_obs.index.name = adata.obs.index.name
# for all columns of type object, convert to string
object_columns = new_obs.select_dtypes(include=['object']).columns
new_obs[object_columns] = new_obs[object_columns].astype(str)
# for all columns of type float64, convert to float32 
float_columns = new_obs.select_dtypes(include=['float64']).columns
new_obs[float_columns] = new_obs[float_columns].astype('float32')

# %%
adata.obs = new_obs.reindex(adata.obs.index)
adata.obs = adata.obs.rename(columns = {'IFTA_y': 'IFTA', 'sex_y': 'sex', 'date_of_biopsy_y': 'date_of_biopsy', 'corticosteroids_y': 'corticosteroids'})


print("INFO: adata.obs columns")
print(list(adata.obs))
# %% 
# Require minimum genes per cell and minimum cells per gene
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# %%
# Crucial for being able to save and read the mudata object later on
# Otherwise you'll get an error message like TypeError: AnnData.__init__() got an unexpected keyword argument 'TAL'
# If one of your cell type views is called EC/TAL
# Or TypeError: AnnData.__init__() got an unexpected keyword argument 'P'
# if one of your cell type  views is called VSM/P
# First save the original annotations
adata.obs['original_labels'] = adata.obs[groupby]

adata.obs[groupby] = adata.obs[groupby].str.replace('/', '_')
adata.obs[groupby] = adata.obs[groupby].str.replace(' ', '')

adata.X = np.round(adata.X)
adata.layers['counts'] = adata.X
# %% Create a multi-view structure
mdata_original = adata_to_views(adata,
                        groupby=groupby,
                        sample_key=sample_key,
                        obs_keys=['disease', 'sex', 'IFTA',
                                  'corticosteroids',
                                  'classification', 'jm_slope', 'immediate_slope'], # add those to mdata.obs
                        min_prop=0.05, # min nnz values (filter features)
                        min_smpls=3, # min samples per view (filter features)
                        min_cells=min_cells, # min cells per view (filter samples)
                        min_counts=100, # min counts per view (filter samples)
                        mode='sum', # mode of aggregation
                        verbose=True,
                        large_n=5, # edgeR-like filtering
                        min_total_count=15,
                        min_count=10,
                        keep_psbulk_stats = True,
                        skip_checks = True # Ignore that counts are not integers (because of SoupX correction)
                        )

# Enforce minimum number of views per sample:
# lambda function to check if sample is in a view
in_view = lambda view: mdata_original.obs_names.isin([*mdata_original.mod[view].obs_names])*1
# Check for two or more views per sample
mask = sum(list(map(in_view, [*mdata_original.mod.keys()]))) >= 2
pass_samples = mdata_original.obs_names[mask]
for current_mod in mdata_original.mod.keys():
    s = mdata_original.mod[current_mod].obs_names.isin(pass_samples)
    mdata_original.mod[current_mod] = mdata_original.mod[current_mod][s, :].copy()
mdata_original.update()
mdata_original.uns['psbulk_stats'] = mdata_original.uns['psbulk_stats'][mask]

print("INFO: debugging mdata")
print(list(mdata_original.obs))
# print("INFO: debugging views")
# print(list(mdata_original.mod['PT'].obs))

# %% Pre-process the pseudobulk profiles
mdata = preprocess_data(mdata_original,
                        target_sum = 10000,
                        remove_mitochondrial = False,
                        exclude_highly_expressed = False) # default number of top genes

if 'snakemake' in locals():
    print('INFO: Writing AnnData to {0}'.format(h5ad_output))
    adata.write(h5ad_output)
    print('INFO: Writing MuData to {0}'.format(h5mu_output))
    mdata.write(h5mu_output)
# %%
