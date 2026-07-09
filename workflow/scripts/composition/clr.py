# %% 
import pandas as pd
import scanpy as sc
import os
import platform
from skbio.stats.composition import clr, multiplicative_replacement

# %% Parse snakemake
if 'snakemake' in locals():
    adata_path = snakemake.input['adata']
    clr_csv = snakemake.output['clr_csv']
    prop_csv = snakemake.output['prop_csv']
    annotation = snakemake.wildcards['annotation']
    subset = snakemake.wildcards['subset']
else:
    adata_path = 'results/seurat/subset_high_quality/adata_20260717.h5ad'
    annotation ='level1'  
    subset = 'atlas' 
    clr_csv = 'results/composition/clr~' + annotation + '~' + subset + '.csv'
    prop_csv = 'results/composition/prop~' + annotation + '~' + subset + '.csv'

# %%
adata = sc.read_h5ad(adata_path)

# %% Subset if necessary
if subset in adata.obs['level1'].unique():
    print("Subsetting Atlas")
    adata = adata[adata.obs['level1'] == subset]
else:
    print("No subsetting performed") 

# %% Regularise annotations
# Assuming adata is your AnnData object and 'annotation' is the key for the annotation column
adata.obs[annotation] = adata.obs[annotation].str.replace('+', 'pos')
adata.obs[annotation] = adata.obs[annotation].str.replace('-', '_')
adata.obs[annotation] = adata.obs[annotation].str.replace('/', '_')
adata.obs[annotation] = adata.obs[annotation].str.replace(' ', '')

# %%
# Group the DataFrame by sample and count the occurrences of each 'level1' annotation within each sample
grouped = adata.obs.groupby('sample')[annotation].value_counts().unstack(fill_value=0)
# Remove rows which are all zero (K68)
grouped = grouped.loc[(grouped != 0).any(axis=1)]
# Calculate the proportions within each group
proportions = grouped.div(grouped.sum(axis=1), axis=0).values

# %%
# Get row and column names
sample_names = grouped.index
annotation_names = grouped.columns

# Create a new DataFrame with proportions and appropriate row and column names
proportions_df = pd.DataFrame(proportions, index=sample_names, columns=annotation_names)

# %% Replace zeros with small positive delta making sure all proportions still add to 1
# Then perform clr transform
clr_mat = clr(multiplicative_replacement(proportions))
clr_df = pd.DataFrame(clr_mat, index=sample_names, columns=annotation_names)

# %%
clr_df.to_csv(clr_csv)
proportions_df.to_csv(prop_csv)

# %%
