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
    adata_path = 'kidney_biopsies_MCFA/data/adata20240621.h5ad'
    annotation ='main_class'# 'subcluster_class' #  
    subset = 'atlas' # 'FIB_VSMC'
    clr_csv = 'kidney_biopsies_MCFA/results/composition/clr~' + annotation + '~' + subset + '.csv'
    prop_csv = 'kidney_biopsies_MCFA/results/composition/prop~' + annotation + '~' + subset + '.csv'

    if platform.system() == 'Linux':
        adata_path = os.path.join('/mnt/sds-hd/sd22b002/projects', adata_path)
        clr_csv = os.path.join('/mnt/sds-hd/sd22b002/projects', clr_csv)
        prop_csv = os.path.join('/mnt/sds-hd/sd22b002/projects', prop_csv)

    elif platform.system() == 'Darwin':
        adata_path = os.path.join('/Users/charlotteboys/GitHub', adata_path)
        clr_csv = os.path.join('/Users/charlotteboys/GitHub', clr_csv)
        prop_csv = os.path.join('/Users/charlotteboys/GitHub', prop_csv)

# %%
adata = sc.read_h5ad(adata_path)

# %% Subset if necessary
if subset in adata.obs['main_class'].unique():
    print("Subsetting Atlas")
    adata = adata[adata.obs['main_class'] == subset]
else:
    print("No subsetting performed") 

# %% Regularise annotations
# Assuming adata is your AnnData object and 'annotation' is the key for the annotation column
adata.obs[annotation] = adata.obs[annotation].str.replace('+', 'pos')
adata.obs[annotation] = adata.obs[annotation].str.replace('-', '_')
adata.obs[annotation] = adata.obs[annotation].str.replace('/', '_')
adata.obs[annotation] = adata.obs[annotation].str.replace(' ', '')

# %%
# Group the DataFrame by 'Sample' and count the occurrences of each 'main_class'
grouped = adata.obs.groupby('Sample')[annotation].value_counts().unstack(fill_value=0)
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
