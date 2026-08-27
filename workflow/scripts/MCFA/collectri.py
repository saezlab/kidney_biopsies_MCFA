# %% 
import decoupler as dc
import pandas as pd
import scanpy as sc
# %%
if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    prior = snakemake.input['prior']
    acts = snakemake.output['acts']
    factor_name = snakemake.params['factor']
else:
    loadings = 'results/MCFA/loadings~hvg~5.csv'
    prior = 'data/prior_knowledge/collectri.csv'
    acts = 'results/MCFA/acts_collectri.csv'
    factor_name = 'Factor1'

# %%
collectri = pd.read_csv(prior)
# %% Collect factor loadings
loadings_all = pd.read_csv(loadings, index_col=0)

factor_loadings = loadings_all[(loadings_all['Factor'] == factor_name) &
                               (~loadings_all['view'].str.contains('&'))]

celltypes = factor_loadings['view'].unique()

# %%
def enrichment(celltype, net):
    mat = factor_loadings[factor_loadings['view'] == celltype].set_index('variable').rename(columns={'Loading': factor_name})[[factor_name]].T
    tf_acts, tf_padj = dc.mt.ulm(data=mat, net=net, tmin=5)
    df = pd.concat(
        [tf_padj.T.rename(columns={factor_name: 'padj'}),
        tf_acts.T.rename(columns={factor_name: 'act'})],
        axis=1)
    df.index.rename('gene')
    df['Factor'] = factor_name
    df['view'] = celltype
    return df

# %% Run through each cell type
res = []
for cell in celltypes:
    res.append(enrichment(cell, collectri))
df  = pd.concat(res, axis=0)
# %%
df.rename_axis('TF').reset_index().to_csv(acts, index=False)
