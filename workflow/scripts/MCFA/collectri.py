# %% 
import decoupler as dc
import pandas as pd
import scanpy as sc
# %%
if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    acts = snakemake.output['acts']
else:
    loadings = 'results/MCFA/loadings~all~5.csv'
    acts = 'results/MCFA/acts_collectri.csv'

# %%
collectri = dc.op.collectri(organism='human')
# # %% Save collectri for reference
# collectri.to_csv('data/PriorKnowledge/collectri.csv')
# %% Collect Factor 2 loadings
loadings_all = pd.read_csv(loadings, index_col=0)

loadings_f2 = loadings_all[(loadings_all['Factor'] == "Factor2") &
                           (~loadings_all['view'].str.contains('&'))]

celltypes = loadings_f2['view'].unique()

# %%
def enrichment(celltype, net):
    mat = loadings_f2[loadings_f2['view'] == celltype].set_index('variable').rename(columns={'Loading':'Factor2'})[['Factor2']].T
    tf_acts, tf_padj = dc.mt.ulm(data=mat, net=net, tmin=5)
    df = pd.concat(
        [tf_padj.T.rename(columns={'Factor2': 'padj'}),
        tf_acts.T.rename(columns={'Factor2': 'act'})],
        axis=1)
    df.index.rename('gene')
    df['Factor'] = 'Factor2'
    df['view'] = celltype
    return df

# %% Run through each cell type
res = []
for cell in celltypes:
    res.append(enrichment(cell, collectri))
df  = pd.concat(res, axis=0)
# %%
df.rename_axis('TF').reset_index().to_csv(acts, index=False)