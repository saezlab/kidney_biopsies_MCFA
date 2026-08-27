# %% Import libraries
import decoupler as dc
import pandas as pd

if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    prior = snakemake.input['prior']
    enrichment_output = snakemake.output['enrichment_output']
    factor_name = snakemake.params['factor']
else:
    loadings = 'results/MCFA/loadings~hvg~5.csv'
    prior = 'data/prior_knowledge/hallmark.csv'
    enrichment_output = 'results/MCFA/hallmark_Factor1.csv'
    factor_name = 'Factor1'

# %%
hallmark = pd.read_csv(prior)
hallmark = hallmark[~hallmark.duplicated(['source', 'target'])]

# %%
loadings_all = pd.read_csv(loadings, index_col=0)

factor_loadings = loadings_all[(loadings_all['Factor'] == factor_name) &
                               (~loadings_all['view'].str.contains('&'))]
# %%
def enrichment(celltype, net, positive = True, n = 200):
    df = factor_loadings[(factor_loadings['view'] == celltype)]
    if df.empty:
        print(f"No loadings found for {celltype}")
        return None
    if positive:
        genes = df.sort_values('Loading', ascending = False).head(n)['variable'].to_numpy()
    else:
        genes = df.sort_values('Loading', ascending = True).head(n)['variable'].to_numpy()
    genes = pd.Index(genes).drop_duplicates()
    universe = pd.Index(net['target']).drop_duplicates()
    ora_input = pd.DataFrame(-1, index=[celltype], columns=universe)
    ora_input.loc[celltype, universe.intersection(genes)] = 1
    n_selected = int((ora_input.loc[celltype] == 1).sum())
    n_up = max(len(universe) - n_selected, 1)
    scores, padj = dc.mt.ora(data=ora_input, net=net, tmin=5, n_up=n_up)
    enr = pd.concat(
        [
            scores.T.rename(columns={celltype: 'Combined score'}),
            padj.T.rename(columns={celltype: 'FDR p-value'})
        ],
        axis=1
    )
    return enr.rename_axis('Term').reset_index().sort_values('FDR p-value')
# %%
#####################
# HALLMARK
#####################
# Initialize lists to store results
enr_list = []
cell_types = factor_loadings['view'].unique()

# Run pathway_activity for each factor and cell type
results = [
    enrichment(celltype=celltype, net=hallmark, positive=True, n=200).assign(view=celltype)
    for celltype in cell_types
]

# %%
# Combine all results into a single DataFrame
enr_combined_up = pd.concat(results, axis=0).assign(direction = 'up')


##### Now for the negative loadings
# %%
# Initialize lists to store results
enr_list_down = []

cell_types = factor_loadings['view'].unique()

# Run pathway_activity for each factor and cell type
results = [
    enrichment(celltype=celltype, net=hallmark, positive=False, n=200).assign(view=celltype)
    for celltype in cell_types
]

# %%
# Combine all results into a single DataFrame
enr_combined_down = pd.concat(results, axis=0).assign(direction = 'down')

# %%
# Combine up and down results
enr_combined = pd.concat([enr_combined_up, enr_combined_down], axis=0)

# %%
enr_combined.to_csv(enrichment_output)
# %%
