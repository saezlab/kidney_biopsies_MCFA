# %% Import libraries
import platform as platform
import os as os
import decoupler as dc
import pandas as pd
import numpy as np

if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    enrichment_output = snakemake.output['enrichment_output']
else:
    import os
    os.chdir('/Users/charlotteboys/GitHub/KuppeBoys2024')
    loadings = 'results/MCFA/loadings~all~5.csv'
    enrichment_output = 'results/downstream/hallmark_Factor2.csv'

# %%
msigdb = dc.get_resource('MSigDB')
hallmark = msigdb[msigdb['collection']=='hallmark']
# Remove duplicated entries
hallmark = hallmark[~hallmark.duplicated(['geneset', 'genesymbol'])]

# %%
reactome = msigdb[msigdb['collection']=='reactome_pathways']
reactome = reactome[~reactome.duplicated(['geneset', 'genesymbol'])]

# %%
go = msigdb[msigdb['collection']=='go_biological_process']
go = go[~go.duplicated(['geneset', 'genesymbol'])]

# %%
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
#%%
msigdb_selected = msigdb[msigdb['collection'].isin(['hallmark', 'reactome_pathways', 'go_biological_process'])]
# %%
loadings_all = pd.read_csv(loadings, index_col=0)

loadings_f2 = loadings_all[(loadings_all['Factor'] == "Factor2") &
                           (~loadings_all['view'].str.contains('&'))]
# %%
def enrichment(celltype, net, positive = True, n = 200):
    df = loadings_f2[(loadings_f2['view'] == celltype)]
    if df.empty:
        print(f"No loadings found for {celltype}")
        return None
    if positive:
        genes = df.sort_values('Loading', ascending = False).head(n)['variable'].to_numpy()
    else:
        genes = df.sort_values('Loading', ascending = True).head(n)['variable'].to_numpy()
    enr = dc.get_ora_df(
        df = genes,
        net = net,
        source = 'geneset',
        target = 'genesymbol'
    )
    return enr.sort_values('FDR p-value')
# %%
#####################
# HALLMARK
#####################
# Initialize lists to store results
enr_list = []
cell_types = loadings_f2['view'].unique()

# Run pathway_activity for each factor and cell type
results = [enrichment(celltype=celltype,
                        net = hallmark,
                        positive = True,
                        n = 200).assign(view = celltype) for celltype in cell_types]

# %%
# Combine all results into a single DataFrame
enr_combined_up = pd.concat(results, axis=0).assign(direction = 'up')


##### Now for the negative loadings
# %%
# Initialize lists to store results
enr_list_down = []

cell_types = loadings_f2['view'].unique()

# Run pathway_activity for each factor and cell type
results = [enrichment(celltype=celltype,
                        net = hallmark,
                        positive = False,
                        n = 200).assign(view = celltype) for celltype in cell_types]

# %%
# Combine all results into a single DataFrame
enr_combined_down = pd.concat(results, axis=0).assign(direction = 'down')

# %%
# Combine up and down results
enr_combined = pd.concat([enr_combined_up, enr_combined_down], axis=0)

# %%
enr_combined.to_csv(enrichment_output)
# %%
