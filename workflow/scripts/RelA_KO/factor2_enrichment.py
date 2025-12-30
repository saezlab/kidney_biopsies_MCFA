# %% 
import decoupler as dc
import pandas as pd

# %% Get working directory
import os
os.chdir('/Users/charlotteboys/GitHub/KuppeBoys2024')
# %% Read in TMM
tmm = pd.read_csv('data/RelA_KO/post_TMM.txt', sep='\t', index_col=0)
# subset EC_AA and EC_AB replicates
tmm = tmm[['EC_AA_REP1', 'EC_AA_REP2', 'EC_AA_REP3', 'EC_AB_REP1', 'EC_AB_REP2', 'EC_AB_REP3']]
# Drop rows where all values are zero
tmm = tmm[(tmm != 0).any(axis=1)]
# %% Z score transform gene expression across all replicates
tmm_zcore = tmm.subtract(tmm.mean(axis=1), axis=0).divide(tmm.std(axis=1), axis=0)
# %% Load factor 2 pt loadings
loadings = pd.read_csv('results/mofaxliana_nogroups/loadings~adata~LIANA_all~all~5.csv', index_col=0)

# %%
loadings_f2_PT = loadings[(loadings['Factor'] == "Factor2") & (loadings['view'] == "PT")]
# %%
net = loadings_f2_PT.rename(columns={'Loading':'weight', 'Factor': 'source', 'variable': 'target'})
mat = tmm_zcore.T
tf_acts, tf_padj = dc.mt.ulm(data=mat, net=net, tmin=5)
df = pd.concat(
    [tf_padj.rename(columns={'Factor2': 'padj'}),
    tf_acts.rename(columns={'Factor2': 'act'})],
    axis=1)
df.index.rename('gene')
df.to_csv('results/RelA_KO/factor2_enrichment.csv')
# %% Plot targets
dc.pl.source_targets(
    data = tmm_zcore,
    net = net,
    x = 'weight',
    y = 'EC_AA_REP1',
    name = 'Factor2')
# %%
# %% Plot targets
dc.pl.source_targets(
    data = tmm_zcore,
    net = net,
    x = 'weight',
    y = 'EC_AB_REP1',
    name = 'Factor2')
# %%
