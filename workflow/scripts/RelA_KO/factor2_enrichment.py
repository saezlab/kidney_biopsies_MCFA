import os
import pandas as pd
import decoupler as dc

# Use snakemake inputs when available, otherwise fall back to repo defaults
if 'snakemake' in locals():
    tmm_path = snakemake.input['tmm']
    loadings_path = snakemake.input['loadings']
    out_path = snakemake.output['factor2_enrichment']
else:
    tmm_path = 'data/RelA_KO/post_TMM.txt'
    loadings_path = 'results/MCFA/loadings~all~5.csv'
    out_path = 'results/RelA_KO/factor2_enrichment.csv'

# Read TMM
tmm = pd.read_csv(tmm_path, sep='\t', index_col=0)

# Subset expected columns if present (EC_AA and EC_AB replicates)
expected_cols = ['EC_AA_REP1', 'EC_AA_REP2', 'EC_AA_REP3', 'EC_AB_REP1', 'EC_AB_REP2', 'EC_AB_REP3']
cols = [c for c in expected_cols if c in tmm.columns]
if not cols:
    raise RuntimeError(f'None of the expected replicate columns found in {tmm_path}')
tmm = tmm[cols]

# Drop rows where all values are zero
tmm = tmm[(tmm != 0).any(axis=1)]
# %% Z score transform gene expression across all replicates
tmm_zcore = tmm.subtract(tmm.mean(axis=1), axis=0).divide(tmm.std(axis=1), axis=0)
# %% Load factor 2 pt loadings
loadings = pd.read_csv(loadings_path, index_col=0)
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
df.to_csv(out_path, index=False)