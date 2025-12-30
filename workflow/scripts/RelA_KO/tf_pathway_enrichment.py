# %% Import libraries
import platform as platform
import os as os
import decoupler as dc
import pandas as pd
import numpy as np

if 'snakemake' in locals():
    loadings = snakemake.input['loadings']
    grn = snakemake.input['grn']
else:
    import os
    os.chdir('/Users/charlotteboys/GitHub/KuppeBoys2024')
    degs_path = 'data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt'
    grn = 'data/ATAC/grn~PT.csv'

# %%
celloracle = pd.read_csv(grn)
# %%
collectri = dc.op.collectri(organism='human')
# %%
progeny = dc.op.progeny(organism='human')
progeny_500 = dc.op.progeny(organism='human', top = 500)
progeny_200 = dc.op.progeny(organism='human', top = 200)
progeny_100 = dc.op.progeny(organism='human', top = 100)
# %%
degs = pd.read_table(degs_path, index_col=0)

mat = degs[['stat']].T
# %%
progeny_acts, progeny_padj = dc.mt.ulm(data=mat, net=progeny, tmin=5)
# %%
celloracle_acts, celloracle_padj = dc.mt.ulm(data=mat, net=celloracle, tmin=5)
# %%
collectri_acts, collectri_padj = dc.mt.ulm(data=mat, net=collectri, tmin=5)

# %%
df3 = pd.concat(
    [progeny_padj.T.rename(columns={'stat': 'padj'}),
    progeny_acts.T.rename(columns={'stat': 'act'})],
    axis=1)
df3.index.rename('gene')
df3.rename_axis('Pathway').reset_index().to_csv('results/RelA_KO/progeny_pathways.csv', index=False)
# %%
df = pd.concat(
    [collectri_padj.T.rename(columns={'stat': 'padj'}),
    collectri_acts.T.rename(columns={'stat': 'act'})],
    axis=1)
df.index.rename('gene')
# %%
df.rename_axis('TF').reset_index().to_csv('results/RelA_KO/collectri_TFs.csv', index=False)

# %%
df2 = pd.concat(
    [celloracle_padj.T.rename(columns={'stat': 'padj'}),
    celloracle_acts.T.rename(columns={'stat': 'act'})],
    axis=1)
df2.index.rename('gene')
df2.rename_axis('TF').reset_index().to_csv('results/RelA_KO/celloracle_TFs.csv', index=False)
# %%
dc.pl.source_targets(
    data=degs,
    x='weight',
    y='stat',
    net=celloracle,
    name='HNF4A',
    top=15,
    figsize=(4, 4)
)
#
#  %%
dc.pl.source_targets(
    data=degs,
    x='weight',
    y='stat',
    net=celloracle,
    name='KLF6',
    top=15,
    figsize=(4, 4)
)

# %%
dc.pl.source_targets(
    data=degs,
    x='weight',
    y='stat',
    net=celloracle,
    name='ELF3',
    top=15,
    figsize=(4, 4)
)
# %%
