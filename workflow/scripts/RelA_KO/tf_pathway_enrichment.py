import os
import pandas as pd
import decoupler as dc


if 'snakemake' in locals():
    degs_path = snakemake.input['degs']
else:
    degs_path = 'data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt'


# Read DEGs (expect index column = gene and a 'stat' column)
degs = pd.read_table(degs_path, index_col=0)
if 'stat' not in degs.columns:
    raise RuntimeError("Input DEGs table must contain a 'stat' column (DESeq2 result).")

mat = degs[['stat']].T

# PROGENy: pathway activities
progeny = dc.op.progeny(organism='human')
progeny_acts, progeny_padj = dc.mt.ulm(data=mat, net=progeny, tmin=5)

df_progeny = pd.concat([
    progeny_padj.T.rename(columns={'stat': 'padj'}),
    progeny_acts.T.rename(columns={'stat': 'act'})
], axis=1)
df_progeny = df_progeny.rename_axis('Pathway').reset_index()

os.makedirs('results/RelA_KO', exist_ok=True)
df_progeny.to_csv('results/RelA_KO/progeny_pathways.csv', index=False)


# Collectri: TF activities
collectri = dc.op.collectri(organism='human')
collectri_acts, collectri_padj = dc.mt.ulm(data=mat, net=collectri, tmin=5)

df_collectri = pd.concat([
    collectri_padj.T.rename(columns={'stat': 'padj'}),
    collectri_acts.T.rename(columns={'stat': 'act'})
], axis=1)
df_collectri = df_collectri.rename_axis('TF').reset_index()
df_collectri.to_csv('results/RelA_KO/collectri_TFs.csv', index=False)

