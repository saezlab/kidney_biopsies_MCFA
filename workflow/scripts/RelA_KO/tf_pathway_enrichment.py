import pandas as pd
import decoupler as dc


if 'snakemake' in locals():
    degs_path = snakemake.input['degs']
    progeny_prior = snakemake.input['progeny_prior']
    collectri_prior = snakemake.input['collectri_prior']
    progeny_out = snakemake.output['progeny']
    collectri_out = snakemake.output['collectri']
else:
    degs_path = 'data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt'
    progeny_prior = 'data/prior_knowledge/progeny.csv'
    collectri_prior = 'data/prior_knowledge/collectri.csv'
    progeny_out = 'results/RelA_KO/progeny_pathways.csv'
    collectri_out = 'results/RelA_KO/collectri_TFs.csv'


# Read DEGs (expect index column = gene and a 'stat' column)
degs = pd.read_table(degs_path, index_col=0)
if 'stat' not in degs.columns:
    raise RuntimeError("Input DEGs table must contain a 'stat' column (DESeq2 result).")

mat = degs[['stat']].T

# PROGENy: pathway activities
progeny = pd.read_csv(progeny_prior)
progeny_acts, progeny_padj = dc.mt.ulm(data=mat, net=progeny, tmin=5)

df_progeny = pd.concat([
    progeny_padj.T.rename(columns={'stat': 'padj'}),
    progeny_acts.T.rename(columns={'stat': 'act'})
], axis=1)
df_progeny = df_progeny.rename_axis('Pathway').reset_index()
df_progeny.to_csv(progeny_out, index=False)


# Collectri: TF activities
collectri = pd.read_csv(collectri_prior)
collectri_acts, collectri_padj = dc.mt.ulm(data=mat, net=collectri, tmin=5)

df_collectri = pd.concat([
    collectri_padj.T.rename(columns={'stat': 'padj'}),
    collectri_acts.T.rename(columns={'stat': 'act'})
], axis=1)
df_collectri = df_collectri.rename_axis('TF').reset_index()
df_collectri.to_csv(collectri_out, index=False)
