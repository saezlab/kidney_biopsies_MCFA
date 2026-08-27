import decoupler as dc


if "snakemake" in locals():
    collectri_out = snakemake.output["collectri"]
    progeny_out = snakemake.output["progeny"]
    hallmark_out = snakemake.output["hallmark"]
else:
    collectri_out = "data/prior_knowledge/collectri.csv"
    progeny_out = "data/prior_knowledge/progeny.csv"
    hallmark_out = "data/prior_knowledge/hallmark.csv"


dc.op.collectri(organism="human").to_csv(collectri_out, index=False)
dc.op.progeny(organism="human").to_csv(progeny_out, index=False)
dc.op.hallmark(organism="human").to_csv(hallmark_out, index=False)
