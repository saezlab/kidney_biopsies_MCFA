
rule factor_enrichment:
    input:
        loadings = 'results/MCFA/loadings~all~5.csv'
    output:
        enrichment_output = 'results/downstream/hallmark_Factor2.csv'
    params:
        gene_sets = config["gene_sets"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/downstream/factor_enrichment.py'
