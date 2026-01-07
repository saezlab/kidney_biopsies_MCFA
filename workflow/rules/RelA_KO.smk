rule tf_pathway_enrichment:
    """Run PROGENy and Collectri TF enrichment from DESeq2 results."""
    input:
        degs="data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt",
    output:
        progeny="results/RelA_KO/progeny_pathways.csv",
        collectri="results/RelA_KO/collectri_TFs.csv"
    singularity:
        config["singularities"]["scanpy_update"]
    script:
        "workflow/scripts/RelA_KO/tf_pathway_enrichment.py"


rule factor2_enrichment:
    """Run Factor2 TF enrichment and write CSV."""
    input:
        tmm="data/RelA_KO/post_TMM.txt",
        loadings="results/MCFA/loadings~all~5.csv"
    output:
        factor2_enrichment = "results/RelA_KO/factor2_enrichment.csv"
    singularity:
        config["singularities"]["scanpy_update"]
    script:
        "workflow/scripts/RelA_KO/factor2_enrichment.py"
