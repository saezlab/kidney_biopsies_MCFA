rule tf_pathway_enrichment:
    """Run PROGENy and Collectri TF enrichment from DESeq2 results."""
    input:
        degs="data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt",
        progeny_prior=ancient("data/prior_knowledge/progeny.csv"),
        collectri_prior=ancient("data/prior_knowledge/collectri.csv")
    output:
        progeny="results/RelA_KO/progeny_pathways.csv",
        collectri="results/RelA_KO/collectri_TFs.csv"
    container:
        config["containers"]["decouplerv2"]
    script:
        "workflow/scripts/RelA_KO/tf_pathway_enrichment.py"


rule factor1_enrichment:
    """Run Factor1 TF enrichment and write CSV."""
    input:
        tmm="data/RelA_KO/post_TMM.txt",
        loadings="results/MCFA/loadings~hvg~5.csv"
    output:
        factor1_enrichment = "results/RelA_KO/factor1_enrichment.csv"
    container:
        config["containers"]["decouplerv2"]
    script:
        "workflow/scripts/RelA_KO/factor1_enrichment.py"
