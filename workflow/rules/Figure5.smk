# Figure 5a-c: Factor 1 Hallmark, PROGENy, and CollecTRI enrichments
rule plot_factor_enrichments:
    input:
        collectri="results/MCFA/acts_collectri.csv",
        progeny="results/MCFA/acts_progeny.csv",
        hallmark="results/MCFA/hallmark_Factor1.csv"
    output:
        pdf="plots/figures/Figure5/factor_enrichments.pdf",
        png="plots/figures/Figure5/factor_enrichments.png"
    params:
        factor="Factor1",
        top_hallmark_sets=10,
        top_collectri_tfs=20,
        views=["PT", "DCT", "PEC", "TAL", "CNT_PC", "Stromal"]
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/factor_enrichments.R"


# Figure 5g: RelA KO Factor 1 enrichment score
rule plot_rela_factor_score:
    input:
        enrichment="results/RelA_KO/factor1_enrichment.csv"
    output:
        pdf="plots/figures/Figure5/rela_factor_score.pdf",
        png="plots/figures/Figure5/rela_factor_score.png"
    params:
        y_label="Factor 1 PT Enrichment Score"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/rela_factor_score.R"


# Figure 5i: RelA KO PROGENy pathway activities
rule plot_rela_pathways:
    input:
        pathways="results/RelA_KO/progeny_pathways.csv"
    output:
        pdf="plots/figures/Figure5/rela_pathways.pdf",
        png="plots/figures/Figure5/rela_pathways.png"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/rela_pathways.R"


# Figure 5j: Factor 1 versus RelA KO CollecTRI TF activities
rule plot_rela_tfs:
    input:
        rela_tfs="results/RelA_KO/collectri_TFs.csv",
        factor_tfs="results/MCFA/acts_collectri.csv"
    output:
        pdf="plots/figures/Figure5/rela_tfs.pdf",
        png="plots/figures/Figure5/rela_tfs.png"
    params:
        view="PT",
        factor_label="Factor 1"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/rela_tfs.R"
