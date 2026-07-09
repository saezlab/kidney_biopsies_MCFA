rule plot_figure3a_toploadings_dotplot:
    input:
        loadings="results/MCFA/loadings~hvg~5.csv",
        centering_means="results/MCFA/centering_means~hvg~5.csv"
    output:
        pdf="plots/figures/Figure3/3a_toploadings_dotplot.pdf",
        png="plots/figures/Figure3/3a_toploadings_dotplot.png"
    params:
        factor="Factor1",
        top_n=10,
        centering_group="group1",
        views=["PT", "DCT", "PEC", "TAL", "CNT_PC", "Stromal"]
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/3a_plot_toploadings_dotplot.R"
