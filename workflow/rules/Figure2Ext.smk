rule plot_ext2a_model_explanatory_power:
    input:
        r2_by_sample="results/MCFA/robustness_r2_by_sample~hvg.csv"
    output:
        pdf="plots/figures/Figure2Ext/ext2a_model_explanatory_power.pdf",
        png="plots/figures/Figure2Ext/ext2a_model_explanatory_power.png",
        csv="results/MCFA/model_explanatory_power~hvg.csv"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/ext2a_model_explanatory_power.R"

rule plot_ext2_loadings_sensitivity_correlations:
    input:
        standard_loadings="results/MCFA/loadings~hvg~5.csv",
        exclude_tin_loadings="results/MCFA/loadings~hvg~5~exclude-TIN.csv",
        exclude_controls_loadings="results/MCFA/loadings~hvg~5~exclude-controls.csv"
    output:
        pdf="plots/figures/Figure2Ext/ext2_loadings_sensitivity_correlations.pdf",
        png="plots/figures/Figure2Ext/ext2_loadings_sensitivity_correlations.png",
        correlations_csv="results/MCFA/loadings_sensitivity_correlations~hvg~5.csv"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/ext2_loadings_sensitivity_correlations.R"

rule plot_ext2_liana_loadings_correlation:
    input:
        standard_loadings="results/MCFA/loadings~hvg~5.csv",
        liana_loadings="results/MCFA/loadings~liana_hvg~5.csv"
    output:
        pdf="plots/figures/Figure2Ext/ext2_liana_loadings_correlation.pdf",
        png="plots/figures/Figure2Ext/ext2_liana_loadings_correlation.png",
        correlations_csv="results/MCFA/loadings_liana_correlations.csv"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/ext2_liana_loadings_correlation.R"
