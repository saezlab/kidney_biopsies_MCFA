rule plot_updated_rna_umap_annotations:
    input:
        umap="results/seurat/subset_high_quality/umap_embeddings.csv",
        metadata="results/seurat/subset_high_quality/metadata.csv",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        plot="plots/figures/Figure1/rna_umap_level1_level2.pdf"
    params:
        cluster_column="SCT_snn_res.0.6"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/plot_rna_umap_level1_level2.R"

rule plot_level1_marker_dotplot:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds",
        markers_csv="results/seurat/subset_high_quality/markers/findallmarkers_level1.csv",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        plot="plots/figures/Figure1/level1_marker_dotplot.pdf"
    params:
        cluster_column="SCT_snn_res.0.6",
        group_column="level1",
        top_n=3
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/plot_level1_marker_dotplot.R"

rule plot_level2_marker_dotplot:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds",
        markers_csv="results/seurat/subset_high_quality/markers/findallmarkers_level2.csv",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        plot="plots/figures/Figure1/level2_marker_dotplot.pdf"
    params:
        cluster_column="SCT_snn_res.0.6",
        level1_column="level1",
        level2_column="level2",
        top_n=3
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/figures/plot_level2_marker_dotplot.R"
