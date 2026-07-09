rule soupx_to_sce:
    input:
        matrix="data/RNA/SoupX/{sample}/matrix.mtx.gz",
        features="data/RNA/SoupX/{sample}/features.tsv.gz",
        barcodes="data/RNA/SoupX/{sample}/barcodes.tsv.gz"
    output:
        rds="results/seurat/sce~{sample}.rds"
    wildcard_constraints:
        sample="|".join(config["rna_samples"])
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/soupx_to_sce.R"

rule run_scDblFinder:
    input:
        rds="results/seurat/sce~{sample}.rds"
    output:
        flag="results/seurat/flags/scDblFinder~{sample}.txt"
    wildcard_constraints:
        sample="|".join(config["rna_samples"])
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/run_scDblFinder.R"

rule qc_plots:
    input:
        rds="results/seurat/sce~{sample}.rds",
        flag="results/seurat/flags/scDblFinder~{sample}.txt"
    output:
        plot="plots/seurat/qc~{sample}.pdf"
    wildcard_constraints:
        sample="|".join(config["rna_samples"])
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/qc_plots.R"

rule merge_pass_qc_seurat:
    input:
        rds=expand("results/seurat/sce~{sample}.rds", sample=config["rna_pass_qc"]),
        flag=expand("results/seurat/flags/scDblFinder~{sample}.txt", sample=config["rna_pass_qc"])
    output:
        rds="results/seurat/merged_seurat.rds"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/merge_pass_qc_seurat.R"


rule run_sct_harmony_clustering:
    input:
        rds="results/seurat/merged_seurat.rds"
    output:
        rds="results/seurat/sct_harmony.rds",
        metadata_csv="results/seurat/metadata.csv",
        pca_embeddings_csv="results/seurat/pca_embeddings.csv",
        harmony_embeddings_csv="results/seurat/harmony_embeddings.csv",
        umap_pre_embeddings_csv="results/seurat/umap_pre_harmony_embeddings.csv",
        umap_embeddings_csv="results/seurat/umap_embeddings.csv",
        pca_sd_csv="results/seurat/pca_stdev.csv",
        clustree_plot="plots/seurat/annotation/clustree_sct_harmony.pdf",
        session_info_csv="results/seurat/session_info.csv"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/run_sct_harmony_clustering.R"

rule find_markers:
    input:
        rds="results/seurat/sct_harmony.rds"
    output:
        csv="results/seurat/markers/findallmarkers_res_{resolution}.csv"
    params:
        resolution="{resolution}"
    wildcard_constraints:
        resolution=r"\d+\.\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/find_markers.R"

rule plot_umap_resolution:
    input:
        umap="results/seurat/umap_embeddings.csv",
        metadata="results/seurat/metadata.csv"
    output:
        plot="plots/seurat/umap_res_{resolution}.pdf"
    params:
        resolution="{resolution}"
    wildcard_constraints:
        resolution=r"\d+\.\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_umap_resolution.R"

rule plot_umap_sct_resolutions:
    input:
        umap="results/seurat/umap_embeddings.csv",
        metadata="results/seurat/metadata.csv"
    output:
        plot="plots/seurat/annotation/umap_sct_resolutions.pdf"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_umap_sct_resolutions.R"

rule plot_umap_main_class_annotations:
    input:
        umap="results/seurat/umap_embeddings.csv",
        annotation_csv="config/biopsydata_withkpmp_final_annotations(in).csv"
    output:
        plot="plots/seurat/annotation/umap_main_class_annotations.pdf"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_umap_main_class_annotations.R"

rule plot_main_class_proportions_by_sct_cluster:
    input:
        metadata="results/seurat/metadata.csv",
        annotation_csv="config/biopsydata_withkpmp_final_annotations(in).csv"
    output:
        plot="plots/seurat/annotation/main_class_proportions_by_sct_res_0_2.pdf"
    params:
        cluster_column="SCT_snn_res.0.2"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_main_class_proportions_by_sct_cluster.R"

rule plot_qc_metrics_by_sct_cluster:
    input:
        metadata="results/seurat/metadata.csv"
    output:
        plot="plots/seurat/annotation/qc_metrics_by_sct_res_0_2.pdf"
    params:
        cluster_column="SCT_snn_res.0.2"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_qc_metrics_by_sct_cluster.R"

rule plot_cluster_proportions_by_metadata:
    input:
        metadata="results/seurat/metadata.csv",
        patient_metadata="data/metadata/patient_metadata_EDA.csv"
    output:
        plot="plots/seurat/exploratory/cluster_5_8_proportions_by_age_egfr.pdf"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_cluster_proportions_by_metadata.R"

rule plot_markers_dotplot:
    input:
        rds="results/seurat/sct_harmony.rds",
        csv="results/seurat/markers/findallmarkers_res_{resolution}.csv"
    output:
        plot="plots/seurat/dotplot_res_{resolution}.pdf"
    params:
        resolution="{resolution}"
    wildcard_constraints:
        resolution=r"\d+\.\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_markers_dotplot.R"

rule run_sct_harmony_clustering_subset:
    input:
        merged_seurat_rds="results/seurat/merged_seurat.rds",
        metadata_csv="results/seurat/metadata.csv"
    output:
        rds="results/seurat/subset_{subset}/sct_harmony.rds",
        metadata_csv="results/seurat/subset_{subset}/metadata.csv",
        pca_embeddings_csv="results/seurat/subset_{subset}/pca_embeddings.csv",
        harmony_embeddings_csv="results/seurat/subset_{subset}/harmony_embeddings.csv",
        umap_pre_embeddings_csv="results/seurat/subset_{subset}/umap_pre_harmony_embeddings.csv",
        umap_embeddings_csv="results/seurat/subset_{subset}/umap_embeddings.csv",
        pca_sd_csv="results/seurat/subset_{subset}/pca_stdev.csv",
        clustree_plot="plots/seurat/annotation/subset_{subset}/clustree_sct_harmony.pdf",
        session_info_csv="results/seurat/subset_{subset}/session_info.csv"
    params:
        subset_resolution=lambda wildcards: config["sct_subsets"][wildcards.subset]["resolution"],
        subset_clusters=lambda wildcards: config["sct_subsets"][wildcards.subset]["clusters"]
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/run_sct_harmony_clustering_subset.R"

rule plot_subset_umap_sct_resolutions:
    input:
        umap="results/seurat/subset_{subset}/umap_embeddings.csv",
        metadata="results/seurat/subset_{subset}/metadata.csv"
    output:
        plot="plots/seurat/subset_{subset}/umap_sct_resolutions.pdf"
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_umap_sct_resolutions.R"


rule plot_subset_umap_main_class_annotations:
    input:
        umap="results/seurat/subset_{subset}/umap_embeddings.csv",
        annotation_csv="config/biopsydata_withkpmp_final_annotations(in).csv"
    output:
        plot="plots/seurat/subset_{subset}/umap_main_class_annotations.pdf"
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_umap_main_class_annotations.R"


rule plot_subset_qc_metrics_by_sct_cluster:
    input:
        metadata="results/seurat/subset_{subset}/metadata.csv"
    output:
        plot="plots/seurat/subset_{subset}/qc_metrics_by_sct_res_{resolution}.pdf"
    params:
        cluster_column=lambda wildcards: "SCT_snn_res." + wildcards.resolution.replace("_", ".")
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)",
        resolution=r"\d+_\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_qc_metrics_by_sct_cluster.R"


rule plot_subset_main_class_proportions_by_sct_cluster:
    input:
        metadata="results/seurat/subset_{subset}/metadata.csv",
        annotation_csv="config/biopsydata_withkpmp_final_annotations(in).csv"
    output:
        plot="plots/seurat/subset_{subset}/main_class_proportions_by_sct_res_0_6.pdf"
    params:
        cluster_column="SCT_snn_res.0.6"
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_main_class_proportions_by_sct_cluster.R"

rule find_subset_markers:
    input:
        rds="results/seurat/subset_{subset}/sct_harmony.rds"
    output:
        csv="results/seurat/subset_{subset}/markers/findallmarkers_res_{resolution}.csv"
    params:
        resolution="{resolution}"
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)",
        resolution=r"\d+\.\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/find_markers.R"

rule plot_subset_markers_dotplot:
    input:
        rds="results/seurat/subset_{subset}/sct_harmony.rds",
        csv="results/seurat/subset_{subset}/markers/findallmarkers_res_{resolution}.csv"
    output:
        plot="plots/seurat/subset_{subset}/dotplot_res_{resolution}.pdf"
    params:
        resolution="{resolution}"
    wildcard_constraints:
        subset="|".join(config.get("sct_subsets", {}).keys()) if config.get("sct_subsets", {}) else "(?!)",
        resolution=r"\d+\.\d+"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/plot_markers_dotplot.R"

rule find_celltype_markers:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        csv="results/seurat/subset_high_quality/markers/findallmarkers_level1.csv"
    params:
        cluster_column="SCT_snn_res.0.6",
        group_column="level1"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/find_celltype_markers.R"

rule find_level2_markers:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        csv="results/seurat/subset_high_quality/markers/findallmarkers_level2.csv"
    params:
        cluster_column="SCT_snn_res.0.6",
        group_column="level2"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/find_celltype_markers.R"

rule find_dpt_markers_excluding_pt_apt:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds",
        metadata="results/seurat/subset_high_quality/metadata_annotations.csv"
    output:
        csv="results/seurat/subset_high_quality/markers/findmarkers_level2_dPT_excluding_PT_aPT.csv"
    params:
        group_column="level2",
        target_group="dPT",
        exclude_groups=["PT", "aPT"]
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/find_target_markers_excluding_groups.R"

rule annotate_high_quality_metadata:
    input:
        metadata="results/seurat/subset_high_quality/metadata.csv",
        annotation_csv="config/RNA_level1_annotations.csv"
    output:
        csv="results/seurat/subset_high_quality/metadata_annotations.csv"
    params:
        cluster_column="SCT_snn_res.0.6"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/annotate_metadata_with_celltypes.R"

rule export_raw_counts:
    input:
        rds="results/seurat/subset_high_quality/sct_harmony.rds"
    output:
        mtx="results/seurat/subset_high_quality/raw_counts.mtx",
        barcodes="results/seurat/subset_high_quality/raw_counts_barcodes.tsv.gz",
        features="results/seurat/subset_high_quality/raw_counts_features.tsv.gz",
        variable_features="results/seurat/subset_high_quality/sct_variable_features.csv"
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/export_raw_counts.R"

rule export_raw_counts_h5ad:
    input:
        mtx="results/seurat/subset_high_quality/raw_counts.mtx",
        barcodes="results/seurat/subset_high_quality/raw_counts_barcodes.tsv.gz",
        features="results/seurat/subset_high_quality/raw_counts_features.tsv.gz",
        variable_features="results/seurat/subset_high_quality/sct_variable_features.csv",
        metadata="results/seurat/subset_high_quality/metadata_annotations.csv",
        harmony="results/seurat/subset_high_quality/harmony_embeddings.csv"
    output:
        h5ad="results/seurat/subset_high_quality/adata_20260717.h5ad"
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/seurat/export_raw_counts_h5ad.py"
