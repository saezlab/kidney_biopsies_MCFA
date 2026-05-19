rule all_seurat:
    input:
        "results/seurat/merged_seurat.rds",
        "results/seurat/sct_harmony.rds",
        "results/seurat/metadata.csv",
        "results/seurat/pca_embeddings.csv",
        "results/seurat/harmony_embeddings.csv",
        "results/seurat/umap_pre_harmony_embeddings.csv",
        "results/seurat/umap_embeddings.csv",
        "results/seurat/pca_stdev.csv",
        "plots/seurat/clustree_sct_harmony.pdf",
        "plots/seurat/umap_sct_resolutions.pdf",
        "results/seurat/session_info.csv"

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
        clustree_plot="plots/seurat/clustree_sct_harmony.pdf",
        session_info_csv="results/seurat/session_info.csv",
    container:
        config["containers"]["seurat"]
    script:
        "../scripts/seurat/run_sct_harmony_clustering.R"