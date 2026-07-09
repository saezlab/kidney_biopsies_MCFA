rule download_gencode:
    output:
        gtf="data/reference/gencode.v38.primary_assembly.annotation.gtf.gz"
    params:
        url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
        )
    script:
        "../scripts/glue/download_gencode.py"


#snakemake --cores 1 results/glue/coembedding.csv --software-deployment-method conda --rerun-triggers mtime --resources gpu=1
rule run_glue:
    input:
        rna_h5ad="results/seurat/subset_high_quality/adata_20260717.h5ad",
        atac_h5ad="data/ATAC/peak_mat2.h5ad",
        atac_lsi_csv="data/ATAC/harmony_lsi.csv",
        gtf="data/reference/gencode.v38.primary_assembly.annotation.gtf.gz"
    output:
        csv="results/glue/coembedding.csv",
        model_dir=directory("results/glue/model")
    params:
        cuda_device=config["glue"]["cuda_device"]
    resources:
        gpu=1
    conda:
        "../envs/glue.yaml"
    script:
        "../scripts/glue/run_glue.py"


rule plot_coembedding:
    input:
        csv="results/glue/coembedding.csv",
        rna_annotations="results/seurat/subset_high_quality/metadata_annotations.csv",
        atac_annotations="data/ATAC/label_transfer_RNA_ATAC_V3_clean.csv",
        clusters="results/glue/coembedding_clusters.csv",
        glue_annotations="config/GLUE_annotations.csv"
    output:
        plot="plots/glue/coembedding_umap.pdf",
        coords="results/glue/coembedding_umap.csv",
        annotated_plot="plots/glue/coembedding_umap_annotations.pdf"
    container:
        config["containers"]["scanpy"]
    script:
        "../scripts/glue/plot_coembedding.py"