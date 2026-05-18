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
