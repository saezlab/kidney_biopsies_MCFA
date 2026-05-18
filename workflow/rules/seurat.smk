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
