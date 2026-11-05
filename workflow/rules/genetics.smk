# (C) Sergio Alejandro Gómez
# Snakemake rules of the genetics analyses: eGFR and CKD enrichment in kidney cell-type chromatin, candidate variants,
# target genes, kidney eQTL/caQTL enrichment and TF motif disruption, for the cohort with and without pt_C5 nuclei.
# Included from workflow/Snakefile as
#     configfile: "config/genetics.yaml"
#     module genetics:
#         snakefile: "rules/genetics.smk"
#         config: config
#     use rule * from genetics as genetics_*
# or run on its own from the repository root:
#     snakemake -s workflow/rules/genetics.smk --configfile config/genetics.yaml --use-conda --cores 8
# Each rule calls one script of workflow/scripts/genetics with --key=value arguments; README.md there describes them.

G = config["genetics"]
IN = G["inputs"]
R = G["results"]
S = "workflow/scripts/genetics"
V = R + "/{version}"
QTL = R + "/qtl_resources"
CKD = R + "/ckd_finemapping"
EXPORT = R + "/full_cohort/peak_matrix_export"
VARIANTS = R + "/egfr_finemapping/egfr_finemapped_variants_hg38.tsv"
ENV = "../envs/genetics.yaml"
ENV_PARQUET = "../envs/genetics_parquet.yaml"

wildcard_constraints:
    version = "with_ptC5|without_ptC5"

def ptc5(wc):
    return "true" if wc.version == "with_ptC5" else "false"

FIGURES = ["egfr_enrichment", "motif_disruption_top30", "motif_pt_vs_tal", "eqtl_enrichment", "caqtl_enrichment",
           "ckd_enrichment", "egfr_enrichment_control_libraries"]


rule all:
    input:
        expand(V + "/figures/{figure}.pdf", version = G["versions"], figure = FIGURES),
        expand(V + "/tables/{table}.tsv", version = G["versions"],
               table = ["control_libraries_13_labels", "candidate_caqtl_summary", "promoter_fraction_by_population",
                        "motif_pt_tal_by_tf", "ckd_enrichment_without_umod"])


# ---- Inputs ----------------------------------------------------------------------------------------------------------

rule egfr_variants:
    input: carma = IN["egfr_carma"], chain = IN["chain_hg19_to_hg38"]
    output: VARIANTS
    log: R + "/logs/egfr_variants.log"
    conda: ENV
    shell: "Rscript {S}/00_prepare_egfr_variants.R --carma={input.carma} --chain={input.chain} "
           "--out_variants={output} > {log} 2>&1"

rule gtex_parquet_to_tsv:
    input: IN["gtex_pairs_parquet"]
    output: QTL + "/gtex_kidney_signif_pairs.tsv.gz"
    log: R + "/logs/gtex_parquet_to_tsv.log"
    conda: ENV_PARQUET
    shell: "python {S}/00_gtex_parquet_to_tsv.py --parquet={input} --out={output} > {log} 2>&1"

rule qtl_resources:
    input: open4gene = IN["open4gene"], gtex_pairs = QTL + "/gtex_kidney_signif_pairs.tsv.gz", gtex_egenes = IN["gtex_egenes"],
           gtf = IN["ensembl_gtf"], basa = IN["basa_caqtl"], chain = IN["chain_hg19_to_hg38"]
    output: open4gene = QTL + "/open4gene_hg38_unique.rds", gtex_evariants = QTL + "/gtex_kidney_unique_evariants_hg38.rds",
            gtex_pairs = QTL + "/gtex_kidney_significant_eqtl_pairs.rds", genes = QTL + "/gene_annotation_hg38_ensembl108.rds",
            caqtl_variants = QTL + "/basa_unique_caqtl_variants_hg38.rds", caqtl_peaks = QTL + "/basa_caqtl_peaks_hg38.rds",
            caqtl_associations = QTL + "/basa_caqtl_associations_hg38.tsv.gz"
    log: R + "/logs/qtl_resources.log"
    conda: ENV
    shell: "Rscript {S}/00_prepare_qtl_resources.R --open4gene={input.open4gene} --gtex_pairs={input.gtex_pairs} "
           "--gtex_egenes={input.gtex_egenes} --gtf={input.gtf} --basa={input.basa} --chain={input.chain} "
           "--out_open4gene={output.open4gene} --out_gtex_evariants={output.gtex_evariants} --out_gtex_pairs={output.gtex_pairs} "
           "--out_gene_annotation={output.genes} --out_caqtl_variants={output.caqtl_variants} --out_caqtl_peaks={output.caqtl_peaks} "
           "--out_caqtl_associations={output.caqtl_associations} > {log} 2>&1"


# ---- Accessibility of the 10 populations -------------------------------------------------------------------------------

rule export_peak_matrix:
    input: IN["h5ad"]
    output: multiext(EXPORT + "/", "pseudobulk_counts_int32_le.bin", "pseudobulk_peaks.tsv", "pseudobulk_groups.tsv",
                     "cells.tsv.gz", "export_summary.json")
    log: R + "/logs/export_peak_matrix.log"
    conda: ENV
    shell: "python {S}/01_export_peak_matrix.py --h5ad={input} --out_dir={EXPORT} > {log} 2>&1"

rule pseudobulk:
    input: rules.export_peak_matrix.output, labels = IN["population_labels"]
    output: R + "/full_cohort/objects/pseudobulk.rds"
    log: R + "/logs/pseudobulk.log"
    conda: ENV
    shell: "Rscript {S}/02_pseudobulk.R --export_dir={EXPORT} --population_labels={input.labels} --out_pseudobulk={output} > {log} 2>&1"

rule accessibility:
    input: rules.pseudobulk.output
    output: peaks = V + "/objects/peak_universe.rds", counts = V + "/objects/accessibility_counts.rds",
            cpm = V + "/objects/accessibility_cpm.rds", cpm1 = V + "/objects/accessibility_cpm_plus1.rds",
            cells = V + "/objects/cell_metadata.rds", composition = V + "/tables/cohort_composition.tsv"
    params: ptc5 = ptc5
    log: V + "/logs/accessibility.log"
    conda: ENV
    shell: "Rscript {S}/03_accessibility.R --pseudobulk={input} --include_ptC5={params.ptc5} --out_peak_universe={output.peaks} "
           "--out_counts={output.counts} --out_cpm={output.cpm} --out_cpm_plus1={output.cpm1} --out_cell_metadata={output.cells} "
           "--out_composition={output.composition} > {log} 2>&1"


# ---- eGFR enrichment and candidate variants -----------------------------------------------------------------------------

rule egfr_enrichment:
    input: peaks = rules.accessibility.output.peaks, cpm1 = rules.accessibility.output.cpm1, variants = VARIANTS
    output: table = V + "/tables/egfr_enrichment.tsv", by_seed = V + "/tables/egfr_enrichment_by_seed.tsv"
    params: ptc5 = ptc5, seeds = G["gchromvar_seeds"]
    log: V + "/logs/egfr_enrichment.log"
    conda: ENV
    shell: "Rscript {S}/04_egfr_enrichment.R --peak_universe={input.peaks} --cpm_plus1={input.cpm1} --variants={input.variants} "
           "--include_ptC5={params.ptc5} --seeds={params.seeds} --out_enrichment={output.table} --out_by_seed={output.by_seed} > {log} 2>&1"

rule specific_peaks:
    input: egfr = rules.egfr_enrichment.output.table, peaks = rules.accessibility.output.peaks, cpm = rules.accessibility.output.cpm
    output: summary = V + "/tables/specific_peaks_summary.tsv", specific = V + "/objects/specific_peaks.rds",
            specificity = V + "/objects/specificity_matrix.rds"
    params: ptc5 = ptc5
    log: V + "/logs/specific_peaks.log"
    conda: ENV
    shell: "Rscript {S}/05_specific_peaks.R --egfr_enrichment={input.egfr} --peak_universe={input.peaks} --cpm={input.cpm} "
           "--include_ptC5={params.ptc5} --out_summary={output.summary} --out_specific_peaks={output.specific} "
           "--out_specificity={output.specificity} > {log} 2>&1"

rule candidate_variants:
    input: variants = VARIANTS, peaks = rules.accessibility.output.peaks, cpm = rules.accessibility.output.cpm,
           specific = rules.specific_peaks.output.specific, specificity = rules.specific_peaks.output.specificity
    output: table = V + "/tables/candidate_variants.tsv", funnel = V + "/tables/candidate_funnel.tsv",
            rds = V + "/objects/candidate_variants.rds", bed = V + "/objects/candidate_variants.bed"
    params: ptc5 = ptc5
    log: V + "/logs/candidate_variants.log"
    conda: ENV
    shell: "Rscript {S}/06_candidate_variants.R --variants={input.variants} --peak_universe={input.peaks} --cpm={input.cpm} "
           "--specific_peaks={input.specific} --specificity={input.specificity} --include_ptC5={params.ptc5} "
           "--out_candidates={output.table} --out_funnel={output.funnel} --out_candidates_rds={output.rds} --out_bed={output.bed} > {log} 2>&1"

rule control_libraries_10_populations:
    input: pseudobulk = rules.pseudobulk.output, peaks = rules.accessibility.output.peaks, variants = VARIANTS
    output: table = V + "/tables/control_libraries_10_populations.tsv", by_seed = V + "/tables/control_libraries_10_populations_by_seed.tsv",
            summary = V + "/tables/control_libraries_10_populations_summary.tsv"
    params: ptc5 = ptc5, seeds = G["gchromvar_seeds"]
    log: V + "/logs/control_libraries_10_populations.log"
    conda: ENV
    shell: "Rscript {S}/07_control_libraries_10_populations.R --pseudobulk={input.pseudobulk} --peak_universe={input.peaks} "
           "--variants={input.variants} --include_ptC5={params.ptc5} --seeds={params.seeds} --out_comparison={output.table} "
           "--out_by_seed={output.by_seed} --out_summary={output.summary} > {log} 2>&1"

rule control_libraries_13_labels:
    input: pseudobulk = rules.pseudobulk.output, labels = IN["labels_13"], variants = VARIANTS, by_seed = rules.egfr_enrichment.output.by_seed
    output: table = V + "/tables/control_libraries_13_labels.tsv", nuclei = V + "/tables/control_libraries_13_labels_nuclei.tsv",
            libraries = V + "/tables/control_libraries_13_labels_libraries.tsv", summary = V + "/tables/control_libraries_13_labels_summary.tsv"
    params: ptc5 = ptc5
    log: V + "/logs/control_libraries_13_labels.log"
    conda: ENV
    shell: "Rscript {S}/08_control_libraries_13_labels.R --pseudobulk={input.pseudobulk} --labels13={input.labels} --variants={input.variants} "
           "--egfr_by_seed={input.by_seed} --include_ptC5={params.ptc5} --out_comparison={output.table} --out_nuclei={output.nuclei} "
           "--out_libraries={output.libraries} --out_summary={output.summary} > {log} 2>&1"


# ---- Target genes, eQTL and caQTL ---------------------------------------------------------------------------------------

rule target_genes:
    input: table = rules.candidate_variants.output.table, rds = rules.candidate_variants.output.rds,
           open4gene = rules.qtl_resources.output.open4gene, gtex_evariants = rules.qtl_resources.output.gtex_evariants,
           gtex_pairs = rules.qtl_resources.output.gtex_pairs, genes = rules.qtl_resources.output.genes
    output: long = V + "/tables/candidate_target_gene_long.tsv", summary = V + "/tables/candidate_target_gene_summary.tsv",
            counts = V + "/tables/candidate_target_gene_evidence_counts.tsv", convergent = V + "/tables/convergent_candidate_gene_links.tsv",
            gtex = V + "/tables/candidate_gtex_exact_eqtl_matches.tsv", open4gene = V + "/tables/candidate_open4gene_matches.tsv",
            context = V + "/tables/candidate_open4gene_atlas_peak_overlap.tsv"
    log: V + "/logs/target_genes.log"
    conda: ENV
    shell: "Rscript {S}/09_target_genes.R --candidates={input.table} --candidates_rds={input.rds} --open4gene={input.open4gene} "
           "--gtex_evariants={input.gtex_evariants} --gtex_pairs={input.gtex_pairs} --gene_annotation={input.genes} "
           "--out_long={output.long} --out_summary={output.summary} --out_counts={output.counts} --out_convergent={output.convergent} "
           "--out_gtex_matches={output.gtex} --out_open4gene_matches={output.open4gene} --out_open4gene_peak_context={output.context} > {log} 2>&1"

rule eqtl_enrichment:
    input: peaks = rules.accessibility.output.peaks, cpm = rules.accessibility.output.cpm, specificity = rules.specific_peaks.output.specificity,
           egfr = rules.egfr_enrichment.output.table, gtex_evariants = rules.qtl_resources.output.gtex_evariants,
           genes = rules.qtl_resources.output.genes
    output: table = V + "/tables/eqtl_enrichment.tsv", coefficients = V + "/tables/eqtl_enrichment_coefficients.tsv",
            counts = V + "/tables/eqtl_peak_counts.tsv", annotation = V + "/objects/peak_annotation_eqtl.rds"
    log: V + "/logs/eqtl_enrichment.log"
    conda: ENV
    shell: "Rscript {S}/10_eqtl_enrichment.R --peak_universe={input.peaks} --cpm={input.cpm} --specificity={input.specificity} "
           "--egfr_enrichment={input.egfr} --gtex_evariants={input.gtex_evariants} --gene_annotation={input.genes} "
           "--out_enrichment={output.table} --out_coefficients={output.coefficients} --out_counts={output.counts} "
           "--out_peak_annotation={output.annotation} > {log} 2>&1"

rule caqtl_enrichment:
    input: annotation = rules.eqtl_enrichment.output.annotation, caqtl_peaks = rules.qtl_resources.output.caqtl_peaks,
           caqtl_variants = rules.qtl_resources.output.caqtl_variants, egfr = rules.egfr_enrichment.output.table,
           eqtl = rules.eqtl_enrichment.output.table
    output: table = V + "/tables/caqtl_enrichment.tsv", coefficients = V + "/tables/caqtl_enrichment_coefficients.tsv",
            counts = V + "/tables/caqtl_peak_counts.tsv", annotation = V + "/objects/peak_annotation_caqtl.rds"
    log: V + "/logs/caqtl_enrichment.log"
    conda: ENV
    shell: "Rscript {S}/11_caqtl_enrichment.R --peak_annotation={input.annotation} --caqtl_peaks={input.caqtl_peaks} "
           "--caqtl_variants={input.caqtl_variants} --egfr_enrichment={input.egfr} --eqtl_enrichment={input.eqtl} "
           "--out_enrichment={output.table} --out_coefficients={output.coefficients} --out_counts={output.counts} "
           "--out_peak_annotation={output.annotation} > {log} 2>&1"

rule candidate_caqtl:
    input: candidates = rules.candidate_variants.output.table, genes = rules.target_genes.output.long,
           motifs = V + "/tables/motif_disruption_by_variant.tsv", caqtl_variants = rules.qtl_resources.output.caqtl_variants,
           caqtl_peaks = rules.qtl_resources.output.caqtl_peaks, caqtl_associations = rules.qtl_resources.output.caqtl_associations
    output: summary = V + "/tables/candidate_caqtl_summary.tsv", exact = V + "/tables/candidate_caqtl_exact_matches.tsv",
            multilayer = V + "/tables/multilayer_candidate_gene_evidence.tsv", high = V + "/tables/multilayer_high_support_examples.tsv"
    log: V + "/logs/candidate_caqtl.log"
    conda: ENV
    shell: "Rscript {S}/12_candidate_caqtl.R --candidates={input.candidates} --target_genes={input.genes} --motif_disruption={input.motifs} "
           "--caqtl_variants={input.caqtl_variants} --caqtl_peaks={input.caqtl_peaks} --caqtl_associations={input.caqtl_associations} "
           "--out_summary={output.summary} --out_exact_matches={output.exact} --out_multilayer={output.multilayer} "
           "--out_multilayer_high={output.high} > {log} 2>&1"

rule promoter_fraction:
    input: peaks = rules.accessibility.output.peaks, specificity = rules.specific_peaks.output.specificity, genes = rules.qtl_resources.output.genes
    output: V + "/tables/promoter_fraction_by_population.tsv"
    log: V + "/logs/promoter_fraction.log"
    conda: ENV
    shell: "Rscript {S}/13_promoter_fraction.R --peak_universe={input.peaks} --specificity={input.specificity} "
           "--gene_annotation={input.genes} --out={output} > {log} 2>&1"


# ---- TF motif disruption ------------------------------------------------------------------------------------------------

rule motif_disruption:
    input: candidates = rules.candidate_variants.output.table, bed = rules.candidate_variants.output.bed, genes = rules.qtl_resources.output.genes
    output: events = V + "/tables/motif_disruption_events.tsv", by_tf = V + "/tables/motif_disruption_by_tf.tsv",
            by_variant = V + "/tables/motif_disruption_by_variant.tsv", top30 = V + "/tables/motif_disruption_top30_tfs.tsv",
            settings = V + "/tables/motifbreakr_settings.tsv", peak_class = V + "/tables/candidate_peak_class_counts.tsv",
            results = V + "/objects/motifbreakr_results.rds"
    # the worker count is passed as a parameter, not as {threads}: Snakemake lowers threads to --cores, and the number of
    # workers sets the chunking of the candidates and with it the row order of the saved motifbreakR object
    params: workers = G["motifbreakr_workers"]
    threads: G["motifbreakr_workers"]
    log: V + "/logs/motif_disruption.log"
    conda: ENV
    shell: "Rscript {S}/14_motif_disruption.R --candidates={input.candidates} --bed={input.bed} --gene_annotation={input.genes} "
           "--n_workers={params.workers} --out_events={output.events} --out_by_tf={output.by_tf} --out_by_variant={output.by_variant} "
           "--out_top30={output.top30} --out_settings={output.settings} --out_peak_class={output.peak_class} "
           "--out_results={output.results} > {log} 2>&1"

rule motif_pt_vs_tal:
    input: candidates = rules.candidate_variants.output.table, by_variant = rules.motif_disruption.output.by_variant,
           by_tf = rules.motif_disruption.output.by_tf
    output: groups = V + "/tables/motif_pt_tal_group_counts.tsv", by_tf = V + "/tables/motif_pt_tal_by_tf.tsv",
            top30 = V + "/tables/motif_pt_tal_top30_tfs.tsv", similarity = V + "/tables/motif_pt_tal_similarity.tsv",
            topk = V + "/tables/motif_pt_tal_topk_overlap.tsv", load = V + "/tables/motif_tf_load_by_group.tsv"
    log: V + "/logs/motif_pt_vs_tal.log"
    conda: ENV
    shell: "Rscript {S}/15_motif_pt_vs_tal.R --candidates={input.candidates} --by_variant={input.by_variant} --by_tf={input.by_tf} "
           "--out_group_counts={output.groups} --out_by_tf={output.by_tf} --out_top30={output.top30} --out_similarity={output.similarity} "
           "--out_topk={output.topk} --out_tf_load={output.load} > {log} 2>&1"


# ---- CKD: fine-mapping (SuSiE-RSS) and enrichment -----------------------------------------------------------------------

rule ckd_loci_and_ld:
    input: gwas = IN["ckd_gwas"], panel = multiext(IN["ld_panel"], ".bed", ".bim", ".fam")
    output: loci = CKD + "/ckd_loci.tsv", variants = CKD + "/ckd_locus_variants.tsv", genotypes = CKD + "/ckd_locus_genotypes.rds"
    params: panel = IN["ld_panel"]
    log: R + "/logs/ckd_loci_and_ld.log"
    conda: ENV
    shell: "Rscript {S}/16_ckd_loci_and_ld.R --ckd_gwas={input.gwas} --ld_panel={params.panel} --out_loci={output.loci} "
           "--out_variants={output.variants} --out_genotypes={output.genotypes} > {log} 2>&1"

rule ckd_finemapping:
    input: loci = rules.ckd_loci_and_ld.output.loci, variants = rules.ckd_loci_and_ld.output.variants,
           genotypes = rules.ckd_loci_and_ld.output.genotypes, chain = IN["chain_hg19_to_hg38"]
    output: variants = CKD + "/ckd_finemapped_variants_hg38.rds", qc = CKD + "/ckd_locus_qc.tsv", pips = CKD + "/ckd_pip_by_locus.tsv"
    threads: G["ckd_finemapping_threads"]
    log: R + "/logs/ckd_finemapping.log"
    conda: ENV
    shell: "Rscript {S}/17_ckd_finemapping.R --loci={input.loci} --variants={input.variants} --genotypes={input.genotypes} "
           "--chain={input.chain} --threads={threads} --out_variants={output.variants} --out_locus_qc={output.qc} "
           "--out_pip_by_locus={output.pips} > {log} 2>&1"

rule ckd_enrichment:
    input: peaks = rules.accessibility.output.peaks, cpm1 = rules.accessibility.output.cpm1,
           variants = rules.ckd_finemapping.output.variants, qc = rules.ckd_finemapping.output.qc
    output: table = V + "/tables/ckd_enrichment.tsv", by_seed = V + "/tables/ckd_enrichment_by_seed.tsv",
            without_umod = V + "/tables/ckd_enrichment_without_umod.tsv", without_umod_by_seed = V + "/tables/ckd_enrichment_without_umod_by_seed.tsv"
    params: seeds = G["gchromvar_seeds"]
    log: V + "/logs/ckd_enrichment.log"
    conda: ENV
    shell: "Rscript {S}/18_ckd_enrichment.R --peak_universe={input.peaks} --cpm_plus1={input.cpm1} --ckd_variants={input.variants} "
           "--ckd_locus_qc={input.qc} --seeds={params.seeds} --out_enrichment={output.table} --out_by_seed={output.by_seed} "
           "--out_without_umod={output.without_umod} --out_without_umod_by_seed={output.without_umod_by_seed} > {log} 2>&1"


# ---- Figures --------------------------------------------------------------------------------------------------------------

rule figures:
    input: rules.egfr_enrichment.output.table, rules.candidate_variants.output.table, rules.motif_disruption.output.top30,
           rules.motif_pt_vs_tal.output.by_tf, rules.motif_pt_vs_tal.output.similarity, rules.eqtl_enrichment.output.table,
           rules.caqtl_enrichment.output.table, rules.ckd_enrichment.output.table, rules.control_libraries_10_populations.output.table,
           rules.control_libraries_10_populations.output.summary
    output: expand(V + "/figures/{figure}.{ext}", figure = FIGURES, ext = ["pdf", "png"], allow_missing = True)
    params: ptc5 = ptc5, tables = V + "/tables", out_dir = V + "/figures"
    log: V + "/logs/figures.log"
    conda: ENV
    shell: "Rscript {S}/19_figures.R --tables={params.tables} --include_ptC5={params.ptc5} --out_dir={params.out_dir} > {log} 2>&1"
