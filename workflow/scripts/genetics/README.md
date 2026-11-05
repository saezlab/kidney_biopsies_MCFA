# Genetics analyses

These scripts ask which kidney cell populations potentially carry the genetic signal of kidney function. They test which of the 10
snATAC-seq populations are enriched for fine-mapped eGFR and CKD variants. They then select the candidate variants that lie in
chromatin specific to the enriched populations. For these candidates they link putative target genes, test kidney eQTL and
caQTL enrichment, and predict TF motif disruption.

## Running

From the repository root:

```
snakemake -s workflow/rules/genetics.smk --configfile config/genetics.yaml --use-conda --cores 8
```

With `--use-conda`, Snakemake creates the two environments of `workflow/envs/` and then runs `genetics.post-deploy.sh`,
which installs the R packages that conda does not provide (`genetics_postinstall.R`: Bioconductor 3.21, gchromVAR and
the development version of susieR). To run the scripts by hand, create the environment once:

```
mamba env create --strict-channel-priority -f workflow/envs/genetics.yaml
mamba env create --strict-channel-priority -f workflow/envs/genetics_parquet.yaml
conda activate genetics && Rscript workflow/envs/genetics_postinstall.R
```

`docker/genetics/Dockerfile` builds the same environment as an image. Every script can also be run on its own: its header
states the aim, the inputs, the outputs and an example command, and every path and parameter is a `--key=value` argument.
`--include_ptC5=true|false` selects the version.

## Steps

| script | aim | main outputs (per version, `results/genetics/<version>/`) |
|---|---|---|
| `00_prepare_egfr_variants.R` | CARMA eGFR variants (GRCh37) lifted to hg38 as single bases, REF/ALT oriented to hg38 | `egfr_finemapping/egfr_finemapped_variants_hg38.tsv` |
| `00_gtex_parquet_to_tsv.py`, `00_prepare_qtl_resources.R` | Open4Gene links, GTEx v10 kidney-cortex eQTLs, bASA caQTLs and Ensembl 108 genes, in hg38 | `qtl_resources/` |
| `01_export_peak_matrix.py`, `02_pseudobulk.R` | cluster × library pseudobulk of the ArchR PeakMatrix (all 63,891 nuclei) | `full_cohort/objects/pseudobulk.rds` |
| `03_accessibility.R` | counts and CPM of the 10 populations in the 37-library cohort; cohort composition | `objects/accessibility_*.rds`, `tables/cohort_composition.tsv` |
| `04_egfr_enrichment.R` | weighted gchromVAR of the eGFR PIPs; seed 2024 reported, seeds 1–5 for spread | `tables/egfr_enrichment.tsv` |
| `05_specific_peaks.R` | cell-type-specific peaks (CPM > 0 and at least twice the mean of the other populations) | `objects/specificity_matrix.rds`, `tables/specific_peaks_summary.tsv` |
| `06_candidate_variants.R` | candidates: PIP > 0.1 in a peak specific to an enriched population | `tables/candidate_variants.tsv`, `tables/candidate_funnel.tsv` |
| `07_control_libraries_10_populations.R` | eGFR enrichment without the five control libraries (37 vs 32 libraries) | `tables/control_libraries_10_populations*.tsv` |
| `08_control_libraries_13_labels.R` | the same on the 13 labels, without control, AKI-control, KIN and BKVAN libraries | `tables/control_libraries_13_labels*.tsv` |
| `09_target_genes.R` | target genes: Open4Gene links and exact GTEx kidney-cortex eQTL matches; nearest gene as proximity only | `tables/candidate_target_gene_*.tsv` |
| `10_eqtl_enrichment.R`, `11_caqtl_enrichment.R` | eQTL and caQTL enrichment of specific peaks: logistic models adjusted for accessibility, GC and TSS distance, with cluster-robust SEs | `tables/eqtl_enrichment.tsv`, `tables/caqtl_enrichment.tsv` |
| `12_candidate_caqtl.R` | caQTL support of each candidate and the evidence layers of each candidate–gene pair | `tables/candidate_caqtl_summary.tsv` |
| `13_promoter_fraction.R` | share of promoter-proximal peaks (TSS ± 2 kb) among specific and other peaks | `tables/promoter_fraction_by_population.tsv` |
| `14_motif_disruption.R` | motifbreakR on the candidates (692 JASPAR 2022 matrices, `method = "ic"`, `filterp = TRUE`, `threshold = 0.85`, strong effects only) | `tables/motif_disruption_*.tsv` |
| `15_motif_pt_vs_tal.R` | motif disruption of PT- versus TAL-attributed candidates (descriptive) | `tables/motif_pt_tal_*.tsv` |
| `16_ckd_loci_and_ld.R`, `17_ckd_finemapping.R` | CKD loci, 1000 Genomes EUR reference genotypes, SuSiE-RSS fine-mapping with locus QC, PIPs lifted to hg38 | `ckd_finemapping/` |
| `18_ckd_enrichment.R` | weighted gchromVAR of the CKD PIPs, all eligible loci and without the UMOD locus | `tables/ckd_enrichment*.tsv` |
| `19_figures.R` | figures of the version | `figures/*.pdf`, `figures/*.png` |

Shared helpers are in `R/` (argument parsing and table I/O in `utils.R`, gchromVAR in `gchromvar.R`, QTL models in `qtl.R`,
motifbreakR in `motif.R`, PLINK reading in `plink.R`, population constants in `accessibility.R`).

## Inputs

Paths are set in `config/genetics.yaml`.

| input | source |
|---|---|
| `peak_mat2.h5ad` | ArchR PeakMatrix of the study (63,891 nuclei of 43 libraries × 447,919 peaks), exported from the ArchR project |
| `resources/genetics/atac_population_labels.csv`, `atac_13_labels.csv` | cluster annotation of the study (in this repository) |
| `finemapping_carma.csv` | CARMA fine-mapping of eGFR in European ancestry, Liu et al. 2025, *Science*, Table S7 (GRCh37), saved as semicolon-separated text with decimal commas |
| `hg19ToHg38.over.chain` | UCSC, `hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz` (unzipped) |
| `Open4Gene_62K_Cells_Peak_to_Gene_Significant_Associations.txt.gz` | Open4Gene kidney multiome peak-to-gene links, Liu et al. 2025, Tables S20–S21 (hg19) |
| `Kidney_Cortex.v10.eQTLs.signif_pairs.parquet`, `Kidney_Cortex.v10.eGenes.txt.gz` | GTEx v10 single-tissue eQTLs (GTEx Portal downloads, hg38) |
| `bASA_RASQUAL_82_Kidney_Significant_Associations.txt.gz` | kidney caQTLs (bASA/RASQUAL, 82 kidneys), Liu et al. 2025, Tables S14–S15 (hg19) |
| `Homo_sapiens.GRCh38.108.gtf` | Ensembl release 108 (unzipped) |
| `CKD_overall_EA_JW_20180223_nstud23.dbgap.txt` | CKDGen CKD meta-analysis, European ancestry, Wuttke et al. 2019, *Nat Genet* (ckdgen.imbi.uni-freiburg.de) |
| `g1000_eur.{bed,bim,fam}` | 1000 Genomes Phase 3 EUR reference in PLINK format (503 individuals), from the MAGMA auxiliary files (`g1000_eur.zip`, ctg.cncr.nl/software/magma) |

The hg38 genome (BSgenome.Hsapiens.UCSC.hg38 1.4.5) and the JASPAR 2022 matrices (MotifDb 1.50.0) come with the environment.

## Supplementary tables

Each supplementary table is a selection of columns from these outputs, both versions side by side:

| table | outputs |
|---|---|
| eGFR enrichment | `egfr_enrichment.tsv`, `egfr_enrichment_by_seed.tsv`, `cohort_composition.tsv`, `specific_peaks_summary.tsv` |
| candidate selection | `candidate_funnel.tsv` |
| candidate variants | `candidate_variants.tsv`, `candidate_target_gene_summary.tsv`, CPM of the candidate peaks from `objects/accessibility_cpm.rds` |
| variant–gene pairs | `candidate_target_gene_long.tsv` |
| TF motif disruption, per TF | `motif_disruption_by_tf.tsv`, `motif_pt_tal_by_tf.tsv` |
| TF motif disruption, per candidate | `motif_disruption_by_variant.tsv` |
| eQTL and caQTL enrichment | `eqtl_enrichment.tsv`, `caqtl_enrichment.tsv` |
| caQTL support of the candidates | `candidate_caqtl_summary.tsv` |
| CKD enrichment | `ckd_enrichment.tsv`, `ckd_enrichment_by_seed.tsv`, `ckd_enrichment_without_umod.tsv`, `ckd_enrichment_without_umod_by_seed.tsv` |

## Reproducibility notes

- Seeds are fixed: gchromVAR background peaks (2024; spread over seeds 1–5), the motifbreakR worker cluster (2024) and
  every SuSiE-RSS fit (2024).
- CKD fine-mapping uses the development version of susieR, 0.16.6 (GitHub commit ef213feed2cb82419677661a8c986e1504df2c73),
  which provides the finite-reference and LD-mismatch corrections (`R_finite`, `R_mismatch = "eb"`).
