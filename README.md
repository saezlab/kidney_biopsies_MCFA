# kidney_biopsies_MCFA — Multicellular Factor Analysis (Snakemake)

This repository contains a Snakemake workflow implementing a Multicellular Factor Analysis (MCFA) of single-nucleus RNA-seq data from a human kidney biopsy dataset.
For the analysis of snATAC-seq data, see: github.com/saezlab/kidney_biopsies_ArchR

The MCFA analysis is organized within this repository as a reproducible Snakemake workflow. Also contained within this repository are the scripts used to produce figures from the manuscript, located in `workflow/scripts/figures`.

## Files necessary to produce manuscript figures

All files listed below as being contained in the `results` folder are produced by the Snakemake workflow within this repository. All required input files in the `data` folder are produced externally to this repository (see note on input data below).

- 1b_modality_heatmap.R
    - 'data/data_overview/modality_with_duplicates~heatmap.csv'
- 1c_1e_followup_boxplot_and_eGFR_vs_slope.R
    - 'data/data_overview/follow_up~boxplot.csv'
    - 'data/metadata/patient_metadata_EDA.csv'
- 1d_plot_endpoints.R
    - 'data/data_overview/endpoint_atlas~heatmap.csv'
    - 'data/data_overview/proteinuria_over_time.csv'
- 1g_ext1e_umap_and_label_transfer.R
    - 'data/data_overview/proteinuria_over_time.csv'
    - 'data/ATAC/coordinates.csv'
- 2b_2c_plot_MDS.R
    - 'results/MCFA/factors_obs~all~5.csv'
- 2d_plot_metadata_associations.R
    - 'results/MCFA/metadata_pearson_correlations.csv'
    - 'results/MCFA/metadata_anova.csv'
    - 'results/MCFA/clr_pearson_correlations.csv'
- 2e_eGFR_vs_Factor1_vs_aPT_to_PT_ratio.R
    - 'results/composition/prop~main_class~atlas.csv'
    - 'results/MCFA/factors_obs~hvg~5.csv'
- 2f_PT_aPT_marker_motifs.R
    - 'data/ATAC/top_disrupted_tfs.csv'
    - 'data/ATAC/enrich_motifs.rds'
- 3a_plot_toploadings.R
    - 'results/MCFA/loadings~all~5.csv'
- 5a_5b_5c_plot_factor2_enrichments.R
    - 'results/MCFA/acts_collectri.csv'
    - 'results/MCFA/acts_progeny.csv'
    - 'results/MCFA/hallmark_Factor2.csv'
- 5g_factor2_enrichment_score.R
    - 'results/RelA_KO/factor2_enrichment.csv'
- 5i_5j_compare_factor2_and_KO.R
    - 'data/RelA_KO/DESEQ2.EC_AB.versus.EC_AA.txt'
    - 'results/RelA_KO/collectri_TFs.csv'
    - 'results/MCFA/acts_collectri.csv'
    - 'results/RelA_KO/progeny_pathways.csv'
    - 'results/MCFA/acts_progeny.csv'
- ext1a_upset_plot.R
    - 'data/data_overview/modality~heatmap.csv'
- ext1b_age_at_biopsy.R
    - 'data/data_overview/age~boxplot.csv'
- ext1d_plot_marker_motifs.R
    - 'data/ATAC/enrich_motifs.rds'
- ext1f_ext1g_cluster_proportions_nuclei_counts.R
    - 'data/data_overview/patient_metadata_including_xenium_and_duplicate_samples.csv'
    - 'data/annotations/biopsydata_withkpmp_final_annotations(in).csv'
    - 'data/annotations/label_transfer_RNA_ATAC_V3_clean.csv'
- ext2a_ext2d_factor_robustness.R
    - 'results/MCFA/factor_robustness.csv'
    - 'results/MCFA/robustness_r2_by_sample~all.csv'
- ext2b_plot_r2.R
    - 'results/MCFA/r2_total~all~5.csv'
- ext2c_factor_score_distributions.R
    - 'results/MCFA/factors_obs~all~5.csv'
- ext2e_clr_compositions_boxplot.R
    - 'results/composition/clr~main_class~atlas.csv'
    - 'data/metadata/patient_metadata_EDA.csv'

## Workflow overview

- Workflow driven by Snakemake, with rules in `workflow/rules/` split into core subtasks: `MCFA`, `composition` analysis and `RelA_KO` analysis.
- Analysis scripts organized under `workflow/scripts/` (Python and R), split into core subtasks: `MCFA/`, `composition/`, `RelA_KO/`.
- Outputs (results and plots) are generated to `results/` and `plots/` (git-ignored).
- Configuration at `config/config.yaml` controls workflow parameters, container images and paths.

## Requirements 
- This analysis uses Snakemake 9.21.0. To reproduce the analysis, first create an environment with Snakemake v9 using `conda env create -f environment.yml`
- Cluster access (Slurm). A `config/slurm/` folder contains example jobscript and utilities for submitting to Slurm, which may need adapting to fulfill job submission requirements of your Slurm system. %TODO: update
- Singularity images in SIF format (see below). %TODO: update
- Data (see below).

## Configuration

Primary configuration files live under `config/`:

- `config/config.yaml` — main workflow parameters (samples, settings, thresholds, etc.).
- `config/slurm/` — helper scripts and `slurm-jobscript.sh` for cluster execution.

Before running, inspect and edit `config/config.yaml` to point to the correct input paths for your system.

## Singularity images

The respository references several Singularity images to ensure full reproducibility of this workflow. This repository will be updated with access links once the files are readily available on Zenodo.

## Data

The repository references several input files contained in `data`. This repository will be updated with access links once the files are readily available on Zenodo.

## Running the workflow (examples)

1) Dry run (see what would be executed):

```zsh
snakemake -n
```

2) Run on a Slurm cluster using the provided Slurm scripts (example):

```zsh
# Adapt config/slurm/config.yaml for your cluster, then:
snakemake --profile config/slurm
```
## Acknowledgements
The authors acknowledge support by the state of Baden-Württemberg through bwHPC and the German Research Foundation (DFG) through grant INST 35/1597-1 FUGG, as well as the data storage service SDS@hd supported by the Ministry of Science, Research and the Arts Baden-Württemberg (MWK) and the German Research Foundation (DFG) through grant INST 35/1503-1 FUGG. Charlotte Boys gratefully acknowledges DFG funding through the Clinical Research Unit 5011 InteraKD (Project ID: 445703531).

## License

This project is distributed under the GNU General Public License v3.0. See `LICENSE` for full details.
