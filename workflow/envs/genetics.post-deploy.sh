#!/bin/sh
# (C) Sergio Alejandro Gómez
# Snakemake (--use-conda) runs this script inside the environment right after creating it from genetics.yaml. It installs
# the R packages that conda does not provide (genetics_postinstall.R). Snakemake starts it from its working directory,
# the repository root; the package list R_package_versions.tsv is written into the environment folder.
set -eu
script="$(pwd)/workflow/envs/genetics_postinstall.R"
cd "${CONDA_PREFIX:-.}"
Rscript "$script"
