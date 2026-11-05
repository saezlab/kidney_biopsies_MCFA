# (C) Sergio Alejandro Gómez
# Image for the genetics analyses: the conda environment of workflow/envs/genetics.yaml plus its post-install step.
# Build from the repository root:
#     docker build -f docker/genetics/Dockerfile -t kidney_genetics:1.0 .
# For a cluster:  singularity build kidney_genetics_1.0.sif docker-daemon://kidney_genetics:1.0
FROM condaforge/miniforge3:24.11.3-0

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 TZ=UTC

# The figures use Arial. Without it a metric-compatible substitute is used and the PNGs differ slightly in their text.
RUN echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
 && apt-get update && apt-get install -y --no-install-recommends software-properties-common \
 && add-apt-repository -y multiverse && apt-get update \
 && apt-get install -y --no-install-recommends ttf-mscorefonts-installer fontconfig && fc-cache -f \
 && rm -rf /var/lib/apt/lists/*

RUN conda config --system --set channel_priority strict
COPY workflow/envs/genetics.yaml /opt/genetics/genetics.yaml
COPY workflow/envs/genetics_parquet.yaml /opt/genetics/genetics_parquet.yaml
RUN mamba env create -f /opt/genetics/genetics.yaml && mamba env create -f /opt/genetics/genetics_parquet.yaml && conda clean -afy

COPY workflow/envs/genetics_postinstall.R /opt/genetics/genetics_postinstall.R
SHELL ["conda", "run", "--no-capture-output", "-n", "genetics", "/bin/bash", "-c"]
RUN cd /opt/genetics && Rscript genetics_postinstall.R \
 && conda env export -n genetics > /opt/genetics/conda_env_export.yaml \
 && conda list -n genetics --explicit > /opt/genetics/conda_explicit.txt

# The code is mounted from the repository at run time; nothing analysis-specific is baked into the image.
ENV PATH=/opt/conda/envs/genetics/bin:$PATH
LABEL org.opencontainers.image.title="kidney_genetics" \
      org.opencontainers.image.description="R 4.5.1 + Bioconductor 3.21 (chromVAR, gchromVAR, motifbreakR), development susieR 0.16.6 and Python 3.10 for the genetics analyses of the kidney biopsy snATAC-seq study"
WORKDIR /work
