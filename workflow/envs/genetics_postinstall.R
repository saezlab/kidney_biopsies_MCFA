# (C) Sergio Alejandro Gómez
# Post-install step of workflow/envs/genetics.yaml. Run once inside the activated environment. Installs, in this order:
#   1. Matrix 1.7-3 (CRAN archive)
#   2. the Bioconductor 3.21 packages, from the frozen 3.21 repository; BiocParallel 1.42.1 from the 3.21 archive
#   3. car 3.1-3 after Matrix, so that lme4 is built against Matrix 1.7-3
#   4. gchromVAR at GitHub commit e4f33cad4115160ee4bdf16fd625c2fcd0bf3910 (version 0.3.2)
#   5. susieR 0.16.6, the development version at GitHub commit ef213feed2cb82419677661a8c986e1504df2c73 (tag 0.16.6),
#      with its CRAN dependencies at the versions used
# GitHub archives are checked against their SHA-256 before installation. The script stops if any installed version differs
# from the one used for the published results, and writes the full package list to R_package_versions.tsv.

CRAN_REPO <- Sys.getenv("GENETICS_CRAN_REPO", "https://cloud.r-project.org")
BIOC <- c(BioCsoft = "https://bioconductor.org/packages/3.21/bioc",
          BioCann = "https://bioconductor.org/packages/3.21/data/annotation")
options(repos = c(BIOC, CRAN = CRAN_REPO), timeout = 3600, Ncpus = max(1L, parallel::detectCores() - 1L))
lib <- .libPaths()[1]
cat("installing into", lib, "\n")

download <- function(url, file = basename(url)) { f <- file.path(tempdir(), file); download.file(url, f, mode = "wb"); f }
install_url <- function(url) install.packages(download(url), repos = NULL, type = "source", lib = lib)

install_github_commit <- function(url, repo, commit, sha256) {
  f <- download(url, sprintf("%s_%s.tar.gz", basename(repo), commit))
  got <- digest::digest(file = f, algo = "sha256")
  if (!identical(got, sha256))
    stop(repo, " archive SHA-256 is ", got, ", expected ", sha256, ". GitHub archives are not always byte-stable: clone the ",
         "repository, check out ", commit, " and install that directory with remotes::install_local().")
  remotes::install_local(f, lib = lib, dependencies = FALSE, upgrade = "never", build = FALSE)
}
install_cran_version <- function(pkg, version) remotes::install_version(pkg, version, lib = lib, upgrade = "never", dependencies = FALSE)

# 1. Matrix
install_url("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.7-3.tar.gz")

# 2. Bioconductor 3.21
BIOC_PKGS <- c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment", "Biostrings",
               "BSgenome", "rtracklayer", "chromVAR", "TFBSTools", "MotifDb", "motifbreakR", "BSgenome.Hsapiens.UCSC.hg38")
install.packages(setdiff(BIOC_PKGS, rownames(installed.packages(lib.loc = lib))), lib = lib, type = "source",
                 dependencies = c("Depends", "Imports", "LinkingTo"))
install_url("https://bioconductor.org/packages/3.21/bioc/src/contrib/Archive/BiocParallel/BiocParallel_1.42.1.tar.gz")

# 3. car after Matrix
install.packages(c("carData", "abind", "Formula", "pbkrtest", "quantreg", "lme4", "scales"), lib = lib, type = "source")
install_url("https://cran.r-project.org/src/contrib/Archive/car/car_3.1-3.tar.gz")

# 4. gchromVAR
install_github_commit("https://codeload.github.com/caleblareau/gchromVAR/tar.gz/e4f33cad4115160ee4bdf16fd625c2fcd0bf3910",
                      "caleblareau/gchromVAR", "e4f33cad4115160ee4bdf16fd625c2fcd0bf3910",
                      "dfa9442c97f3c4876eeee0864438e22662bc336a27684a1720aa6b0bea853132")

# 5. susieR (development version) and its dependencies
for (p in c("Rcpp", "RcppParallel", "RcppZiggurat")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, type = "source")
SUSIER_DEPS <- c(cpp11 = "0.5.3", cpp11armadillo = "0.5.4", Rfast = "2.1.5.2", matrixStats = "1.5.0", mixsqp = "0.3-54",
                 reshape = "0.8.10", crayon = "1.5.3")
for (p in names(SUSIER_DEPS)) install_cran_version(p, SUSIER_DEPS[[p]])
install_github_commit("https://github.com/stephenslab/susieR/archive/ef213feed2cb82419677661a8c986e1504df2c73.tar.gz",
                      "stephenslab/susieR", "ef213feed2cb82419677661a8c986e1504df2c73",
                      "f36dd4045a448c0b07a1ea8429dc531d573e31fbe822bfe61842244ddd5ddd88")

# ---- versions used ----
VERSIONS <- c(Matrix = "1.7.3", IRanges = "2.42.0", GenomeInfoDb = "1.44.3", GenomicRanges = "1.60.0", SummarizedExperiment = "1.38.1",
              Biostrings = "2.76.0", BSgenome = "1.76.0", rtracklayer = "1.68.0", chromVAR = "1.30.1", BiocParallel = "1.42.1",
              TFBSTools = "1.46.0", MotifDb = "1.50.0", motifbreakR = "2.22.0", BSgenome.Hsapiens.UCSC.hg38 = "1.4.5",
              gchromVAR = "0.3.2", car = "3.1.3", susieR = "0.16.6", Rfast = "2.1.5.2", mixsqp = "0.3.54", matrixStats = "1.5.0",
              data.table = "1.17.8", digest = "0.6.38", jsonlite = "2.0.0", dplyr = "1.2.1", TFMPvalue = "1.0.0",
              ggplot2 = "4.0.3", R.utils = "2.13.0")
got <- vapply(names(VERSIONS), function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "not installed", "")
write.table(installed.packages()[, c("Package", "Version", "LibPath", "Built")], "R_package_versions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
if (as.character(getRversion()) != "4.5.1") stop("R ", getRversion(), " instead of 4.5.1")
bad <- got != VERSIONS
if (any(bad)) stop("versions differ: ", paste(sprintf("%s %s (expected %s)", names(VERSIONS)[bad], got[bad], VERSIONS[bad]), collapse = "; "))
cat("all", length(VERSIONS), "package versions as expected; R_package_versions.tsv written\n")
