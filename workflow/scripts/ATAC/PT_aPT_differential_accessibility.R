if (exists("snakemake")) {
  project_path <- snakemake@input[["project_path"]]
  marker_peaks_rds <- snakemake@output[["marker_peaks"]]
  apt_open_enrich_motifs_rds <- snakemake@output[["apt_open_enrich_motifs"]]
  pt_open_enrich_motifs_rds <- snakemake@output[["pt_open_enrich_motifs"]]
  group_by <- snakemake@params[["group_by"]]
  use_group <- snakemake@params[["use_group"]]
  background_group <- snakemake@params[["background_group"]]
  threads <- snakemake@threads
} else {
  project_path <- "data/ATAC/proj_peaks2"
  marker_peaks_rds <- "results/ATAC/PT_aPT_differential_accessibility_marker_peaks.rds"
  apt_open_enrich_motifs_rds <- "results/ATAC/aPT_open_motif_enrichment.rds"
  pt_open_enrich_motifs_rds <- "results/ATAC/PT_open_motif_enrichment.rds"
  group_by <- "cluster"
  use_group <- "pt_C6"
  background_group <- "pt_C3"
  threads <- 1L
}

suppressPackageStartupMessages({
  library(ArchR)
})

threads <- max(1L, as.integer(threads))

cat("project_path: ", project_path, "\n", sep = "")
cat("group_by: ", group_by, "\n", sep = "")
cat("use_group: ", use_group, "\n", sep = "")
cat("background_group: ", background_group, "\n", sep = "")
cat("marker_peaks_rds: ", marker_peaks_rds, "\n", sep = "")
cat("apt_open_enrich_motifs_rds: ", apt_open_enrich_motifs_rds, "\n", sep = "")
cat("pt_open_enrich_motifs_rds: ", pt_open_enrich_motifs_rds, "\n", sep = "")
cat("threads: ", threads, "\n", sep = "")

addArchRGenome("hg38")
addArchRThreads(threads = threads)

proj <- loadArchRProject(path = project_path)

if (!"Motif" %in% names(proj@peakAnnotation)) {
  stop(
    "ArchR project is missing required peak annotation 'Motif'. ",
    "Run the motifAnnotations rule before differential accessibility. ",
    "Available peak annotations: ",
    paste(names(proj@peakAnnotation), collapse = ", "),
    call. = FALSE
  )
}

cell_col_data <- getCellColData(
  ArchRProj = proj,
  select = group_by,
  drop = FALSE
)

if (!group_by %in% colnames(cell_col_data)) {
  stop("ArchR project metadata is missing group_by column: ", group_by, call. = FALSE)
}

group_values <- as.character(cell_col_data[[group_by]])
missing_groups <- setdiff(c(use_group, background_group), unique(group_values))

if (length(missing_groups) > 0) {
  stop(
    "Requested group(s) not found in ",
    group_by,
    ": ",
    paste(missing_groups, collapse = ", "),
    call. = FALSE
  )
}

cat("Cell counts by ", group_by, ":\n", sep = "")
print(table(group_values, useNA = "ifany"))

marker_peaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = group_by,
  useGroups = use_group,
  bgdGroups = background_group,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(marker_peaks, file = marker_peaks_rds)

apt_open_enrich_motifs <- peakAnnoEnrichment(
  seMarker = marker_peaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

pt_open_enrich_motifs <- peakAnnoEnrichment(
  seMarker = marker_peaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

saveRDS(apt_open_enrich_motifs, file = apt_open_enrich_motifs_rds)
saveRDS(pt_open_enrich_motifs, file = pt_open_enrich_motifs_rds)

cat("Differential accessibility analysis complete.\n")
