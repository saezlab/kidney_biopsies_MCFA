library(magrittr)
library(ComplexHeatmap)
library(circlize)
seEnrich <- readRDS("data/ATAC/enrich_motifs.rds")

# Adapt function ArchR::plotEnrichHeatmap() without having to have ArchR installed
n = 5 # Number of markers to plot per celltype
cutOff = 20 # Lower bound for -log10(padj) to be considered a marker motif
pMax = Inf
clusterCols = TRUE
binaryClusterRows = TRUE
labelRows = TRUE
rastr = TRUE
transpose = TRUE


mat <- assays(seEnrich)[["mlog10Padj"]]
mat_unfiltered <- assays(seEnrich)[["mlog10Padj"]]
#write.csv(mat_unfiltered, "results/ATAC/motif_enrichment_mlog10Padj.csv")

# List of TF motifs to plot
# Top n per cell type
# Then check for adjusted p value
keep <- lapply(seq_len(ncol(mat)), function(x){
  idx <- head(order(mat[, x], decreasing = TRUE), n)
  rownames(mat)[idx[which(mat[idx,x] > cutOff)]]
}) %>% unlist %>% unique

mat <- mat_unfiltered[keep, ,drop = FALSE]

# Not sure the purpose of this yet
passMat <- lapply(seq_len(nrow(mat)), function(x){
  (mat[x, ] >= 0.9*max(mat[x, ])) * 1
}) %>% Reduce("rbind", .) %>% data.frame
colnames(passMat) <- colnames(mat)

mat[mat > pMax] <- pMax

.rowScale <- function(mat = NULL, min = NULL, max = NULL){
  if(!is.null(min)){
    rMin <- min
  }else{
    rMin <- matrixStats::rowMins(mat)
  }
  if(!is.null(max)){
    rMax <- max
  }else{
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff/rScale
  out <- list(mat=matScale, min=rMin, max=rMax)
  return(out)
}

mat <- .rowScale(as.matrix(mat), min = 0)
if(pMax != 100){
  rownames(mat[[1]]) <- paste0(rownames(mat[[1]]), " (",round(mat$max),")")
  rownames(passMat) <- rownames(mat[[1]])
}

.binarySort <- function(m = NULL, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){
  
  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }
  
  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  return(list(mat = m, hclust = hc))
  
}


mat2 <- mat[[1]] * 100

if(nrow(mat2) > 1 & ncol(mat2) > 1){
  if(binaryClusterRows){
    #cn <- order(colMeans(mat2), decreasing=TRUE)
    bS <- .binarySort(mat2, lmat = passMat[rownames(mat2), colnames(mat2)], clusterCols = clusterCols)
    mat2 <- bS[[1]][,colnames(mat2)]
    clusterRows <- FALSE
    clusterCols <- bS[[2]]
  }else{
    clusterRows <- TRUE
  }
}else{
  clusterCols <- NULL
  clusterRows <- FALSE
}

if(nrow(mat2) > 100){
  borderColor <- FALSE
}else{
  borderColor <- TRUE
}

## Recreate ArchR color palette

# Example definition of ArchRPalettes
ArchRPalettes <- list(
  comet = c("1"="#E6E7E8", "2"="#3A97FF", "3"="#8816A7", "4"="black") 
)

# Function to create a continuous color palette
paletteContinuous <- function(
    set = "comet", 
    n = 256, 
    reverse = FALSE
){
  
  # Validate inputs
  if (!is.character(set)) {
    stop("set must be a character vector or name of a palette")
  }
  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive integer")
  }
  if (!is.logical(reverse)) {
    stop("reverse must be a logical value")
  }
  
  # Get the palette from ArchRPalettes
  pal <- ArchRPalettes[[set]]
  if (is.null(pal)) {
    stop(paste("Palette", set, "not found in ArchRPalettes"))
  }
  
  # Create the continuous color palette
  palOut <- colorRampPalette(pal)(n)
  
  # Reverse the palette if requested
  if (reverse) {
    palOut <- rev(palOut)
  }
  
  return(palOut)
}

# Generate the palette with the specified parameters
pal <- paletteContinuous(set = "comet", n = 100)

mat = as.matrix(mat2)
scale = FALSE
limits = c(0, max(mat2))
color = pal 
clusterCols = FALSE 
clusterRows = FALSE
labelRows = TRUE
useRaster = rastr
fontSizeCols = 6
borderColor = borderColor
customColLabel = seq_len(ncol(mat2))
showRowDendrogram = FALSE
draw = FALSE
name = "Norm. Enrichment -log10(P-adj) [0-Max]"
colAnnoPerRow = 4
fontSizeLabels = 8

#Z-score
if (scale) {
  message("Scaling Matrix..")
  mat <- .rowZscores(mat, limit = FALSE)
  name <- paste0(name," Z-Scores")
}

#Prepare Limits if needed
breaks <- NULL
if (!is.null(limits)) {
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
}else{
  limits <- c(round(min(mat),2), round(max(mat),2))
}

#Scale Values 0 - 1
mat <- (mat - min(limits)) / (max(limits) - min(limits))
breaks <- seq(0, 1, length.out = length(color))
color <- circlize::colorRamp2(breaks, color)

if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
  anno_check_version_rows <- ComplexHeatmap::anno_mark
  anno_check_version_cols <- ComplexHeatmap::anno_mark
}else{
  anno_check_version_rows <- ComplexHeatmap::row_anno_link
  anno_check_version_cols <- ComplexHeatmap::column_anno_link
}

message("Adding Annotations..")
customColLabelIDs <- colnames(mat)[customColLabel]
ht1Anno <- HeatmapAnnotation(
  df = colData,
  #col = colorMap, 
  #show_legend = showLegend,
  show_annotation_name = TRUE,
  gp = gpar(col = "NA"),
  annotation_legend_param =
    list(
      nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
    ),
  foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels))
)

ht1 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(mat) %>% t(),
  name = name,
  col = circlize::colorRamp2(breaks, pal),
  heatmap_legend_param = 
    list(
      at = c(0, 1),
      labels = c(round(min(limits),2), round(max(limits),2)),
      color_bar = "continuous", 
      legend_direction = "horizontal",
      legend_width = unit(3, "cm")
    ), 
  rect_gp = gpar(col = borderColor),
  #Column Options
  show_column_names = TRUE,
  cluster_columns = TRUE, 
  show_column_dend = FALSE,
  clustering_method_columns = "ward.D2",
  column_names_gp = gpar(fontsize = 11), 
  column_names_max_height = unit(100, "mm"),
  #Row Options
  show_row_names = labelRows,
  row_names_gp = gpar(fontsize = 11), 
  cluster_rows = TRUE, 
  show_row_dend = TRUE, 
  clustering_method_rows = "ward.D2",
  split = NULL)
draw(ht1)
