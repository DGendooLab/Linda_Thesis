library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

set.seed(50)

Samples <- c("s3_mouse", "s8_mouse", "s9_mouse")
RawData <- list()
#RawData[[1]] <- s3_mouse_data
for(i in 1:length(Samples)) {
  RawData[[i]] <- Read10X(data.dir = paste0("Data/", Samples[i]))
  colnames(RawData[[i]]) <- paste(Samples[i], colnames(RawData[[i]]), sep="_")
}
names(RawData) <- Samples

SeuratList <- list()
for(i in 1:length(Samples)){
  SeuratList[[i]] <- CreateSeuratObject(counts=RawData[[i]], project=Samples[i], min.cells=3, min.features=200)
  SeuratList[[i]][["percent.mt"]] <- PercentageFeatureSet(SeuratList[[i]], pattern = "^mt-")
  SeuratList[[i]] <- subset(SeuratList[[i]], subset = nFeature_RNA > 500 & percent.mt < 20)
}

# Loop through SeuratList to get dimensions for each Seurat object
for(i in 1:length(SeuratList)) {
  dims <- dim(SeuratList[[i]])
  cat("Sample:", Samples[i], "- Genes:", dims[1], "Cells:", dims[2], "\n")
}

# pre-processing steps (scaling automatically done in integration)
SeuratList <- lapply(SeuratList, NormalizeData)
SeuratList <- lapply(SeuratList, FindVariableFeatures, selection.method="vst", nfeatures=1500)
names(SeuratList) <- Samples

# Plotting the HVGs
labeled_plots <- list()

# Loop over each Seurat object in SeuratList
for (i in seq_along(SeuratList)) {
  # Generate variable feature plot for the current Seurat object
  plot1 <- VariableFeaturePlot(SeuratList[[i]])
  
  # Extract top variable features for the current Seurat object
  top_variable_features <- head(VariableFeatures(SeuratList[[i]]), 50)
  
  # Label the top variable features on the plot
  plot2 <- LabelPoints(plot = plot1, points = top_variable_features, repel = TRUE, xnudge = 0, ynudge = 0)
  
  # Add sample title to the plot using ggtitle()
  plot2 <- plot2 + ggtitle(paste("Sample:", Samples[i]))  # Add title to the plot
  
  # Store the labeled plot in the list
  labeled_plots[[i]] <- plot2
}


for (plot in labeled_plots) {
  plot.new()   # Create a new plot window
  print(plot)  # Draw the current plot
}


# Plotting PCAs for each sample
PCA_SeuratList <- lapply(SeuratList, ScaleData)
PCA_SeuratList <- lapply(PCA_SeuratList, RunPCA)

pca_plots <- list()

#plotting PCA
for (i in seq_along(PCA_SeuratList)) {
  print(pca_plots[[i]] <- DimPlot(PCA_SeuratList[[i]], reduction = "pca"))
}


plot_elbow <- function(pca_result, sample_name) {
  ElbowPlot(pca_result)
}

# Generate and display the elbow plots for each sample
for (i in seq_along(PCA_SeuratList)) {
  print(plot_elbow(PCA_SeuratList[[i]], sample_name = Samples[i]))
}

#########################################################


#finding pairs of cells in the samples, which helps in creating a single reference
Anchors <- FindIntegrationAnchors(SeuratList)

#creating an integrated dataset using the anchors -- further used for integrated analysis
NormEpi <- IntegrateData(Anchors)

NormEpi <- ScaleData(NormEpi, verbose=FALSE)

dimUsed <- 15
NormEpi <- RunPCA(NormEpi, npcs=dimUsed, verbose=FALSE)
NormEpi <- RunUMAP(NormEpi, dims=1:dimUsed, seed.use=2021)
DimPlot(NormEpi, reduction = "umap", group.by="orig.ident")

NormEpi <- FindNeighbors(NormEpi, dims=1:dimUsed)
NormEpi <- FindClusters(NormEpi, resolution=0.2)

table(Idents(NormEpi))

DimPlot(NormEpi, reduction = "umap", label=TRUE)

# ################# PAIR-WISE COMPARISON ##################
# 
# clusters <- unique(Idents(NormEpi))
# pairwise <- combn(clusters, 2)
# 
# results <- list()
# 
# for(i in 1:ncol(pairwise)) {
#   ident1 <- pairwise[1, i]
#   ident2 <- pairwise[2, i]
#   
#   # Perform differential expression analysis between the two clusters
#   markers <- FindMarkers(NormEpi, ident.1 = ident1, ident.2 = ident2, 
#                          min.pct = 0.2, 
#                          logfc.threshold = 0.25)
#   
#   markers$gene <- rownames(markers)
#   # Add the comparison information to the results
#   comparisons <- pairwise[, i]
#   markers$comparison <- paste(comparisons[1], comparisons[2], sep = '_')
#   
#   # Store the results in the list
#   results[[i]] <- markers
# }
# 
# results.df <- do.call(rbind, results)
# write.csv(results.df, "top100markers1v1-parameters-posneg.csv")
# 
# #Creating an expression heatmap
# top_markers <- results.df %>%
#   group_by(comparison) %>%
#   dplyr::filter(avg_log2FC > 0.8) %>%
#   slice_head(n = 10) %>%
#   ungroup()
# 
# # Extract unique genes from the top markers
# unique_genes <- unique(top_markers$gene)
# 
# # Generating a heatmap
# DoHeatmap(NormEpi, features = unique_genes) + NoLegend()

########### ONE VS REST COMPARISON ##################

NormEpi.markers <- FindAllMarkers(NormEpi,  min.pct = 0.20, logfc.threshold = 0.25)

topMarkers <- split(NormEpi.markers, NormEpi.markers$cluster)
top100 <- lapply(topMarkers, head, n=100)
top100 <- do.call("rbind", top100)
top100
# write.csv(top100,"top100markergenes.csv")