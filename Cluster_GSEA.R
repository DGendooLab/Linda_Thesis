# Load necessary libraries
library(Seurat)
library(dplyr)
library(fgsea)
library(GSEABase)
library(ggplot2)
library(pheatmap)


NormEpi.markers <- FindAllMarkers(NormEpi,  min.pct = 0.20, logfc.threshold = 0.25)
topMarkers <- split(NormEpi.markers, NormEpi.markers$cluster)
# topMarkers_onevall <- do.call("rbind", topMarkers)

top_markers <- NormEpi.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.8) %>%
  slice_head(n = 10) %>%
  ungroup()

ranks_list <- list()


for (cluster in names(topMarkers)) {
  print(paste("Processing cluster:", cluster))  # Debugging statement
  
  # Rank genes and handle any potential errors
  ranked_genes <- tryCatch({
    topMarkers[[cluster]] %>%
      dplyr::arrange(desc(avg_log2FC)) %>%  # Arrange genes by log fold-change
      dplyr::select(gene, avg_log2FC) %>%  # Select the gene and the metric to rank by
      dplyr::rename(ranking = avg_log2FC)  # Rename the metric column to 'ranking'
  }, error = function(e) {
    print(paste("Error in ranking genes for cluster:", cluster, " - ", e$message))
    NULL
  })
  
  # If ranking was successful, convert to named vector for fgsea
  if (!is.null(ranked_genes)) {
    ranks <- setNames(ranked_genes$ranking, ranked_genes$gene)
    ranks_list[[cluster]] <- ranks
  }
}

####################################

# Load pathways from MSigDB or other relevant databases
pathways <- gmtPathways("/rds/projects/g/gendood-preclinomics/Linda_Thesis/Data/mh.all.v2023.2.Mm.symbols.gmt")

# Create a list to store fgsea results for each cluster
fgsea_results <- list()

# Run GSEA for each cluster using fgseaMultilevel
for (cluster in names(ranks_list)) {
  fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks_list[[cluster]])
  
  # Convert fgsea results to data frame and remove any list columns
  fgsea_res_df <- as.data.frame(fgsea_res)
  fgsea_res_df <- lapply(fgsea_res_df, function(x) {
    if (is.list(x)) {
      return(sapply(x, paste, collapse = ", "))
    } else {
      return(x)
    }
  })
  fgsea_res_df <- as.data.frame(fgsea_res_df)
  
  # Save results for each cluster
  write.csv(fgsea_res_df, paste0("fgsea_res_cluster_", cluster, ".csv"))
  
  # Save the result into the list
  fgsea_results[[cluster]] <- fgsea_res_df
}

# Initialize an empty list to store NES scores
nes_list <- list()

# Extract NES scores and create a matrix
for (cluster in names(fgsea_results)) {
  fgsea_res_df <- fgsea_results[[cluster]]
  
  # Extract pathway and NES columns
  nes_scores <- fgsea_res_df %>% select(pathway, NES)
  
  # Add cluster identifier to pathway names
  colnames(nes_scores)[2] <- paste("NES", cluster, sep = "_")
  
  # Convert to data frame
  nes_df <- as.data.frame(nes_scores)
  
  # Append to the list
  nes_list[[cluster]] <- nes_df
}

# Merge all NES data frames by the 'pathway' column
nes_matrix <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), nes_list)

# Replace NA values with 0 (or any other appropriate value)
nes_matrix[is.na(nes_matrix)] <- 0

# Set row names to pathways and remove the pathway column
rownames(nes_matrix) <- nes_matrix$pathway
nes_matrix <- nes_matrix %>% select(-pathway)

# Create the heatmap
outputpath <- ("/rds/projects/g/gendood-preclinomics/Linda_Thesis/plots/heatmap_GSEA.pdf")
pdf(outputpath)
pheatmap(as.matrix(nes_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "NES Scores Heatmap")
dev.off()
