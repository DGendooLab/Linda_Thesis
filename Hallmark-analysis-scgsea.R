library(Seurat)
library(GSEABase)
library(fgsea)
library(dplyr)

NormEpi <- readRDS("/rds/homes/l/lsm221/thesis/NormEpi.rds")

set.seed(12)

# Load the Hallmark gene sets
gene_sets <- gmtPathways("/rds/homes/l/lsm221/thesis/Data/mh.all.v2023.2.Mm.symbols.gmt")

# Filter the gene sets to include only those with genes in the Seurat object
seurat_genes <- rownames(NormEpi@assays$integrated)

filtered_gene_sets <- lapply(gene_sets, function(genes) {
  intersect(genes, seurat_genes)
})

# Remove any empty gene sets
filtered_gene_sets <- filtered_gene_sets[sapply(filtered_gene_sets, length) > 0]

# Check the lengths of the filtered gene sets
gene_set_lengths <- sapply(filtered_gene_sets, length)
summary(gene_set_lengths)

############ without parallel processing ############ 
# Function to perform fgsea enrichment analysis
perform_fgsea <- function(seurat_obj, gene_sets) {
  data <- as.matrix(seurat_obj@assays$integrated@scale.data)
  fgsea_results <- list()
  
  for (i in 1:ncol(data)) {
    message(paste("Processing column:", i, "of", ncol(data)))
    ranks <- rank(-data[, i])
    names(ranks) <- rownames(data)
    fgsea_res <- fgsea(pathways = gene_sets, stats = ranks, minSize = 1, maxSize = Inf)
    fgsea_results[[colnames(data)[i]]] <- fgsea_res
  }
  
  return(fgsea_results)
}

# Run fgsea enrichment analysis on the subset with parallel processing
fgsea_results <- perform_fgsea(NormEpi, filtered_gene_sets)

# Save results to an RDS file
saveRDS(fgsea_results, file = "/rds/homes/l/lsm221/thesis/fgsea_results_hallmark_v5.rds")

################################

# Convert fgsea results to a data frame
fgsea_df <- do.call(rbind, lapply(fgsea_results, function(res) {
  as.data.frame(res)
}))

# Add rownames as a new column
fgsea_df <- cbind(cell = rownames(fgsea_df), fgsea_df)

# Ensure the leadingEdge column is properly formatted
fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, function(x) paste(unlist(x), collapse = ","))

rownames(fgsea_df) <- NULL

# Write to CSV with the new column
write.csv(fgsea_df, file = "/rds/homes/l/lsm221/thesis/fgsea_results_hallmark_v5_scgsea.csv", row.names = FALSE)
