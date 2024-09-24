library(dplyr)
library(tidyr)
library(pheatmap)
library(ComplexHeatmap)

# Step 1: Load the Data
fgsea_df <- read.csv("/rds/homes/l/lsm221/thesis/fgsea_results_v4_scgsea_kegg.csv")

# Step 2: Split the cell identifier to separate cell ID and pathway number and Filter Data
fgsea_df_results <- fgsea_df %>%
  mutate(cell_id = sub("\\..*", "", cell),
         pathway_num = sub(".*\\.", "", cell))

# fgsea_df_results$pathway <- sub("^[^_]*_", "", fgsea_df_results$pathway)

filtered_df <- fgsea_df_results %>% filter(padj <= 0.05)

# Step 3: Reshape Data
# Choose the metric for the heatmap, e.g., NES (Normalized Enrichment Score)
heatmap_data <- filtered_df %>%
  select(cell_id, pathway, NES) %>%
  spread(key = pathway, value = NES)

write.csv(heatmap_data,"/rds/homes/l/lsm221/thesis/cell_NES_heatmap.csv")

# Convert to matrix and handle NAs (e.g., replace with 0 or a specific value)
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$cell_id
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Step 4: Generate the Heatmap
# Using pheatmap
output_path <- ("/rds/homes/l/lsm221/thesis/cell_NES_heatmap_TEST.pdf")
pdf(output_path)
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = TRUE)
dev.off()
# Or using ComplexHeatmap
Heatmap(heatmap_matrix,
        name = "NES",
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8))
