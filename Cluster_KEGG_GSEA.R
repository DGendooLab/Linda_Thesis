# Load necessary libraries
library(Seurat)
library(dplyr)
library(fgsea)
library(GSEABase)
library(ggplot2)
library(pheatmap)
library(biomaRt)

# Set seed for reproducibility
set.seed(12)

# Find markers for all clusters
NormEpi.markers <- readRDS("/rds/projects/g/gendood-preclinomics/Linda_Thesis/NormEpi_markers.rds")

# Function to map gene names to Entrez IDs
get_entrez_ids <- function(gene_names) {
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  gene_map <- getBM(attributes = c("mgi_symbol", "entrezgene_id"),
                    filters = "mgi_symbol",
                    values = gene_names,
                    mart = mart)
  return(gene_map)
}

# Get unique gene names from the markers
unique_genes <- unique(NormEpi.markers$gene)

# Map gene names to Entrez IDs
gene_map <- get_entrez_ids(unique_genes)

# Merge Entrez IDs into the markers data frame
NormEpi.markers <- left_join(NormEpi.markers, gene_map, by = c("gene" = "mgi_symbol"))

# Filter markers with Entrez IDs
NormEpi.markers <- NormEpi.markers %>% filter(!is.na(entrezgene_id))

# Split the data after the Entrez ID has been added
topMarkers <- split(NormEpi.markers, NormEpi.markers$cluster)

# Create ranked gene lists for fgsea
ranks_list <- list()

for (cluster in names(topMarkers)) {
  print(paste("Processing cluster:", cluster))  # Debugging statement
  
  # Ensure entrezgene_id is present
  if ("entrezgene_id" %in% colnames(topMarkers[[cluster]])) {
    ranked_genes <- tryCatch({
      topMarkers[[cluster]] %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        dplyr::select(entrezgene_id, avg_log2FC) %>%
        dplyr::rename(ranking = avg_log2FC)
    }, error = function(e) {
      print(paste("Error in ranking genes for cluster:", cluster, " - ", e$message))
      NULL
    })
    
    # If ranking was successful, convert to named vector for fgsea
    if (!is.null(ranked_genes)) {
      print(head(ranked_genes))  # Print the ranked genes for debugging
      ranks <- setNames(ranked_genes$ranking, ranked_genes$entrezgene_id)
      ranks_list[[cluster]] <- ranks
    }
  } else {
    print(paste("Error: entrezgene_id not found in cluster:", cluster))
  }
}

# Load pathways from MSigDB or other relevant databases
pathways <- gmtPathways("/rds/homes/l/lsm221/thesis/Data/mouse_kegg_genesets.gmt")

###############removing unwanted pathways################
# List of pathways to exclude
exclude_pathways <- c("mmu01522_Endocrine_resistance", "mmu04070_Phosphatidylinositol_signaling_system", "mmu04071_Sphingolipid_signaling_pathway", "mmu04350_TGF-beta_signaling_pathway", "mmu04931_Insulin_resistance", "mmu05223_Non-small_cell_lung_cancer", "mmu01230_Biosynthesis_of_amino_acids", 
                      "mmu01521_EGFR_tyrosine_kinase_inhibitor_resistance", "mmu04657_IL-17_signaling_pathway", "mmu04912_GnRH_signaling_pathway", "mmu04916_Melanogenesis", "mmu04917_Prolactin_signaling_pathway", "mmu05212_Pancreatic_cancer", "mmu05218_Melanoma", "mmu05220_Chronic_myeloid_leukemia", 
                      "mmu05330_Allograft_rejection", "mmu05332_Graft-versus-host_disease", "mmu04142_Lysosome", "mmu04150_mTOR_signaling_pathway", "mmu04520_Adherens_junction", "mmu04664_Fc_epsilon_RI_signaling_pathway", "mmu04910_Insulin_signaling_pathway", "mmu05235_PD-L1_expression_and_PD-1_checkpoint_pathway_in_cancer", 
                      "mmu01200_Carbon_metabolism", "mmu04141_Protein_processing_in_endoplasmic_reticulum", "mmu04217_Necroptosis", "mmu05211_Renal_cell_carcinoma", "mmu05213_Endometrial_cancer", "mmu05230_Central_carbon_metabolism_in_cancer", "mmu00010_Glycolysis_/_Gluconeogenesis", 
                      "mmu00564_Glycerophospholipid_metabolism", "mmu01524_Platinum_drug_resistance", "mmu04140_Autophagy", "mmu04211_Longevity_regulating_pathway", "mmu04623_Cytosolic_DNA-sensing_pathway", "mmu04930_Type_II_diabetes_mellitus", "mmu05216_Thyroid_cancer", 
                      "mmu00590_Arachidonic_acid_metabolism", "mmu03320_PPAR_signaling_pathway", "mmu04120_Ubiquitin_mediated_proteolysis", "mmu04213_Longevity_regulating_pathway", "mmu04330_Notch_signaling_pathway", "mmu04370_VEGF_signaling_pathway", "mmu00534_Glycosaminoglycan_biosynthesis", 
                      "mmu04622_RIG-I-like_receptor_signaling_pathway", "mmu04672_Intestinal_immune_network_for_IgA_production", "mmu05217_Basal_cell_carcinoma", "mmu00561_Glycerolipid_metabolism", "mmu00051_Fructose_and_mannose_metabolism", "mmu00240_Pyrimidine_metabolism", "mmu00480_Glutathione_metabolism", 
                      "mmu00562_Inositol_phosphate_metabolism", "mmu00565_Ether_lipid_metabolism", "mmu01232_Nucleotide_metabolism", "mmu04216_Ferroptosis", "mmu05204_Chemical_carcinogenesis", "mmu00030_Pentose_phosphate_pathway", "mmu00514_Other_types_of_O-glycan_biosynthesis", "mmu00591_Linoleic_acid_metabolism", 
                      "mmu02010_ABC_transporters", "mmu04137_Mitophagy", "mmu00140_Steroid_hormone_biosynthesis", "mmu00190_Oxidative_phosphorylation", "mmu00512_Mucin_type_O-glycan_biosynthesis", "mmu00513_Various_types_of_N-glycan_biosynthesis", "mmu00532_Glycosaminoglycan_biosynthesis", "mmu00600_Sphingolipid_metabolism", 
                      "mmu00601_Glycosphingolipid_biosynthesis", "mmu00830_Retinol_metabolism", "mmu03460_Fanconi_anemia_pathway", "mmu04215_Apoptosis", "mmu04340_Hedgehog_signaling_pathway", "mmu05219_Bladder_cancer", "mmu00250_Alanine_aspartate_and_glutamate_metabolism", "mmu00270_Cysteine_and_methionine_metabolism", 
                      "mmu00280_Valine_leucine_and_isoleucine_degradation", "mmu00310_Lysine_degradation", "mmu00350_Tyrosine_metabolism", "mmu00510_N-Glycan_biosynthesis", "mmu00592_alpha-Linolenic_acid_metabolism", "mmu00760_Nicotinate_and_nicotinamide_metabolism", "mmu01240_Biosynthesis_of_cofactors", 
                      "mmu03083_Polycomb_repressive_complex", "mmu04392_Hippo_signaling_pathway", "mmu00380_Tryptophan_metabolism", "mmu00410_beta-Alanine_metabolism", "mmu00520_Amino_sugar_and_nucleotide_sugar_metabolism", "mmu00620_Pyruvate_metabolism", "mmu00650_Butanoate_metabolism", 
                      "mmu01250_Biosynthesis_of_nucleotide_sugars", "mmu01523_Antifolate_resistance", "mmu03013_Nucleocytoplasmic_transport", "mmu03018_RNA_degradation", "mmu03040_Spliceosome", "mmu03050_Proteasome", "mmu03440_Homologous_recombination", "mmu04146_Peroxisome", "mmu00100_Steroid_biosynthesis", 
                      "mmu00511_Other_glycan_degradation", "mmu00515_Mannose_type_O-glycan_biosynthesis")

######################################

# Filter out the unwanted pathways
pathways <- pathways[!names(pathways) %in% exclude_pathways]

# Add small random noise to the ranking to break ties
set.seed(123)
ranks_list <- lapply(ranks_list, function(ranks) {
  ranks <- ranks + rnorm(length(ranks), mean = 0, sd = 1e-6)
  ranks
})

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
  # write.csv(fgsea_res_df, paste0("fgsea_res_cluster_", cluster, ".csv"))
  
  # Save the result into the list
  fgsea_results[[cluster]] <- fgsea_res_df
}

# Initialize an empty list to store NES scores
nes_list <- list()

# Extract NES scores and create a matrix
for (cluster in names(fgsea_results)) {
  fgsea_res_df <- fgsea_results[[cluster]]
  
  # Proceed only if it is a data.frame
  if (is.data.frame(fgsea_res_df)) {
    # Extract pathway and NES columns
    nes_scores <- dplyr::select(fgsea_res_df, pathway, NES)
    
    # Add cluster identifier to pathway names
    colnames(nes_scores)[2] <- paste("Cluster", cluster, sep = " ")
    
    # Convert to data frame (optional, already a data.frame)
    nes_df <- as.data.frame(nes_scores)
    
    # Append to the list
    nes_list[[cluster]] <- nes_df
  } 
}


# Merge all NES data frames by the 'pathway' column
nes_matrix <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), nes_list)

# Replace NA values with 0 (or any other appropriate value)
nes_matrix[is.na(nes_matrix)] <- 0

# Set row names to pathways and remove the pathway column
rownames(nes_matrix) <- nes_matrix$pathway
nes_matrix <- nes_matrix %>% dplyr::select(-pathway)

# Remove the "mmuXXXX_" prefix from row names
cleaned_rownames <- sub("^mmu[0-9]+_", "", rownames(nes_matrix))

# Make row names unique
rownames(nes_matrix) <- make.unique(cleaned_rownames)


# # Select top pathways based on a threshold or top N pathways
top_n <- 20  # You can adjust this number
top_pathways <- nes_matrix[order(rowMeans(nes_matrix), decreasing = TRUE)[1:top_n], ]

# Create the heatmap
outputpath <- ("/rds/projects/g/gendood-preclinomics/Linda_Thesis/plots/heatmap_GSEA_KEGG_allpathways.pdf")
pdf(outputpath)
pheatmap(as.matrix(nes_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row", 
         fontsize_row = 6,  # Reduce the text size for row names
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Top NES Scores Heatmap")
dev.off()
