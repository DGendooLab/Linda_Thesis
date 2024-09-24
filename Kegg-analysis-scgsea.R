library(Seurat)
library(biomaRt)
library(GSEABase)
library(fgsea)
library(BiocParallel)
library(dplyr)

NormEpi <- readRDS("/rds/projects/g/gendood-preclinomics/Linda_Thesis/NormEpi.rds")

set.seed(12)


# Initialize the biomaRt object for mouse
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve mapping of gene symbols to Entrez Gene IDs
gene_conversion <- getBM(
  filters = "mgi_symbol",
  attributes = c("mgi_symbol", "entrezgene_id"),
  values = rownames(NormEpi@assays$integrated),
  mart = ensembl
)

# Remove any rows with missing Entrez IDs
gene_conversion <- gene_conversion[!is.na(gene_conversion$entrezgene_id), ]

# Create a mapping from gene symbols to Entrez Gene IDs
symbol_to_entrez <- setNames(gene_conversion$entrezgene_id, gene_conversion$mgi_symbol)

# Map the gene symbols in the Seurat object to Entrez Gene IDs
seurat_genes_entrez <- symbol_to_entrez[rownames(NormEpi@assays$integrated)]

# Remove any genes that could not be mapped
seurat_genes_entrez <- seurat_genes_entrez[!is.na(seurat_genes_entrez)]

# Ensure only mapped genes are included in the Seurat object
mapped_genes <- names(seurat_genes_entrez)

# Subset the Seurat object to keep only the mapped genes
NormEpi_filtered <- subset(NormEpi, features = mapped_genes)

# Update the rownames with the mapped Entrez Gene IDs in the integrated assay
rownames(NormEpi_filtered@assays$integrated@data) <- seurat_genes_entrez
rownames(NormEpi_filtered@assays$integrated@scale.data) <- seurat_genes_entrez

# Load the KEGG Pathway gene sets from the GMT file
gene_sets <- gmtPathways("/rds/projects/g/gendood-preclinomics/Linda_Thesis/Data/mouse_kegg_genesets.gmt")


# Filter the gene sets to include only those with genes in the Seurat object
filtered_gene_sets <- lapply(gene_sets, function(genes) {
  intersect(genes, seurat_genes_entrez)
})

# Remove any empty gene sets
filtered_gene_sets <- filtered_gene_sets[sapply(filtered_gene_sets, length) > 0]

###############removing unwanted pathways################
# List of pathways to exclude
exclude_pathways <- c("mmu00052_Galactose_metabolism", "mmu00220_Arginine_biosynthesis", 
                      "mmu00061_Fatty_acid_biosynthesis", "mmu00230_Purine_metabolism", 
                      "mmu00260_Glycine_serine_and_threonine_metabolism", "mmu00524_Neomycin_kanamycin_and_gentamicin_biosynthesis", 
                      "mmu00500_Starch_and_sucrose_metabolism", 
                      "mmu00450_Selenocompound_metabolism", 
                      "mmu00340_Histidine_metabolism", 
                      "mmu00330_Arginine_and_proline_metabolism", 
                      "mmu00430_Taurine_and_hypotaurine_metabolism", "mmu00630_Glyoxylate_and_dicarboxylate_metabolism", 
                      "mmu00640_Propanoate_metabolism", 
                      "mmu00740_Riboflavin_metabolism", 
                      "mmu00750_Vitamin_B6_metabolism", 
                      "mmu00785_Lipoic_acid_metabolism", 
                      "mmu00670_One_carbon_pool_by_folate", "mmu00910_Nitrogen_metabolism", 
                      "mmu00920_Sulfur_metabolism", 
                      "mmu00900_Terpenoid_backbone_biosynthesis", 
                      "mmu00860_Porphyrin_metabolism", 
                      "mmu00980_Metabolism_of_xenobiotics_by_cytochrome_P450", 
                      "mmu00982_Drug_metabolism", 
                      "mmu00983_Drug_metabolism", "mmu03250_Viral_life_cycle", 
                      "mmu03260_Virion", 
                      "mmu03265_Virion", 
                      "mmu03266_Virion", 
                      "mmu03267_Virion", "mmu04061_Viral_protein_interaction_with_cytokine_and_cytokine_receptor", "mmu04114_Oocyte_meiosis", "mmu04260_Cardiac_muscle_contraction",
                      "mmu04261_Adrenergic_signaling_in_cardiomyocytes", 
                      "mmu04270_Vascular_smooth_muscle_contraction", 
                      "mmu04380_Osteoclast_differentiation", 
                      "mmu04614_Renin-angiotensin_system", "mmu04710_Circadian_rhythm",
                      "mmu04713_Circadian_entrainment",
                      "mmu04714_Thermogenesis",
                      "mmu04720_Long-term_potentiation",
                      "mmu04721_Synaptic_vesicle_cycle",
                      "mmu04722_Neurotrophin_signaling_pathway",
                      "mmu04723_Retrograde_endocannabinoid_signaling",
                      "mmu04724_Glutamatergic_synapse",
                      "mmu04725_Cholinergic_synapse",
                      "mmu04726_Serotonergic_synapse",
                      "mmu04727_GABAergic_synapse",
                      "mmu04728_Dopaminergic_synapse",
                      "mmu04730_Long-term_depression",
                      "mmu04740_Olfactory_transduction",
                      "mmu04742_Taste_transduction",
                      "mmu04744_Phototransduction",
                      "mmu04820_Cytoskeleton_in_muscle_cells",
                      "mmu04913_Ovarian_steroidogenesis",
                      "mmu04914_Progesterone-mediated_oocyte_maturation",
                      "mmu04918_Thyroid_hormone_synthesis",
                      "mmu04920_Adipocytokine_signaling_pathway",
                      "mmu04921_Oxytocin_signaling_pathway",
                      "mmu04922_Glucagon_signaling_pathway",
                      "mmu04923_Regulation_of_lipolysis_in_adipocytes",
                      "mmu04924_Renin_secretion",
                      "mmu04925_Aldosterone_synthesis_and_secretion",
                      "mmu04926_Relaxin_signaling_pathway",
                      "mmu04927_Cortisol_synthesis_and_secretion",
                      "mmu04928_Parathyroid_hormone_synthesis",
                      "mmu04929_GnRH_secretion", "mmu04936_Alcoholic_liver_disease",
                      "mmu04940_Type_I_diabetes_mellitus",
                      "mmu04950_Maturity_onset_diabetes_of_the_young",
                      "mmu04960_Aldosterone-regulated_sodium_reabsorption",
                      "mmu04961_Endocrine_and_other_factor-regulated_calcium_reabsorption",
                      "mmu04962_Vasopressin-regulated_water_reabsorption",
                      "mmu04964_Proximal_tubule_bicarbonate_reclamation",
                      "mmu04966_Collecting_duct_acid_secretion",
                      "mmu04970_Salivary_secretion",
                      "mmu04971_Gastric_acid_secretion",
                      "mmu04972_Pancreatic_secretion",
                      "mmu04973_Carbohydrate_digestion_and_absorption",
                      "mmu04974_Protein_digestion_and_absorption",
                      "mmu04975_Fat_digestion_and_absorption",
                      "mmu04976_Bile_secretion",
                      "mmu04977_Vitamin_digestion_and_absorption",
                      "mmu04978_Mineral_absorption",
                      "mmu05010_Alzheimer_disease","mmu04932_Non-alcoholic_fatty_liver_disease",
                      "mmu04934_Cushing_syndrome","mmu05012_Parkinson_disease",
                      "mmu05014_Amyotrophic_lateral_sclerosis",
                      "mmu05016_Huntington_disease",
                      "mmu05017_Spinocerebellar_ataxia",
                      "mmu05020_Prion_disease",
                      "mmu05022_Pathways_of_neurodegeneration",
                      "mmu05030_Cocaine_addiction",
                      "mmu05031_Amphetamine_addiction",
                      "mmu05032_Morphine_addiction",
                      "mmu05033_Nicotine_addiction",
                      "mmu05034_Alcoholism",
                      "mmu05100_Bacterial_invasion_of_epithelial_cells",
                      "mmu05132_Salmonella_infection",
                      "mmu05133_Pertussis",
                      "mmu05134_Legionellosis",
                      "mmu05135_Yersinia_infection",
                      "mmu05140_Leishmaniasis",
                      "mmu05142_Chagas_disease","mmu05143_African_trypanosomiasis",
                      "mmu05144_Malaria",
                      "mmu05145_Toxoplasmosis",
                      "mmu05146_Amoebiasis",
                      "mmu05150_Staphylococcus_aureus_infection",
                      "mmu05152_Tuberculosis",
                      "mmu05160_Hepatitis_C",
                      "mmu05161_Hepatitis_B",
                      "mmu05162_Measles",
                      "mmu05164_Influenza_A",
                      "mmu05165_Human_papillomavirus_infection",
                      "mmu05166_Human_T-cell_leukemia_virus_1_infection",
                      "mmu05167_Kaposi_sarcoma-associated_herpesvirus_infection",
                      "mmu05168_Herpes_simplex_virus_1_infection",
                      "mmu05169_Epstein-Barr_virus_infection",
                      "mmu05170_Human_immunodeficiency_virus_1_infection",
                      "mmu05171_Coronavirus_disease",
                      "mmu05163_Human_cytomegalovirus_infection","mmu05310_Asthma",
                      "mmu05320_Autoimmune_thyroid_disease",
                      "mmu05321_Inflammatory_bowel_disease",
                      "mmu05322_Systemic_lupus_erythematosus",
                      "mmu05323_Rheumatoid_arthritis",
                      "mmu05330_Allograft_rejection",
                      "mmu05332_Graft-versus-host_disease",
                      "mmu05340_Primary_immunodeficiency",
                      "mmu05410_Hypertrophic_cardiomyopathy",
                      "mmu05412_Arrhythmogenic_right_ventricular_cardiomyopathy",
                      "mmu05414_Dilated_cardiomyopathy",
                      "mmu05415_Diabetic_cardiomyopathy",
                      "mmu05416_Viral_myocarditis",
                      "mmu05417_Lipid_and_atherosclerosis",
                      "mmu05418_Fluid_shear_stress_and_atherosclerosis"
)

######################################

# Assuming filtered_gene_sets is your list of gene sets
filtered_gene_sets <- filtered_gene_sets[!names(filtered_gene_sets) %in% exclude_pathways]


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


# Run fgsea enrichment analysis on a subset with parallel processing
fgsea_results_parallel <- perform_fgsea(NormEpi_filtered, filtered_gene_sets)

saveRDS(fgsea_results_parallel, file = "/rds/projects/g/gendood-preclinomics/Linda_Thesis/fgsea_results_v5.rds")

################################


# Convert fgsea results to a data frame
fgsea_df <- do.call(rbind, lapply(fgsea_results_parallel, function(res) {
  as.data.frame(res)
}))


# Add rownames as a new column
fgsea_df <- cbind(cell = rownames(fgsea_df), fgsea_df)

# Ensure the leadingEdge column is properly formatted
fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, function(x) paste(unlist(x), collapse = ","))

rownames(fgsea_df) <- NULL

# Write to CSV with the new column
write.csv(fgsea_df, file = "/rds/projects/g/gendood-preclinomics/Linda_Thesis/fgsea_results_v5_scgsea_kegg.csv", row.names = FALSE)

