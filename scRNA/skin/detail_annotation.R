suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(SingleR)
  library(SingleCellExperiment)
  library(dplyr)
})

# ==============================================================================
# 1. SETUP & LOAD DATA
# ==============================================================================
out_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/haniffa_analysis/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading mathematically cleaned fibroblast object...")
FB.subgroup <- readRDS("/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_subset_processed.rds")

# Ensure Default Assay is RNA
DefaultAssay(FB.subgroup) <- "RNA"

# ==============================================================================
# 2. THE DICTIONARIES (Separated for Mathematical Accuracy)
# ==============================================================================

# BIOLOGICAL DICTIONARY: Used ONLY for SingleR Auto-Annotation. No stress markers allowed!
biological_dictionary <- list(
  "F1_Superficial_Papillary" = c("COL13A1", "COL18A1", "COL23A1", "APCDD1", "WIF1", "NKD2", "CRABP1", "CYP26B1", "WNT5A", "COMP", "AXIN2", "RSPO1", "SFRP2", "IGFBP2", "TNFRSF21", "PDE4B", "GUCY1A1"),
  "F2_Universal_Reticular" = c("PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PCOLCE2", "CTHRC1", "SLPI", "CD70", "LGR5"),
  "F2_F3_Perivascular" = c("PPARG", "CXCL12", "APOC1", "C7", "PLA2G2A", "APOE", "EFEMP1", "MYOC", "GDF10"),
  "F3_FRC_like" = c("CCL19", "CXCL12", "CH25H", "IL33", "IL15", "TNFSF13B", "VCAM1", "CD74", "HLA-DRA", "CXCL9", "ADAMDEC1", "IRF8", "HLA-DRB1"),
  "F4_HairFollicle_Associated" = c("ASPN", "COL11A1", "DPEP1", "MKX", "TNMD", "CORIN", "HHIP", "RSPO3", "LEF1", "MEF2C", "MYL4", "TNN", "COCH", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "BMP7", "INHBA", "PTCH1"),
  "F5_Schwann_like" = c("SCN7A", "FMO2", "FGFBP2", "OLFML2A", "NGFR", "ITGA6", "EBF2", "RAMP1", "PEAR1", "RELN", "PLEKHA6", "IGFBP2", "SFRP1", "CDH19", "CLDN1"),
  "Myofibroblast_Disease_Signature" = c("ACTA2", "COL3A1", "COL5A1", "COL8A1", "POSTN", "CTHRC1", "LRRC15", "SFRP4", "ASPN", "RUNX2", "SCX"),
  "F6_Inflammatory_Myofibroblast" = c("IL11", "IL24", "IL7R", "CXCL5", "CXCL8", "CXCL13", "CCL11", "MMP1", "MMP3", "CSF3", "TDO2", "WWC1", "CHI3L1", "STAT4", "CCL5", "CCL3", "FAM167A"),
  "F7_Myofibroblast" = c("PIEZO2", "COL1A1", "COL3A1", "POSTN", "WNT2", "COL10A1", "LAMP5", "NRG1", "OGN", "TAGLN", "KIF26B", "ZNF469", "SULF1", "ADAM12", "CDH2", "LRRC17", "KCNMA1", "ADAM19", "CREB3L1", "CCN4", "FABP5", "C1QTNF3", "CADM1"),
  "F8_Fascia_like_Myofibroblast" = c("ACAN", "THBS4", "ITGA10", "FGF18", "PRG4", "CRTAC1"),
  "Prenatal_LTo_like" = c("CCL21", "CXCL13", "MADCAM1", "FDCSP", "TNFSF11")
)

# VALIDATION DICTIONARY: Used for generating plots. Re-attaches the Stress Markers so you can see them.
validation_dictionary <- biological_dictionary
validation_dictionary[["QC_Stress_Markers"]] <- c("MT2A", "MT1M", "MT1X", "HSP90AA1", "JUNB", "GADD45B", "IER3")

# ==============================================================================
# 3. THE F1-F8 SINGLE-R AUTO-ANNOTATION ENGINE
# ==============================================================================
message("Running SingleR Auto-Annotation for F1-F8 Taxonomy...")

build_custom_singler_reference <- function(marker_dictionary, gene_universe) {
  ref_mat <- sapply(names(marker_dictionary), function(label) {
    genes <- intersect(marker_dictionary[[label]], gene_universe)
    vec <- setNames(rep(0, length(gene_universe)), gene_universe)
    vec[genes] <- 1
    vec
  })
  ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(ref_mat)))
  colData(ref_sce)$label <- colnames(ref_mat)
  return(ref_sce)
}

gene_universe <- rownames(FB.subgroup)
ref_sce <- build_custom_singler_reference(biological_dictionary, gene_universe) # Uses pure biology!
test_sce <- as.SingleCellExperiment(FB.subgroup, assay = "RNA")

# Run SingleR
singler_res <- SingleR(
  test = test_sce,
  ref = ref_sce,
  labels = ref_sce$label,
  clusters = FB.subgroup$seurat_clusters, 
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  prune = FALSE 
)

# Map labels back
singler_df <- as.data.frame(singler_res)
label_map <- setNames(singler_df$labels, rownames(singler_df))
FB.subgroup$SingleR_label <- unname(label_map[as.character(FB.subgroup$seurat_clusters)])

# Set the active identity
Idents(FB.subgroup) <- "SingleR_label"

message("\n--- Biological Clusters Assigned ---")
print(table(Idents(FB.subgroup)))
message("------------------------------------\n")

p_annotated_umap <- DimPlot(FB.subgroup, label = TRUE, repel = TRUE) + 
  ggtitle("SingleR Auto-Annotated F1-F8 Fibroblasts")
ggsave(file.path(out_dir, "UMAP_Annotated_F1_F8.pdf"), plot = p_annotated_umap, width = 10, height = 7)
# ==============================================================================
# 4. HANIFFA F1-F8 MARKER VISUALIZATION (Automated Loop)
# ==============================================================================
message("Generating Crash-Proof Validation Plots...")

for (category_name in names(validation_dictionary)) {
  genes_to_plot <- validation_dictionary[[category_name]]
  valid_genes <- intersect(genes_to_plot, rownames(FB.subgroup))
  
  if (length(valid_genes) == 0) next
  
  # FIX: Filter out genes that have literally 0 expression across ALL cells
  # This prevents the exact 'FDCSP' crash you just experienced.
  exp_data <- GetAssayData(FB.subgroup, assay = "RNA", slot = "data")
  valid_genes <- valid_genes[rowSums(exp_data[valid_genes, , drop = FALSE]) > 0]
  
  if (length(valid_genes) == 0) {
    message(paste("   -> Skipping", category_name, "(All genes have 0 expression)"))
    next
  }
  
  # DotPlot wrapped in tryCatch so a math error won't stop the whole pipeline
  tryCatch({
    p_dot <- DotPlot(FB.subgroup, features = valid_genes) + 
      RotatedAxis() + 
      ggtitle(paste("Fibroblast -", category_name))
    ggsave(file.path(out_dir, paste0("DotPlot_", category_name, ".pdf")), plot = p_dot, width = max(6, length(valid_genes)*0.5 + 2), height = 6)
  }, error = function(e) message(paste("   -> Error generating DotPlot for", category_name)))
  
  # VlnPlot wrapped in tryCatch
  tryCatch({
    p_vln <- VlnPlot(FB.subgroup, features = valid_genes, stack = TRUE, flip = TRUE) +
      theme(legend.position = "none") +
      ggtitle(paste("Fibroblast -", category_name))
    ggsave(file.path(out_dir, paste0("VlnPlot_", category_name, ".pdf")), plot = p_vln, width = 8, height = max(6, length(valid_genes)*1.5))
  }, error = function(e) message(paste("   -> Error generating VlnPlot for", category_name)))
}

# ==============================================================================
# 5. PUBLICATION MASTER DOTPLOT
# ==============================================================================
message("Generating Publication-Ready Master DotPlot...")

paper_signature_genes <- c(
  "APCDD1", "COL18A1", "CRABP1", "CYP26B1", "WNT5A",
  "PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PPARG", "CXCL12",
  "CCL19", "CD74", "CXCL9", "DPEP1", "TNN", "CORIN", "RAMP1", "NGFR",
  "IL11", "CXCL8", "MMP1", "ACTA2", "COL11A1", "COMP", "LRRC15", "ACAN", "ITGA10", "THBS4"
)
valid_paper_genes <- intersect(paper_signature_genes, rownames(FB.subgroup))

expected_order <- c(
  "F1_Superficial_Papillary", "F2_Universal_Reticular", "F2_F3_Perivascular",
  "F3_FRC_like", "F4_HairFollicle_Associated", "F5_Schwann_like", 
  "F6_Inflammatory_Myofibroblast", "F7_Myofibroblast", "F8_Fascia_like_Myofibroblast", 
  "Myofibroblast_Disease_Signature", "Prenatal_LTo_like", "QC_Stress_Markers", "Unassigned"
)

current_levels <- unique(as.character(Idents(FB.subgroup)))
ordered_levels <- intersect(expected_order, current_levels)
Idents(FB.subgroup) <- factor(Idents(FB.subgroup), levels = rev(ordered_levels))

master_dotplot <- DotPlot(FB.subgroup, features = valid_paper_genes, dot.scale = 6) +
  theme_minimal() +
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
  labs(
    title = "Fibroblast Subpopulation Signatures in PMH",
    x = "Key Marker Genes",
    y = "Identified Subclusters"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Publication_Master_DotPlot.pdf"), plot = master_dotplot, width = 14, height = 8)

message("Saving final annotated RDS...")
saveRDS(FB.subgroup, file.path(dirname(out_dir), "processed", "fibroblast_annotated_final.rds"))

message("=== Annotation Complete! ===")