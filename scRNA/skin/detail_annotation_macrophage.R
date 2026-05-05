suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(SingleR)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(42)

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS (SNAKEMAKE INTEGRATION)
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the subset RDS object"),
  make_option(c("-c", "--celltype"), type="character", help="The cell type being processed (e.g., 'fibroblast' or 'macrophage')"),
  make_option(c("-o", "--outdir"), type="character", help="Base subset_cluster output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Create highly organized subdirectories dynamically
base_out_dir <- file.path(opt$outdir, opt$celltype, "haniffa_analysis")
dirs <- list(
  umaps_global  = file.path(base_out_dir, "1_UMAPs", "Global"),
  umaps_sample  = file.path(base_out_dir, "1_UMAPs", "Per_Sample"),
  val_vln       = file.path(base_out_dir, "2_Validation_Plots", "VlnPlots"),
  val_dot       = file.path(base_out_dir, "2_Validation_Plots", "DotPlots"),
  master_sum    = file.path(base_out_dir, "3_Summary_Plots"),
  mucin         = file.path(base_out_dir, "4_Mucin_ECM"),
  proportions   = file.path(base_out_dir, "5_Proportions")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

message(paste("Loading mathematically cleaned", toupper(opt$celltype), "object..."))
seu_obj <- readRDS(opt$input)
DefaultAssay(seu_obj) <- "RNA"

save_dual_format <- function(plot_obj, file_path_no_ext, w = 10, h = 8) {
  pdf(paste0(file_path_no_ext, ".pdf"), width = w, height = h)
  print(plot_obj)
  dev.off()
  png(paste0(file_path_no_ext, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot_obj)
  dev.off()
}

# ==============================================================================
# 2. DICTIONARIES & COLORS (FIBROBLAST & MACROPHAGE)
# ==============================================================================
biological_dictionary <- list(
  "F1_Superficial_Healthy" = c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "COMP", "NKD2", "RSPO1", "AXIN2"),
  "F1_Superficial_Disease" = c("CRABP1", "CYP26B1", "TNFRSF21", "CXCL1", "WNT5A", "COL18A1", "COL23A1", "COL13A1", "NKD2", "AXIN2", "RSPO1"),
  "F2_Reticular_Healthy" = c("CD34", "PI16", "MFAP5", "DPP4", "PCOLCE2", "LGR5", "SLPI", "CD70", "CTHRC1"),
  "F2_Reticular_Disease" = c("PI16", "DPP4", "PCOLCE2", "MFAP5", "CD70", "LGR5"),
  "F2_F3_Perivascular_Healthy" = c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10"),
  "F2_F3_Perivascular_Disease" = c("CXCL12", "APOE", "C7", "PLA2G2A", "EFEMP1", "GDF10", "MYOC"),
  "F3_FRC_like_Healthy" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "IRF8", "IL15", "VCAM1", "HLA-DRB1", "HLA-DRA"),
  "F3_FRC_like_Disease" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "HLA-DRA", "IRF8", "COX4I2", "RBP5", "ADAMDEC1", "CXCL9", "CXCL10", "APOE", "CXCL12"),
  "F4_HF_DPEP_Healthy" = c("DPEP1", "MYL4", "MEF2C", "COL11A1"),
  "F4_HF_DPEP_Disease" = c("MEF2C", "MYL4", "COL11A1", "POSTN", "DPEP1"),
  "F4_HF_TNN_Healthy" = c("TNN", "COCH", "TNMD", "MKX", "NRG3", "SLITRK6"),
  "F4_HF_TNN_Disease" = c("COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD"),
  "F4_HF_DP_Healthy" = c("CORIN", "HHIP", "BMP7", "WNT5A", "LEF1"),
  "F4_HF_DP_Disease" = c("CRABP1", "COL24A1", "RSPO4", "RSPO3", "BMP7", "WNT5A", "LEF1", "SOX18", "HHIP"),
  "F5_Schwann_NGFR_Healthy" = c("NGFR", "ITGA6", "SCN7A", "CDH19", "CLDN1", "SFRP4", "TENM2"),
  "F5_Schwann_NGFR_Disease" = c("NGFR", "TM4SF1", "SFRP4", "ANGPTL7", "ITGA6", "CDH19", "CLDN1", "EBF2", "OLFM2", "SCN7A"),
  "F5_Schwann_RAMP1_Healthy" = c("RAMP1", "RELN", "PLEKHA6", "IGFBP2", "FGFBP2", "SCN7A"),
  "F5_Schwann_RAMP1_Disease" = c("RAMP1", "IGFBP2", "RELN", "COL26A1", "PLEKHA6", "FMO2", "FGFBP2"),
  "F6_Inflammatory_Myo_Disease" = c("IL11", "IL24", "CXCL5", "CXCL6", "CXCL8", "MMP9", "WNT2", "COL10A1", "MMP1", "MMP3"),
  "F7_Myofibroblast" = c("ACTA2", "TAGLN", "CTHRC1", "RUNX2", "KIF26B", "SULF1", "ADAM12", "COL8A1", "LRRC15", "CCN4", "ASPN", "POSTN", "TNC", "COL3A1", "WNT2", "COL10A1"),
  "F8_Fascia_like_Myo_Disease" = c("ACAN", "ITGA10", "CDH2", "DPP4", "CCN3"),
  "F_Fascia_Disease" = c("ITGA10", "CCN3", "DPP4", "CDH13", "PRG4", "CRTAC1", "PCOLCE2", "LGR5")
)

f1_f8_order <- c("F1_Superficial_Healthy", "F1_Superficial_Disease", "F2_Reticular_Healthy", "F2_Reticular_Disease", "F2_F3_Perivascular_Healthy", "F2_F3_Perivascular_Disease", "F3_FRC_like_Healthy", "F3_FRC_like_Disease", "F4_HF_DPEP_Healthy", "F4_HF_DPEP_Disease", "F4_HF_TNN_Healthy", "F4_HF_TNN_Disease", "F4_HF_DP_Healthy", "F4_HF_DP_Disease", "F5_Schwann_NGFR_Healthy", "F5_Schwann_NGFR_Disease", "F5_Schwann_RAMP1_Healthy", "F5_Schwann_RAMP1_Disease", "F6_Inflammatory_Myo_Disease", "F7_Myofibroblast", "F8_Fascia_like_Myo_Disease", "F_Fascia_Disease")

f1_f8_colors <- c("F1_Superficial_Healthy" = "#A6CEE3", "F1_Superficial_Disease" = "#1F77B4", "F2_Reticular_Healthy" = "#B2DF8A", "F2_Reticular_Disease" = "#33A02C", "F2_F3_Perivascular_Healthy" = "#FDBF6F", "F2_F3_Perivascular_Disease" = "#FF7F00", "F3_FRC_like_Healthy" = "#CAB2D6", "F3_FRC_like_Disease" = "#6A3D9A", "F4_HF_DPEP_Healthy" = "#FFFF99", "F4_HF_DPEP_Disease" = "#B15928", "F4_HF_TNN_Healthy"  = "#FFEDA0", "F4_HF_TNN_Disease"  = "#D95F0E", "F4_HF_DP_Healthy"   = "#E31A1C", "F4_HF_DP_Disease"   = "#bd0026", "F5_Schwann_NGFR_Healthy"  = "#FBB4AE", "F5_Schwann_NGFR_Disease"  = "#E7298A", "F5_Schwann_RAMP1_Healthy" = "#F1EEF6", "F5_Schwann_RAMP1_Disease" = "#CE1256", "F6_Inflammatory_Myo_Disease" = "#006D2C", "F7_Myofibroblast"    = "#c73333", "F8_Fascia_like_Myo_Disease"  = "#810F7C", "F_Fascia_Disease" = "#8B0000")

# --- MACROPHAGE DICTIONARY ---
macrophage_dictionary <- list(
  "M0_Non_Polarized" = c("CD68", "MRC1"),
  "M1_Inflammatory"  = c("CD80", "IL6", "CXCL9", "CXCL10"),
  "M2_Wound_Healing" = c("CD163", "MRC1", "FOLR2", "CD209", "IL10", "CCL18")
)
macrophage_order <- c("M0_Non_Polarized", "M1_Inflammatory", "M2_Wound_Healing")
macrophage_colors <- c("M0_Non_Polarized" = "#A6CEE3", "M1_Inflammatory"  = "#E31A1C", "M2_Wound_Healing" = "#33A02C")

# ==============================================================================
# 3. DYNAMIC AUTO-ANNOTATION ENGINE
# ==============================================================================
message("Running Single-Cell Level SingleR Auto-Annotation...")

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

gene_universe <- rownames(seu_obj)

# Setup variables based on celltype
if (opt$celltype == "fibroblast") {
  active_dict <- biological_dictionary
  active_colors <- f1_f8_colors
  active_order <- f1_f8_order
} else if (opt$celltype == "macrophage") {
  active_dict <- macrophage_dictionary
  active_colors <- macrophage_colors
  active_order <- macrophage_order
} else {
  message(paste("No custom dictionary for", opt$celltype, "- exiting detail annotation."))
  quit(save = "no", status = 0)
}

ref_sce <- build_custom_singler_reference(active_dict, gene_universe)
test_sce <- as.SingleCellExperiment(seu_obj, assay = "RNA")

singler_res <- SingleR(
  test = test_sce, ref = ref_sce, labels = ref_sce$label,
  assay.type.test = "logcounts", assay.type.ref = "logcounts", prune = FALSE 
)

seu_obj$SingleR_label <- singler_res$labels
active_levels <- intersect(active_order, unique(seu_obj$SingleR_label))
seu_obj$SingleR_label <- factor(seu_obj$SingleR_label, levels = active_levels)
Idents(seu_obj) <- "SingleR_label"

# --- GLOBAL UMAP ---
p_annotated_umap <- DimPlot(seu_obj, label = TRUE, repel = TRUE, cols = active_colors) + 
  ggtitle(paste("Single-Cell Annotated", toupper(opt$celltype)))
save_dual_format(p_annotated_umap, file.path(dirs$umaps_global, "UMAP_Annotated_Global"), w = 9, h = 7 )

# ==============================================================================
# 4. UNIVERSAL VALIDATION PLOTS (VlnPlots & DotPlots)
# ==============================================================================
message("Generating Validation Plots...")

for (category_name in names(active_dict)) {
  genes_to_plot <- active_dict[[category_name]]
  valid_genes <- intersect(genes_to_plot, rownames(seu_obj))
  if (length(valid_genes) == 0) next
  
  exp_data <- GetAssayData(seu_obj, assay = "RNA", layer = "data")
  valid_genes <- valid_genes[rowSums(exp_data[valid_genes, , drop = FALSE]) > 0]
  if (length(valid_genes) == 0) next
  
  try({
    p_dot <- DotPlot(seu_obj, features = valid_genes) + RotatedAxis() + ggtitle(paste(toupper(opt$celltype), "-", category_name))
    ggsave(file.path(dirs$val_dot, paste0("DotPlot_", category_name, ".pdf")), plot = p_dot, width = max(6, length(valid_genes)*0.5 + 2), height = 6)
  }, silent = TRUE)
  
  try({
    p_vln <- VlnPlot(seu_obj, features = valid_genes, stack = TRUE, flip = TRUE, cols = active_colors) + theme(legend.position = "none") + ggtitle(paste(toupper(opt$celltype), "-", category_name))
    ggsave(file.path(dirs$val_vln, paste0("VlnPlot_", category_name, ".pdf")), plot = p_vln, width = 8, height = max(6, length(valid_genes)*1.5))
  }, silent = TRUE)
}

# ==============================================================================
# 5. FIBROBLAST-SPECIFIC DOWNSTREAM TASKS (Bypassed for Macrophages)
# ==============================================================================
if (opt$celltype == "fibroblast") {
  message("Running Fibroblast-specific downstream analysis (Mucin, Highlights, Lineages)...")
  
  # Publication Master DotPlot
  paper_signature_genes <- c("APCDD1", "COL18A1", "CRABP1", "CYP26B1", "WNT5A", "PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PPARG", "CXCL12", "CCL19", "CD74", "CXCL9", "DPEP1", "TNN", "CORIN", "RAMP1", "NGFR", "IL11", "CXCL8", "MMP1", "ACTA2", "COL11A1", "COMP", "LRRC15", "ACAN", "ITGA10", "THBS4")
  valid_paper_genes <- intersect(paper_signature_genes, rownames(seu_obj))
  master_dotplot <- DotPlot(seu_obj, features = valid_paper_genes, dot.scale = 6) + theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) + labs(title = "Fibroblast Subpopulation Signatures in PMH", x = "Key Marker Genes", y = "Identified Subclusters") + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"), axis.text.y = element_text(size = 11, face = "bold"), legend.position = "right")
  save_dual_format(master_dotplot, file.path(dirs$master_sum, "Publication_Master_DotPlot"), w = 14, h = 8)

  # Mucin & ECM Gene Analysis
  mucin_ecm_genes <- c("MUC1", "HAS1", "HAS2", "MMP1", "MUC12", "HAS3", "VCAN", "FN1", "CEMIP", "HYAL1", "HYAL2", "CTGF", "TGFBI", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "SPARC", "POSTN", "ACTA2", "TAGLN", "LOX", "LOXL2")
  available_mucin <- intersect(mucin_ecm_genes, rownames(seu_obj))
  if(length(available_mucin) > 0) {
    mucin_dotplot <- DotPlot(seu_obj, features = available_mucin, dot.scale = 8) + theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) + labs(title = "Mucin & ECM Production by Fibroblast Subtype", x = "Target Genes", y = "Fibroblast Subcluster") + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(face = "italic", color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
    save_dual_format(mucin_dotplot, file.path(dirs$mucin, "Mucin_DotPlot_Summary"), w = 10, h = 7)
  }

  # Macro-Lineages
  seu_obj$Macro_Lineage <- dplyr::case_when(
    grepl("^F1", seu_obj$SingleR_label) ~ "F1_Superficial", grepl("^F2_Reticular", seu_obj$SingleR_label) ~ "F2_Reticular", grepl("^F2_F3", seu_obj$SingleR_label) ~ "F2/F3_Perivascular", grepl("^F3", seu_obj$SingleR_label) ~ "F3_FRC_like", grepl("^F4", seu_obj$SingleR_label) ~ "F4_HairFollicle", grepl("^F5", seu_obj$SingleR_label) ~ "F5_Schwann", grepl("^F6", seu_obj$SingleR_label) ~ "F6_Inflammatory_Myo", grepl("^F7", seu_obj$SingleR_label) ~ "F7_Myofibroblast", grepl("^F8", seu_obj$SingleR_label) ~ "F8_Fascia_like_Myo", grepl("^F_Fascia", seu_obj$SingleR_label) ~ "F_Fascia", TRUE ~ "Unknown"
  )
  macro_colors <- c("F1_Superficial" = "#1F77B4", "F2_Reticular" = "#2CA02C", "F2/F3_Perivascular" = "#FF7F00", "F3_FRC_like" = "#9467BD", "F4_HairFollicle" = "#E31A1C", "F5_Schwann" = "#E7298A", "F6_Inflammatory_Myo" = "#17BECF", "F7_Myofibroblast" = "#7F7F7F", "F8_Fascia_like_Myo" = "#8C564B", "F_Fascia" = "#8B0000")
  
  p_macro_umap <- DimPlot(seu_obj, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) + ggtitle("Global UMAP: Major Fibroblast Lineages") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  save_dual_format(p_macro_umap, file.path(dirs$umaps_global, "UMAP_MacroLineage_Global"), w = 9, h = 7)
}

# ==============================================================================
# 6. SAVE DETAILED ANNOTATED RDS
# ==============================================================================
message("Saving detailed annotated RDS...")

final_rds_path <- file.path(opt$outdir, opt$celltype, "processed", paste0(opt$celltype, "_detailed_annotated.rds"))

# This prevents the "cannot open the connection" crash!
dir.create(dirname(final_rds_path), recursive = TRUE, showWarnings = FALSE)

saveRDS(seu_obj, final_rds_path)
message("=== Detail Annotation Complete! ===")