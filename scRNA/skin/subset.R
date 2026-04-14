suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(clustree)
  library(harmony)
  library(ggrepel)
  library(SingleR)
  library(SingleCellExperiment)
})

# ==============================================================================
# 1. Nature Immunology 2025 Atlas Dictionary
# ==============================================================================
nature_atlas_dictionary <- list(
  # Non-Fibroblast Outgroups
  "Macrophage" = c("CD14", "CD68", "CD163", "FCGR3A", "MRC1", "C1QA", "C1QB", "C1QC", "TREM2", "APOE", "F13A1", "LYVE1"),
  "Monocyte_DC" = c("FCN1", "S100A8", "S100A9", "VCAN", "CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "HLA-DQA1"),
  "T_Cells_General" = c("CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A"),
  "B_Plasma_Mast" = c("CD19", "MS4A1", "SDC1", "CD38", "JCHAIN", "KIT", "TPSAB1", "CPA3"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "KDR", "ACKR1", "PROX1", "LYVE1", "PDPN"),
  "Keratinocyte_Melanocyte" = c("KRT5", "KRT14", "TP63", "KRT1", "KRT10", "FLG", "MLANA", "PMEL", "MITF"),
  
  # Healthy Fibroblasts
  "F1_Superficial" = c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "COMP", "NKD2", "RSPO1", "AXIN2", "WIF1", "SFRP2", "IGFBP2", "TNFRSF21", "PDE4B", "GUCY1A1"),
  "F2_Universal" = c("CD34", "PI16", "MFAP5", "DPP4", "PCOLCE2", "CTHRC1", "SLPI", "CD70", "LGR5"),
  "F2_3_Perivascular" = c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10"),
  "F3_FRC_like" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "IRF8", "IL15", "VCAM1", "HLA-DRA", "HLA-DRB1"),
  "F4_DS_DPEP1" = c("ASPN", "COL11A1", "MEF2C", "DPEP1", "MYL4"),
  "F4_TNN_COCH" = c("TNN", "COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD"),
  "F4_DP_HHIP" = c("CORIN", "BMP7", "WNT5A", "LEF1", "HHIP", "RSPO3", "INHBA", "PTCH1"),
  "F5_RAMP1" = c("SCN7A", "FMO2", "FGFBP2", "OLFML2A", "PEAR1", "RAMP1", "RELN", "PLEKHA6", "IGFBP2", "SFRP1", "EBF2"),
  "F5_NGFR" = c("NGFR", "SFRP4", "ITGA6", "CDH19", "CLDN1"),
  
  # Disease-Adapted & Specific Fibroblasts
  "F1_Disease_Adapted" = c("CRABP1", "CYP26B1", "CXCL1", "WNT5A", "COL7A1", "GGT5"),
  "F3_Disease_Adapted" = c("CXCL9", "ADAMDEC1"),
  "F6_Inflammatory_Myofibroblast" = c("IL11", "IL24", "IL7R", "CXCL5", "CXCL6", "CXCL8", "CXCL13", "CCL11", "MMP1", "MMP3", "CSF3", "TDO2", "WWC1", "CHI3L1", "STAT4", "CCL5", "FAM167A", "HIF1A", "WNT2", "COL10A1"),
  "F7_Myofibroblast" = c("ACTA2", "TAGLN", "ASPN", "COMP", "COL11A1", "KIF26B", "ZNF469", "RUNX2", "SULF1", "ADAM12", "ADAM19", "COL8A1", "LRRC15", "CREB3L1", "CTHRC1", "CCN4", "FABP5", "CDH2", "C1QTNF3", "CADM1", "LRRC17", "PIEZO2", "SFRP4", "KCNMA1", "SCX"),
  "F8_Fascia_like_Myofibroblast" = c("ACAN", "ITGA10", "THBS4", "FGF18", "PRG4", "CRTAC1"),
  
  # Stress Markers (No LUM/DCN to prevent signal drowning)
  "QC_Stress_Markers" = c("MT2A", "MT1M", "MT1X", "HSP90AA1", "JUNB", "GADD45B", "IER3") 
)

# ==============================================================================
# 2. SingleR Reference Builder
# ==============================================================================
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

# ==============================================================================
# 3. Main Detailed Annotation & Cleaner Function
# ==============================================================================
process_detailed_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir, 
                                    resolutions = seq(0.2, 1.2, by = 0.2), default_res = 0.2,
                                    clusters_to_remove = NULL) {
  
  message(paste("\n========================================================"))
  message(paste("=== Cleaning & Detailed Annotation:", toupper(prefix), "==="))
  message(paste("========================================================"))
  
  # --- Setup Output Directories ---
  subset_dir <- file.path(out_base_dir, prefix)
  out_dirs <- list(
    processed = file.path(subset_dir, "processed"),
    plots = file.path(subset_dir, "plots"),
    tables = file.path(subset_dir, "tables")
  )
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Subset & Standard Processing
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  DefaultAssay(sub_obj) <- "RNA"
  
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = rownames(sub_obj), verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "RNA", verbose = FALSE)
  pcs_to_use <- min(20, ncol(Embeddings(sub_obj, "harmony")))
  
  # Handle "Imposter" clusters if provided
  if (!is.null(clusters_to_remove)) {
    sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
    sub_obj <- FindClusters(sub_obj, resolution = default_res, verbose = FALSE)
    sub_obj <- subset(sub_obj, idents = clusters_to_remove, invert = TRUE)
  }
  
  # Final Dimensionality Reduction
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = default_res, verbose = FALSE)
  sub_obj$seurat_clusters <- as.character(Idents(sub_obj))
  
  # Assign Conditions for downstream Volcano Plots
  sub_obj$Detailed_Condition <- case_when(
    grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Healthy",
    grepl("AC", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Acute",
    grepl("CH", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Chronic",
    TRUE ~ "Unknown"
  )
  
  # 2. Apply SingleR Detailed Annotation (Nature Atlas)
  message("Running Detailed SingleR Annotation...")
  gene_universe <- rownames(sub_obj)
  ref_sce <- build_custom_singler_reference(nature_atlas_dictionary, gene_universe)
  test_sce <- as.SingleCellExperiment(sub_obj, assay = "RNA")
  
  singler_res <- SingleR(
    test = test_sce, ref = ref_sce, labels = ref_sce$label,
    clusters = sub_obj$seurat_clusters, assay.type.test = "logcounts", prune = TRUE
  )
  

  # Map annotations back to Seurat object
  singler_df <- as.data.frame(singler_res)
  singler_df$cluster <- rownames(singler_df)
  
  # 1. Grab the un-pruned labels (The ones we actually want to use)
  label_map <- setNames(singler_df$labels, singler_df$cluster)
  sub_obj$SingleR_label <- unname(label_map[as.character(sub_obj$seurat_clusters)])
  
  # 2. Grab the pruned labels
  pruned_map <- setNames(ifelse(is.na(singler_df$pruned.labels), "Unassigned", singler_df$pruned.labels), singler_df$cluster)
  sub_obj$SingleR_pruned_label <- unname(pruned_map[as.character(sub_obj$seurat_clusters)])

  # 3. Save Output Table and Object
  write.csv(singler_df, file.path(out_dirs$tables, paste0(prefix, "_detailed_annotations.csv")), row.names = FALSE)
  saveRDS(sub_obj, file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds")))
  
  message(paste("Success! Object saved to:", file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds"))))
  return(sub_obj)
}

# ==============================================================================
# 4. Execution
# ==============================================================================
base_out_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster_correct"
pmh_obj <- readRDS("/home/johan/output/skin_pmh_harmony_sctransform2/TN.combined_dim30.rds")
pmh_obj <- subset(pmh_obj, subset = orig.ident2 != "HTY244")
Idents(pmh_obj) <- "seurat_clusters"

# Run the full cleaning & detailed annotation pipeline
fibro_obj <- process_detailed_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "3", "8"), 
  prefix = "fibroblast", 
  out_base_dir = base_out_dir,
  default_res = 0.2,
  clusters_to_remove = c("6")
)

