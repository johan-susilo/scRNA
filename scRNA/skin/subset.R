suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(clustree)
  library(harmony)
  library(ggrepel)
})

# ==============================================================================
# 1. Define the Enhanced Universal Subset Cleaner Function
# ==============================================================================
process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir, 
                           resolutions = seq(0.2, 1.2, by = 0.2), default_res = 0.2,
                           clusters_to_remove = NULL) {
  
  message(paste("\n========================================================"))
  message(paste("=== Starting Cleaner Pipeline for:", toupper(prefix), "==="))
  message(paste("========================================================"))
  
  # --- Setup Output Directories ---
  subset_dir <- file.path(out_base_dir, prefix)
  out_dirs <- list(
    processed = file.path(subset_dir, "processed"),
    plots = file.path(subset_dir, "plots"),
    tables = file.path(subset_dir, "tables"),
    dge = file.path(subset_dir, "dge")
  )
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # Helper: Save both PDF and High-Res PNG (Matched to subset.R)
  safe_save_plot <- function(plot_obj, base_filename, w = 10, h = 6) {
    pdf_path <- file.path(out_dirs$plots, paste0(base_filename, ".pdf"))
    tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) {})
    
    png_path <- file.path(out_dirs$plots, paste0(base_filename, ".png"))
    tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
  }

  # --- 1. Initial Subset ---
  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  
  # --- 2. Processing & Harmony ---
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData() %>% 
    RunPCA(verbose = FALSE)
  
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "RNA", verbose = FALSE)
  
  pcs_to_use <- min(20, ncol(sub_obj) - 1)
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolutions, verbose = FALSE)
  
  # --- 3. Clustree Visualization (From subset.R) ---
  p_tree <- clustree(sub_obj, prefix = "RNA_snn_res.") +
    ggtitle(paste(prefix, "- Clustree Resolution Tracker")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_tree, paste0(prefix, "_clustree"), w = 12, h = 10)

  # Set Initial Identity
  res_col <- paste0("RNA_snn_res.", default_res)
  Idents(sub_obj) <- res_col
  sub_obj$seurat_clusters <- sub_obj[[res_col]]

  # --- 4. Validation DotPlot (The Imposter Check) ---
  proof_genes <- c("LUM", "DCN", "COL1A1", "KRT14", "KRT1", "PTPRC", "CD68", "PECAM1", "VWF")
  valid_features <- intersect(proof_genes, rownames(sub_obj))
  p_val <- DotPlot(sub_obj, features = valid_features) + RotatedAxis() + 
    ggtitle(paste(prefix, "Lineage Validation (Pre-Cleaning)"))
  safe_save_plot(p_val, paste0(prefix, "_validation_dotplot"), w = 10, h = 5)

  # --- 5. The Cleaner (Pass 2) ---
  if (!is.null(clusters_to_remove)) {
    message(paste("CLEANING: Removing imposter clusters:", paste(clusters_to_remove, collapse = ", ")))
    sub_obj <- subset(sub_obj, idents = clusters_to_remove, invert = TRUE)
    # Re-run UMAP/Neighbors on the pure population
    sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
    sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  }

  # --- 6. Define Conditions ---
  sub_obj$Detailed_Condition <- case_when(
    grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Healthy",
    grepl("AC", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Acute",
    grepl("CH", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Chronic",
    TRUE ~ "Unknown"
  )
  sub_obj$Condition <- ifelse(sub_obj$Detailed_Condition == "Healthy", "Healthy", "PMH")

  # --- 7. Enhanced UMAP Visualizations (Matched to subset.R) ---
  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 5) +
    ggtitle(paste(prefix, "- Split by Condition (Res:", default_res, ")"))
  safe_save_plot(p_cond, paste0(prefix, "_umap_condition"), w = 12, h = 6)
  
  p_id2 <- DimPlot(sub_obj, split.by = "orig.ident2", label = TRUE, ncol = 3)
  safe_save_plot(p_id2, paste0(prefix, "_umap_origident2"), w = 15, h = 10)

  # --- 8. Enhanced Proportional Barplots (Matched to subset.R) ---
  # Condition Plot with labels
  prop_cond <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$Detailed_Condition), margin = 2))
  colnames(prop_cond) <- c("Cluster", "Condition", "Proportion")
  
  p_bar <- ggplot(prop_cond %>% filter(Proportion > 0), aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5)) +
    scale_size_continuous(range = c(2, 6), guide = "none") + 
    theme_minimal(base_size = 14) +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste(prefix, "Composition by Disease Stage"))
  safe_save_plot(p_bar, paste0(prefix, "_stage_proportions"), w = 8, h = 7)

  # --- 1. Calculate Proportions by Individual Sample ---
prop_sample <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$orig.ident2), margin = 2))
colnames(prop_sample) <- c("Cluster", "Sample", "Proportion")

# --- 2. Generate the Per-Sample Bar Plot ---
p_bar_sample <- ggplot(prop_sample %>% filter(Proportion > 0), aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5)) +
  scale_size_continuous(range = c(2, 6), guide = "none") + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilts the sample names so they are readable
  scale_y_continuous(labels = scales::percent) +
  ggtitle(paste(tools::toTitleCase(prefix), "Composition by Sample"))
  safe_save_plot(p_bar_sample, paste0(prefix, "_sample_proportions"), w = 10, h = 7)

  # --- 9. DGE & Heatmap (Matched to subset.R) ---
  sub_obj <- JoinLayers(sub_obj)
  all_markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  if (nrow(all_markers) > 0) {
    write.csv(all_markers, file.path(out_dirs$dge, paste0(prefix, "_markers.csv")), row.names = FALSE)
    top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    p_heat <- DoHeatmap(sub_obj, features = top10$gene) + NoLegend() + 
      ggtitle(paste(prefix, "Top 10 Markers (Res:", default_res, ")"))
    safe_save_plot(p_heat, paste0(prefix, "_heatmap_top10"), w = 15, h = 15)
  }

  # --- 10. Tables & Save ---
  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")))
  
  saveRDS(sub_obj, file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds")))
  return(sub_obj)
}

# ==============================================================================
# 2. Execution
# ==============================================================================
base_out_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster_strict"
pmh_obj <- readRDS("/home/johan/output/skin_pmh_harmony_sctransform2/TN.combined_dim30.rds")
pmh_obj <- subset(pmh_obj, subset = orig.ident2 != "HTY244")
Idents(pmh_obj) <- "seurat_clusters"

# Run Fibroblasts
fibro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "3", "8"), 
  prefix = "fibroblast", 
  out_base_dir = base_out_dir,
  default_res = 0.6
)


## Run Pipeline 3: Macrophages
macro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("6"), 
  prefix = "macrophage", 
  out_base_dir = base_out_dir,
  default_res = 0.2
)
#
#
#fibro_macro_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("1", "3", "8","6", "14"), 
#  prefix = "fibro_macrophage", 
#  out_base_dir = base_out_dir,
#  default_res = 0.2,
#  clusters_to_remove = c("6")
#)
#
#krt_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("0", "2", "7", "10","16"), 
#  prefix = "keratinocyte", 
#  out_base_dir = base_out_dir,
#  default_res = 0.2
#)

