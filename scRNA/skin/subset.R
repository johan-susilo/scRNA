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
# 1. CONFIGURATION
# ==============================================================================
BASE_OUT_DIR <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"
INPUT_RDS    <- "/home/johan/output/skin_pmh_harmony_sctransform2/TN.combined_dim30.rds"
PROOF_GENES  <- c("LUM", "DCN", "COL1A1", "KRT14", "KRT1", "PTPRC", "CD68", "PECAM1", "VWF")

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# Safely save both PDF and high-resolution PNG
save_plot <- function(plot_obj, out_dir, base_filename, w = 10, h = 6) {
  pdf_path <- file.path(out_dir, paste0(base_filename, ".pdf"))
  png_path <- file.path(out_dir, paste0(base_filename, ".png"))
  
  tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) {})
  tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
}

# Generate clean, uniform proportion barplots
create_proportion_barplot <- function(seurat_obj, group_col, title) {
  prop_df <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj[[group_col]][[1]]), margin = 2))
  colnames(prop_df) <- c("Cluster", "Group", "Proportion")
  
  ggplot(prop_df %>% filter(Proportion > 0), aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5)) +
    scale_size_continuous(range = c(2, 6), guide = "none") + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(labels = scales::percent) +
    ggtitle(title)
}

# ==============================================================================
# 3. CORE PROCESSING MODULES
# ==============================================================================

# Handles Normalization, PCA, Harmony, UMAP, and Clustering
run_dim_reduction <- function(sub_obj, resolutions) {
  
  DefaultAssay(sub_obj) <- "RNA"
  
  sub_obj <- SCTransform(
    sub_obj, 
    method = "glmGamPoi", 
    vst.flavor = "v2", 
    verbose = FALSE
  )
  
  sub_obj <- RunPCA(sub_obj, assay = "SCT", verbose = FALSE)
  
  sub_obj <- RunHarmony(
    sub_obj, 
    group.by.vars = "orig.ident2", 
    assay.use = "SCT", 
    verbose = FALSE
  )
  
  pcs_to_use <- min(20, ncol(sub_obj) - 1)
  
  # downstream steps use the Harmony reduction
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolutions, verbose = FALSE)
  
  return(sub_obj)
}

# Handles Metadata Assignment
assign_metadata <- function(sub_obj) {
  sub_obj$Detailed_Condition <- case_when(
    grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Healthy",
    grepl("AC", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Acute",
    grepl("CH", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Chronic",
    TRUE ~ "Unknown"
  )
  sub_obj$Condition <- ifelse(sub_obj$Detailed_Condition == "Healthy", "Healthy", "PMH")
  return(sub_obj)
}

# ==============================================================================
# 4. MAIN PIPELINE WRAPPER
# ==============================================================================
process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir, 
                           resolutions = seq(0.2, 1.2, by = 0.2), default_res = 0.2) {
  
  message(sprintf("\n=== Starting Pipeline for: %s ===", toupper(prefix)))
  
  # setup directories
  out_dirs <- list(
    processed = file.path(out_base_dir, prefix, "processed"),
    plots     = file.path(out_base_dir, prefix, "plots"),
    tables    = file.path(out_base_dir, prefix, "tables"),
    dge       = file.path(out_base_dir, prefix, "dge")
  )
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # initial subset
  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  
  # dimensional reduction
  sub_obj <- run_dim_reduction(sub_obj, resolutions)
  
  # set initial identity
  cluster_prefix <- ifelse(paste0("SCT_snn_res.", default_res) %in% colnames(sub_obj@meta.data), 
                           "SCT_snn_res.", 
                           "RNA_snn_res.")
  
  res_col <- paste0(cluster_prefix, default_res)
  Idents(sub_obj) <- res_col
  sub_obj$seurat_clusters <- sub_obj[[res_col]]
  
  # clustree
  p_tree <- clustree(sub_obj, prefix = cluster_prefix) +
    ggtitle(paste(prefix, "- Clustree Resolution Tracker")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p_tree, out_dirs$plots, paste0(prefix, "_clustree"), w = 12, h = 10)

  # validation dotplot
  valid_features <- intersect(PROOF_GENES, rownames(sub_obj))
  p_val <- DotPlot(sub_obj, features = valid_features) + RotatedAxis() + 
    ggtitle(paste(prefix, "Lineage Validation (Pre-Cleaning)"))
  save_plot(p_val, out_dirs$plots, paste0(prefix, "_validation_dotplot"), w = 10, h = 5)

  # metadata assignment
  sub_obj <- assign_metadata(sub_obj)

  # visuals
  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 5) +
    ggtitle(paste(prefix, "- Split by Condition (Res:", default_res, ")"))
  save_plot(p_cond, out_dirs$plots, paste0(prefix, "_umap_condition"), w = 12, h = 6)
  
  p_id2 <- DimPlot(sub_obj, split.by = "orig.ident2", label = TRUE, ncol = 3)
  save_plot(p_id2, out_dirs$plots, paste0(prefix, "_umap_origident2"), w = 15, h = 10)

  p_bar_stage <- create_proportion_barplot(sub_obj, "Detailed_Condition", paste(prefix, "Composition by Disease Stage"))
  save_plot(p_bar_stage, out_dirs$plots, paste0(prefix, "_stage_proportions"), w = 8, h = 7)

  p_bar_sample <- create_proportion_barplot(sub_obj, "orig.ident2", paste(tools::toTitleCase(prefix), "Composition by Sample"))
  save_plot(p_bar_sample, out_dirs$plots, paste0(prefix, "_sample_proportions"), w = 10, h = 7)

  # dge using RNA assay
  DefaultAssay(sub_obj) <- "RNA"

  try({ sub_obj <- JoinLayers(sub_obj) }, silent = TRUE)
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)

  all_markers <- FindAllMarkers(
    sub_obj, 
    assay = "RNA", 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25, 
    verbose = FALSE
  )

  if (nrow(all_markers) > 0) {
    write.csv(all_markers, file.path(out_dirs$dge, paste0(prefix, "_markers.csv")), row.names = FALSE)
    
    # isolate the top 10 genes per cluster
    top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    
    # scale just the top 10 genes for the heatmap to save memory/time
    sub_obj <- ScaleData(sub_obj, features = top10$gene, verbose = FALSE)
    
    p_heat <- DoHeatmap(sub_obj, features = top10$gene, assay = "RNA") + NoLegend() + 
      ggtitle(paste(prefix, "Top 10 Markers (Res:", default_res, ")"))
    
    save_plot(p_heat, out_dirs$plots, paste0(prefix, "_heatmap_top10"), w = 15, h = 15)
  }

  # tables 
  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")))
  saveRDS(sub_obj, file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds")))
  
  message(sprintf("=== Finished Pipeline for: %s ===\n", toupper(prefix)))
  return(sub_obj)
}

# ==============================================================================
# 5. EXECUTION BLOCK
# ==============================================================================

# Load Master Object
pmh_obj <- readRDS(INPUT_RDS)
pmh_obj <- subset(pmh_obj, subset = orig.ident2 != "HTY244")
Idents(pmh_obj) <- "seurat_clusters"

# Run Fibroblasts
fibro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "3", "8"), 
  prefix = "fibroblast", 
  out_base_dir = BASE_OUT_DIR,
  default_res = 1
)

#macro_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("6"), 
#  prefix = "macrophage", 
#  out_base_dir = BASE_OUT_DIR,
#  default_res = 0.2
#)

#mast_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("13"), 
#  prefix = "mast_cell", 
#  out_base_dir = BASE_OUT_DIR,
#  default_res = 0.2
#)

#
#
#fibro_macro_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("1", "3", "8","6", "14"), 
#  prefix = "fibro_macrophage", 
#  out_base_dir = BASE_OUT_DIR,
#  default_res = 0.2,
#)
#
#krt_obj <- process_subset(
#  seurat_obj = pmh_obj, 
#  subset_clusters = c("0", "2", "7", "10","16"), 
#  prefix = "keratinocyte", 
#  out_base_dir = BASE_OUT_DIR,
#  default_res = 0.2
#)

