suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(clustree)
  library(harmony)
})

# ==============================================================================
# 1. Define the Universal Subset Processing Function
# ==============================================================================
process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir, 
                           resolutions = seq(0.2, 1.2, by = 0.2), default_res = 0.4) {
  
  message(paste("\n========================================================"))
  message(paste("=== Starting Pipeline for:", prefix, "==="))
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
  
  # Helper function to safely save BOTH PDF and High-Res PNG
  safe_save_plot <- function(plot_obj, base_filename, w = 10, h = 6) {
    pdf_path <- file.path(out_dirs$plots, paste0(base_filename, ".pdf"))
    tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) {})
    
    png_path <- file.path(out_dirs$plots, paste0(base_filename, ".png"))
    tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
  }

  # --- 1. Subset the Object ---
  message(paste("Subsetting global clusters:", paste(subset_clusters, collapse = ", ")))
  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  
  # --- 2. Re-process RNA Assay & Run HARMONY ---
  message("Re-processing RNA assay...")
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, verbose = FALSE)
  
  message("Running Harmony Integration...")
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "RNA", verbose = FALSE)
  
  # Dynamically calculate dimensions to prevent crashes on very small subsets
  pcs_to_use <- min(20, ncol(sub_obj) - 1)
  
  message("Running UMAP and FindNeighbors on Harmony reduction...")
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  
  # --- 3. MULTI-RESOLUTION CLUSTERING & CLUSTREE ---
  message("Running multiple clustering resolutions...")
  sub_obj <- FindClusters(sub_obj, resolution = resolutions, verbose = FALSE)
  
  message("Generating Clustree visualization...")
  p_tree <- clustree(sub_obj, prefix = "RNA_snn_res.") +
    ggtitle(paste(prefix, "- Clustree Resolution Tracker")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_tree, paste0(prefix, "_clustree"), w = 12, h = 10)
  
  # Set active identity to the requested default resolution
  default_res_col <- paste0("RNA_snn_res.", default_res)
  if(default_res_col %in% colnames(sub_obj@meta.data)) {
    Idents(sub_obj) <- default_res_col
    sub_obj$seurat_clusters <- sub_obj[[default_res_col]]
    message(paste("Active identity set to default resolution:", default_res))
  } else {
    warning("Default resolution not found. Using highest calculated resolution.")
  }

  # --- 4. Define Clinical Conditions ---
  message("Assigning clinical conditions...")
  sub_obj$Condition <- ifelse(grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE), "Healthy", "PMH")
  
  # --- 5. Generate UMAP Visualizations ---
  message("Generating UMAPs...")
  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 5, pt.size = 0.8) +
    ggtitle(paste("Figure 1:", prefix, "- Split by Condition (Res:", default_res, ")"))
  safe_save_plot(p_cond, paste0(prefix, "_umap_condition"), w = 12, h = 6)
  
  p_id1 <- DimPlot(sub_obj, split.by = "orig.ident1", label = TRUE, label.size = 5, pt.size = 0.8, ncol = 2)
  safe_save_plot(p_id1, paste0(prefix, "_umap_origident1"), w = 12, h = 10)
  
  p_id2 <- DimPlot(sub_obj, split.by = "orig.ident2", label = TRUE, label.size = 5, pt.size = 0.8, ncol = 2)
  safe_save_plot(p_id2, paste0(prefix, "_umap_origident2"), w = 12, h = 10)
  
  # --- 6. Proportional Composition Plots (Borderless + Resized Text) ---
  message("Generating stacked barplots...")
  
  # Condition Plot
  prop_cond <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$Condition), margin = 2))
  colnames(prop_cond) <- c("Cluster", "Condition", "Proportion")
  prop_cond_plot <- prop_cond %>% filter(Proportion > 0)
  
  p_bar_cond <- ggplot(prop_cond_plot, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5), color = "black") +
    scale_size_continuous(range = c(2, 6), guide = "none") + 
    theme_minimal(base_size = 14) +
    labs(title = paste(prefix, "- Proportions by Condition"), y = "Percentage of Total Cells") +
    scale_y_continuous(labels = scales::percent)
  safe_save_plot(p_bar_cond, paste0(prefix, "_stacked_barplot_condition"), w = 8, h = 6)
  
  # Sample Plot
  sample_vector <- as.character(sub_obj@meta.data[["orig.ident2"]])
  prop_samp <- as.data.frame(prop.table(table(Idents(sub_obj), sample_vector), margin = 2))
  colnames(prop_samp) <- c("Cluster", "Sample", "Proportion")
  prop_samp_plot <- prop_samp %>% filter(Proportion > 0)
  
  p_bar_samp <- ggplot(prop_samp_plot, aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5), color = "black") +
    scale_size_continuous(range = c(2, 6), guide = "none") + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = paste(prefix, "- Proportions per Sample"), y = "Percentage of Total Cells") +
    scale_y_continuous(labels = scales::percent)
  safe_save_plot(p_bar_samp, paste0(prefix, "_stacked_barplot_sample"), w = 10, h = 6)
  
  # --- 7. AUTOMATED DIFFERENTIAL GENE EXPRESSION (DGE) ---
  message("Joining Seurat v5 layers for DGE...")
  # THIS IS THE FIX: Merge the sample layers back together
  sub_obj <- JoinLayers(sub_obj)
  
  message("Running DGE (FindAllMarkers) for default resolution...")
  all_markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  # Safety check: Only try to group and plot if FindAllMarkers actually found genes!
  if (nrow(all_markers) > 0) {
    write.csv(all_markers, file.path(out_dirs$dge, paste0(prefix, "_res", default_res, "_FindAllMarkers.csv")), row.names = FALSE)
    
    # Extract top 10 markers per cluster for a quick heatmap overview
    top10_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    p_heat <- DoHeatmap(sub_obj, features = top10_markers$gene) + NoLegend() + 
      ggtitle(paste("Top 10 Markers per Cluster (Res:", default_res, ")"))
    safe_save_plot(p_heat, paste0(prefix, "_heatmap_top10_res", default_res), w = 15, h = 15)
  } else {
    message("Warning: No significant DE genes found at this resolution. Skipping heatmap.")
  }
  
  # --- 8. Tracking Tables & Save Object ---
  message("Saving reference tables and final RDS...")
  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")), row.names = TRUE)
  
  rds_out <- file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds"))
  saveRDS(sub_obj, rds_out)
  message(paste("=== Successfully processed and saved:", rds_out, "==="))
  
  return(sub_obj)
}

# ==============================================================================
# 2. Execution Block
# ==============================================================================

# Set paths
base_out_dir <- "/home/johan/output/skin_pmh/subset_cluster"
combined_rds_path <- "/home/johan/output/skin_pmh/TN.combined_dim30.rds"

# Load Parent Object
message("Loading main integrated Seurat object...")
pmh_obj <- readRDS(combined_rds_path)

# VERY IMPORTANT: Remove the chemotherapy-contaminated healthy sample (HTY244)
message("Removing chemotherapy sample HTY244 to ensure a clean baseline...")
pmh_obj <- subset(pmh_obj, subset = orig.ident2 != "HTY244")

# Ensure global clusters are set correctly before subsetting
Idents(pmh_obj) <- "seurat_clusters"
pmh_obj$global_cluster <- Idents(pmh_obj)

# Run Pipeline 1: Fibroblasts
fibro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "5", "7", "10"), 
  prefix = "fibroblast", 
  out_base_dir = base_out_dir,
  default_res = 0.4
)

# Run Pipeline 2: Keratinocytes
krt_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("0", "2", "8", "11"), 
  prefix = "keratinocyte", 
  out_base_dir = base_out_dir,
  default_res = 0.2
)

# Run Pipeline 3: Macrophages
macro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("6"), 
  prefix = "macrophage", 
  out_base_dir = base_out_dir,
  default_res = 0.4 
)

message("\nAll Gold Standard subset pipelines executed successfully!")