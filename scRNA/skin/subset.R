suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ==============================================================================
# 1. Define the Universal Subset Processing Function
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ==============================================================================
# 1. Define the Universal Subset Processing Function
# ==============================================================================
process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir, resolution = 0.4) {
  
  if (!dir.exists(out_base_dir)) {
    dir.create(out_base_dir, recursive = TRUE)
  }

  message(paste("\n========================================================"))
  message(paste("=== Starting Processing for:", prefix, "==="))
  message(paste("========================================================"))
  
  # --- Setup Output Directories ---
  subset_dir <- file.path(out_base_dir, prefix)
  out_dirs <- list(
    processed = file.path(subset_dir, "processed"),
    plots = file.path(subset_dir, "plots"),
    tables = file.path(subset_dir, "tables")
  )
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # --- NEW: Helper function to safely save BOTH PDF and PNG ---
  safe_save_plot <- function(plot_obj, base_filename, w = 10, h = 6) {
    # 1. Save PDF
    pdf_path <- file.path(out_dirs$plots, paste0(base_filename, ".pdf"))
    tryCatch({
      pdf(pdf_path, width = w, height = h)
      print(plot_obj)
      dev.off()
    }, error = function(e) {
      message("Failed to save PDF: ", pdf_path, " : ", e$message)
      if (length(dev.list()) > 0) dev.off()
    })
    
    # 2. Save High-Res PNG
    png_path <- file.path(out_dirs$plots, paste0(base_filename, ".png"))
    tryCatch({
      # units="in" and res=300 ensures publication-quality PNGs
      png(png_path, width = w, height = h, units = "in", res = 300)
      print(plot_obj)
      dev.off()
    }, error = function(e) {
      message("Failed to save PNG: ", png_path, " : ", e$message)
      if (length(dev.list()) > 0) dev.off()
    })
  }

  # --- 1. Subset the Object ---
  message(paste("Subsetting clusters:", paste(subset_clusters, collapse = ", ")))
  
  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }
  
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  
  # --- 2. Re-process the RNA Assay ---
  message("Re-processing RNA assay (Normalize, FindVar, Scale, PCA, UMAP)...")
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, dims = 1:20, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = 1:20, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolution, verbose = FALSE)
  
  Idents(sub_obj) <- "seurat_clusters"
  
  # --- 3. Define Clinical Conditions ---
  message("Assigning clinical conditions...")
  sub_obj$Condition <- ifelse(grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE), 
                              "Healthy", "PMH")
  
  # --- 4. Generate UMAP Visualizations ---
  message("Generating UMAPs...")
  
  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 5, pt.size = 0.8) +
    ggtitle(paste("Figure 1:", prefix, "- Split by Condition"))
  safe_save_plot(p_cond, paste0(prefix, "_umap_condition"), w = 12, h = 6)
  
  p_id1 <- DimPlot(sub_obj, split.by = "orig.ident1", label = TRUE, label.size = 5, pt.size = 0.8, ncol = 2) +
    ggtitle(paste("Identification of PMH-Specific", prefix, "Population"))
  safe_save_plot(p_id1, paste0(prefix, "_umap_origident1"), w = 12, h = 10)
  
  p_id2 <- DimPlot(sub_obj, split.by = "orig.ident2", label = TRUE, label.size = 5, pt.size = 0.8, ncol = 2) +
    ggtitle(paste("Identification of PMH-Specific", prefix, "Population"))
  safe_save_plot(p_id2, paste0(prefix, "_umap_origident2"), w = 12, h = 10)
  
  # --- 5. Proportional Composition Plots ---
  message("Generating stacked barplots...")
  
  # A. Proportions by Condition (Healthy vs PMH)
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
  
  # B. Proportions by Individual Sample (HTY131, UA258, AC259, CH260)
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
  safe_save_plot(p_bar_samp, paste0(prefix, "_stacked_barplot_sample"), w = 8, h = 6)
  
  # --- 6. Tracking Tables ---
  message("Saving reference tables...")
  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")), row.names = TRUE)
  
  # --- 7. Save Processed Object ---
  rds_out <- file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds"))
  saveRDS(sub_obj, rds_out)
  message(paste("=== Successfully processed and saved:", rds_out, "==="))
  
  return(sub_obj)
}

# ==============================================================================
# 2. Execution Block
# ==============================================================================

# Set paths
base_out_dir <- "/home/johan/output/skin_pmh/subset"
combined_rds_path <- file.path("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

# Load Parent Object
message("Loading main integrated Seurat object...")
pmh_obj <- readRDS(combined_rds_path)

# Ensure global clusters are set correctly before any subsetting
Idents(pmh_obj) <- "seurat_clusters"
pmh_obj$global_cluster <- Idents(pmh_obj)

# Run Pipeline 1: Fibroblasts Only
fibro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "5", "7", "10"), 
  prefix = "fibroblast", 
  out_base_dir = base_out_dir
)

# Run Pipeline 2: Macrophages Only
macro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("6"), 
  prefix = "macrophage", 
  out_base_dir = base_out_dir
)

# Run Pipeline 3: Macrophages + Fibroblasts
macro_fibro_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = c("1", "5", "6", "7", "10"), 
  prefix = "macro_fibro", 
  out_base_dir = base_out_dir
)

message("\nAll subset pipelines executed successfully!")