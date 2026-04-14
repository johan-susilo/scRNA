suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(slingshot)
  library(SingleCellExperiment)
  library(viridis)
})

message("\n==================================================================")
message("=== Starting Automated Intra-Lineage Trajectory Analysis ===")
message("==================================================================")

base_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"
cell_types <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

for (cell_type in cell_types) {
  
  message(paste("\n-> Running Trajectory Analysis for:", toupper(cell_type)))
  
  # 1. Setup paths
  rds_path <- file.path(base_dir, cell_type, "processed", paste0(cell_type, "_subset_processed.rds"))
  
  if (!file.exists(rds_path)) {
    message(paste("   - SKIPPING: Could not find RDS file at", rds_path))
    next
  }
  
  out_dir <- file.path(base_dir, cell_type, "trajectory")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  safe_save_plot <- function(plot_obj, base_filename, w = 8, h = 6) {
    png_path <- file.path(out_dir, paste0(base_filename, ".png"))
    tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
  }

  # 2. Load Object
  sub_obj <- readRDS(rds_path)
  
  # Ensure Detailed_Condition exists (from your cleaner script)
  if (!"Detailed_Condition" %in% colnames(sub_obj@meta.data)) {
    message("   - SKIPPING: 'Detailed_Condition' missing from metadata. Ensure the cleaner script was run.")
    next
  }
  
  # 3. Auto-Detect the "Root" Cluster (Baseline / Time 0)
  # We assume the cluster with the highest percentage of "Healthy" cells is the origin
  message("   - Auto-detecting baseline root cluster...")
  prop_table <- prop.table(table(sub_obj$seurat_clusters, sub_obj$Detailed_Condition), margin = 1)
  
  if ("Healthy" %in% colnames(prop_table)) {
    root_cluster <- names(which.max(prop_table[, "Healthy"]))
    message(paste("   - Root cluster selected:", root_cluster, "(Highest proportion of Healthy cells)"))
  } else {
    # Fallback if no healthy cells exist in this subset
    root_cluster <- levels(Idents(sub_obj))[1]
    message(paste("   - Warning: No 'Healthy' condition found. Defaulting root to cluster:", root_cluster))
  }

  # 4. Run Slingshot
  message("   - Converting to SingleCellExperiment and running Slingshot...")
  sce <- as.SingleCellExperiment(sub_obj)
  
  # Run Slingshot using the UMAP embeddings
  sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP', start.clus = root_cluster)
  
  # Extract the primary pseudotime curve (Curve 1)
  sub_obj$Pseudotime <- sce$slingPseudotime_1

  # 5. Visualizations
  message("   - Generating trajectory plots...")
  
  # Plot A: Pseudotime painted on the UMAP
  p_traj <- FeaturePlot(sub_obj, features = "Pseudotime", pt.size = 0.6) +
    scale_color_viridis(option = "C", na.value = "grey90") +
    ggtitle(paste(toupper(cell_type), "- Pseudotime Trajectory")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_traj, paste0(cell_type, "_UMAP_Pseudotime"))

  # Plot B: Disease Progression Density Plot (Satisfies Meeting Protocol)
  meta_df <- sub_obj@meta.data %>% filter(!is.na(Pseudotime))
  meta_df$Detailed_Condition <- factor(meta_df$Detailed_Condition, levels = c("Healthy", "Acute", "Chronic"))
  
  p_density <- ggplot(meta_df, aes(x = Pseudotime, fill = Detailed_Condition, color = Detailed_Condition)) +
    geom_density(alpha = 0.4, size = 1) +
    theme_minimal(base_size = 14) +
    scale_fill_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
    scale_color_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
    labs(
      title = paste(toupper(cell_type), "- Disease Progression Map"),
      subtitle = "Cellular transition across clinical stages",
      x = "Pseudotime (0 = Resting State, High = Active Disease State)",
      y = "Density of Cells"
    )
  safe_save_plot(p_density, paste0(cell_type, "_Density_Pseudotime_by_Stage"))

  # 6. Save Updated Object
  message("   - Saving RDS with Pseudotime metadata...")
  saveRDS(sub_obj, rds_path) # Overwrites the processed RDS to include the new Pseudotime column
}

message("\n==================================================================")
message("=== All Trajectory Analyses Completed! ===")
message("==================================================================")