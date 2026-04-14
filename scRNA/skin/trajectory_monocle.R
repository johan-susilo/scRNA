suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(monocle3)
  library(viridis)
})

message("\n==================================================================")
message("=== Starting Automated Monocle 3 Trajectory Analysis ===")
message("==================================================================")

base_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"
cell_types <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

for (cell_type in cell_types) {
  
  message(paste("\n-> Running Monocle 3 Trajectory for:", toupper(cell_type)))
  
  # 1. Setup paths
  rds_path <- file.path(base_dir, cell_type, "processed", paste0(cell_type, "_subset_processed.rds"))
  
  if (!file.exists(rds_path)) {
    message(paste("   - SKIPPING: Could not find RDS file at", rds_path))
    next
  }
  
  out_dir <- file.path(base_dir, cell_type, "trajectory_monocle3")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  safe_save_plot <- function(plot_obj, base_filename, w = 8, h = 6) {
    png_path <- file.path(out_dir, paste0(base_filename, ".png"))
    tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
  }

  # 2. Load Object
  sub_obj <- readRDS(rds_path)
  
  if (!"Detailed_Condition" %in% colnames(sub_obj@meta.data)) {
    message("   - SKIPPING: 'Detailed_Condition' missing. Ensure cleaner script was run.")
    next
  }
  
  # 3. Convert Seurat to Monocle 3 CDS (Cell Data Set)
  message("   - Converting Seurat object to Monocle 3 CDS...")
  
  # Extract counts, metadata, and gene info
  expression_matrix <- GetAssayData(sub_obj, assay = "RNA", slot = "counts")
  cell_metadata <- sub_obj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
  rownames(gene_annotation) <- rownames(expression_matrix)
  
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  # CRITICAL: Transfer your Harmony/Seurat UMAP coordinates directly into Monocle
  # so it doesn't recalculate and ruin your structure
  reducedDims(cds)[["UMAP"]] <- sub_obj@reductions$umap@cell.embeddings
  
  # 4. Monocle 3 Graph Learning
  message("   - Clustering cells and learning trajectory graph...")
  # Monocle needs to build its own partitions based on your imported UMAP
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
  
  # 5. Auto-Detect the "Root" Cluster (Baseline / Time 0)
  message("   - Auto-detecting baseline root cluster...")
  prop_table <- prop.table(table(sub_obj$seurat_clusters, sub_obj$Detailed_Condition), margin = 1)
  
  if ("Healthy" %in% colnames(prop_table)) {
    root_cluster <- names(which.max(prop_table[, "Healthy"]))
    message(paste("   - Root cluster selected:", root_cluster, "(Highest proportion of Healthy cells)"))
  } else {
    root_cluster <- levels(Idents(sub_obj))[1]
    message(paste("   - Warning: No 'Healthy' condition found. Defaulting root to cluster:", root_cluster))
  }
  
  # Extract the actual cell barcodes belonging to the root cluster
  root_barcodes <- rownames(sub_obj@meta.data[sub_obj$seurat_clusters == root_cluster, ])

  # 6. Calculate Pseudotime
  message("   - Ordering cells by Pseudotime...")
  cds <- order_cells(cds, root_cells = root_barcodes)
  
  # Extract pseudotime and add it back to Seurat
  # Note: Monocle assigns 'Inf' to cells on disconnected islands.
  sub_obj$Pseudotime <- pseudotime(cds)

  # 7. Visualizations
  message("   - Generating trajectory plots...")
  
  # Plot A: The Classic Monocle Graph Plot (Draws the actual trajectory branches)
  p_monocle <- plot_cells(cds, 
                          color_cells_by = "pseudotime", 
                          label_cell_groups = FALSE, 
                          label_leaves = FALSE, 
                          label_branch_points = FALSE, 
                          graph_label_size = 1.5) +
    ggtitle(paste(toupper(cell_type), "- Monocle 3 Trajectory Graph")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_monocle, paste0(cell_type, "_Monocle_Trajectory_Graph"))

  # Plot B: Seurat FeaturePlot mapping Pseudotime
  p_traj <- FeaturePlot(sub_obj, features = "Pseudotime", pt.size = 0.6) +
    scale_color_viridis(option = "C", na.value = "grey90") +
    ggtitle(paste(toupper(cell_type), "- Pseudotime Heatmap")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_traj, paste0(cell_type, "_UMAP_Pseudotime_Heatmap"))

  # Plot C: Disease Progression Density Plot (Satisfies Meeting Protocol)
  # Filter out NA and 'Inf' (disconnected cells Monocle couldn't map)
  meta_df <- sub_obj@meta.data %>% filter(is.finite(Pseudotime))
  meta_df$Detailed_Condition <- factor(meta_df$Detailed_Condition, levels = c("Healthy", "Acute", "Chronic"))
  
  p_density <- ggplot(meta_df, aes(x = Pseudotime, fill = Detailed_Condition, color = Detailed_Condition)) +
    geom_density(alpha = 0.4, size = 1) +
    theme_minimal(base_size = 14) +
    scale_fill_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
    scale_color_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
    labs(
      title = paste(toupper(cell_type), "- Disease Progression Map (Monocle 3)"),
      subtitle = "Cellular transition across clinical stages",
      x = "Pseudotime (0 = Resting State, High = Active Disease State)",
      y = "Density of Cells"
    )
  safe_save_plot(p_density, paste0(cell_type, "_Density_Pseudotime_by_Stage"))

  # 8. Save Updated Object
  message("   - Saving RDS with Monocle 3 Pseudotime metadata...")
  saveRDS(sub_obj, rds_path) 
}

message("\n==================================================================")
message("=== All Monocle 3 Trajectories Completed! ===")
message("==================================================================")
