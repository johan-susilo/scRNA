#!/usr/bin/env Rscript
# Pseudotime Analysis using Monocle3
# Usage: Rscript ./pseudotime.R -r TN.combined_dim30.rds -c 0,2,3,6,15 -o output/pseudotime
# Example: Rscript ./pseudotime.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -c 0,1,2,3 -o /home/johan/output/skin_pmh/pseudotime
# Example: Rscript ./pseudotime.R -r TN.combined_dim30.rds --all_clusters -o output/pseudotime

Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to RDS file (TN.combined_dim30.rds)"),
  make_option(c("-c", "--clusters"), type = "character", default = NULL,
              help = "Comma-separated list of cluster IDs to include (e.g., '0,2,3,6,15')"),
  make_option(c("--all_clusters"), action = "store_true", default = FALSE,
              help = "Use all clusters for pseudotime analysis"),
  make_option(c("-o", "--output"), type = "character", default = "output/pseudotime",
              help = "Output directory [default: output/pseudotime]"),
  make_option(c("--method"), type = "character", default = "monocle3",
              help = "Trajectory inference method: monocle3 [default: monocle3]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Create output directory
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(output_dir, paste0("pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")

message("============================================================")
message("Pseudotime Analysis Pipeline")
message("============================================================")
message("Log file: ", log_file)
message("Started at: ", Sys.time())

# Helper function to safely save PDF plots
safe_save_pdf <- function(plot_obj, filepath, w = 15, h = 15) {
  tryCatch({
    pdf(filepath, width = w, height = h)
    on.exit(dev.off(), add = TRUE)
    print(plot_obj)
    message("Saved plot: ", filepath)
  }, error = function(e) {
    message("Warning: Failed to save plot ", filepath, ": ", conditionMessage(e))
    if (length(dev.list()) > 0) dev.off()
  })
}

# Load and subset data ------------------------------------------------------
message("\n============================================================")
message("Loading Seurat object")
message("============================================================")

if (is.null(opt$rds)) {
  stop("RDS file path must be specified with --rds")
}

TN.combined <- readRDS(file = opt$rds)
message("Loaded Seurat object: ", ncol(TN.combined), " cells, ", nrow(TN.combined), " features")
message("Total clusters: ", length(unique(Idents(TN.combined))))

# Determine which clusters to use
if (opt$all_clusters) {
  message("Using ALL clusters for pseudotime analysis")
  TN.subset <- TN.combined
  selected_clusters <- sort(unique(Idents(TN.combined)))
} else if (!is.null(opt$clusters)) {
  # Parse cluster IDs from comma-separated string
  cluster_ids <- as.numeric(strsplit(opt$clusters, ",")[[1]])
  message("Subsetting to clusters: ", paste(cluster_ids, collapse = ", "))

  TN.subset <- subset(TN.combined, idents = cluster_ids)
  selected_clusters <- cluster_ids

  if (ncol(TN.subset) == 0) {
    stop("No cells found in the specified clusters. Check cluster IDs.")
  }

  message("Subset contains ", ncol(TN.subset), " cells")
} else {
  stop("Must specify either --clusters or --all_clusters")
}

# Plot the subset
message("\nGenerating UMAP plot of selected clusters...")
p_umap <- DimPlot(TN.subset, reduction = "umap", label = TRUE, pt.size = 0.8) +
  ggtitle(paste0("Selected Clusters for Pseudotime: ", paste(selected_clusters, collapse = ", ")))
safe_save_pdf(p_umap, file.path(output_dir, "subset_umap.pdf"))

# Rename clusters for clarity
message("\nRenaming cluster identities...")
cluster_names <- paste0("C", selected_clusters)
names(cluster_names) <- levels(TN.subset)
TN.subset <- RenameIdents(TN.subset, cluster_names)

# Update metadata cluster labels
if ("seurat_clusters" %in% colnames(TN.subset@meta.data)) {
  original_levels <- levels(TN.subset@meta.data$seurat_clusters)
  new_levels <- paste0("C", selected_clusters)
  # Create mapping for all possible levels
  full_mapping <- setNames(paste0("C", 0:(length(original_levels)-1)), original_levels)
  levels(TN.subset@meta.data$seurat_clusters) <- full_mapping[original_levels]
}

message("Cluster renaming complete")

# Convert to Monocle3 format ------------------------------------------------
message("\n============================================================")
message("Converting to Monocle3 CDS object")
message("============================================================")

# Set default assay to RNA
DefaultAssay(TN.subset) <- "RNA"

# Join layers to get proper data matrix
message("Joining RNA layers...")
Joined_subset <- JoinLayers(TN.subset)

# Extract expression data
message("Extracting expression data...")
data <- as(as.matrix(Joined_subset@assays$RNA$data), 'sparseMatrix')

# Create AnnotatedDataFrame for phenotype data (cell metadata)
message("Creating phenotype data...")
pd <- new('AnnotatedDataFrame', data = Joined_subset@meta.data)

# Create AnnotatedDataFrame for feature data (gene metadata)
message("Creating feature data...")
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# Create CellDataSet object for Monocle3
message("Creating Monocle3 CellDataSet...")
cds <- new_cell_data_set(data,
                         cell_metadata = Joined_subset@meta.data,
                         gene_metadata = fData)

message("CDS object created successfully")
message("Cells: ", ncol(cds))
message("Genes: ", nrow(cds))

# Preprocessing and dimensionality reduction --------------------------------
message("\n============================================================")
message("Preprocessing and dimensionality reduction")
message("============================================================")

message("Preprocessing CDS...")
cds <- preprocess_cds(cds, num_dim = 50)

message("Reducing dimensions with UMAP...")
cds <- reduce_dimension(cds)

message("Plotting CDS UMAP colored by cluster...")
p_cds_cluster <- plot_cells(cds, color_cells_by = "seurat_clusters",
                            label_cell_groups = TRUE,
                            label_leaves = FALSE,
                            label_branch_points = FALSE,
                            graph_label_size = 3) +
  ggtitle("Monocle3 UMAP - Colored by Cluster")

safe_save_pdf(p_cds_cluster, file.path(output_dir, "monocle3_umap_clusters.pdf"))

# Cluster cells and learn graph --------------------------------------------
message("\n============================================================")
message("Learning trajectory graph")
message("============================================================")

message("Clustering cells...")
cds <- cluster_cells(cds)

message("Learning trajectory graph...")
cds <- learn_graph(cds)

message("Plotting trajectory...")
p_trajectory <- plot_cells(cds,
                          color_cells_by = "seurat_clusters",
                          label_groups_by_cluster = FALSE,
                          label_leaves = TRUE,
                          label_branch_points = TRUE,
                          graph_label_size = 3) +
  ggtitle("Learned Trajectory")

safe_save_pdf(p_trajectory, file.path(output_dir, "trajectory_by_cluster.pdf"))

# Order cells in pseudotime -------------------------------------------------
message("\n============================================================")
message("Ordering cells in pseudotime")
message("============================================================")

message("NOTE: Automatic root cell selection...")
message("Monocle3 will select root automatically based on the earliest cluster")

# Order cells - Monocle3 will pick root automatically
cds <- order_cells(cds)

message("Pseudotime ordering complete")

# Plot pseudotime
message("\nPlotting pseudotime...")
p_pseudotime <- plot_cells(cds,
                           color_cells_by = "pseudotime",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           label_roots = TRUE,
                           graph_label_size = 3) +
  ggtitle("Cells Ordered by Pseudotime")

safe_save_pdf(p_pseudotime, file.path(output_dir, "pseudotime.pdf"))

# Combined trajectory plots
message("\nGenerating combined trajectory plots...")
p_combined <- plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster = TRUE,
                         label_leaves = TRUE,
                         label_branch_points = TRUE,
                         label_roots = TRUE,
                         graph_label_size = 3) +
  ggtitle("Trajectory with Cluster Labels")

safe_save_pdf(p_combined, file.path(output_dir, "trajectory_combined.pdf"))

# Split by sample if orig.ident1 exists
if ("orig.ident1" %in% colnames(colData(cds))) {
  message("\nPlotting trajectory split by sample group...")
  p_by_sample <- plot_cells(cds,
                            color_cells_by = "orig.ident1",
                            label_cell_groups = FALSE,
                            label_leaves = FALSE,
                            label_branch_points = TRUE,
                            graph_label_size = 3) +
    ggtitle("Trajectory Colored by Sample")

  safe_save_pdf(p_by_sample, file.path(output_dir, "trajectory_by_sample.pdf"))
}

# Save results --------------------------------------------------------------
message("\n============================================================")
message("Saving results")
message("============================================================")

# Save CDS object
message("Saving Monocle3 CDS object...")
saveRDS(cds, file.path(output_dir, "monocle3_cds.rds"))

# Extract pseudotime values and save
message("Extracting pseudotime values...")
pseudotime_df <- data.frame(
  cell = colnames(cds),
  pseudotime = pseudotime(cds),
  cluster = colData(cds)$seurat_clusters
)

# Add sample info if available
if ("orig.ident1" %in% colnames(colData(cds))) {
  pseudotime_df$sample <- colData(cds)$orig.ident1
}

write.csv(pseudotime_df, file.path(output_dir, "pseudotime_values.csv"), row.names = FALSE)
message("Pseudotime values saved")

# Summary statistics
message("\n============================================================")
message("Pseudotime Summary Statistics")
message("============================================================")
message("Min pseudotime: ", round(min(pseudotime(cds), na.rm = TRUE), 3))
message("Max pseudotime: ", round(max(pseudotime(cds), na.rm = TRUE), 3))
message("Mean pseudotime: ", round(mean(pseudotime(cds), na.rm = TRUE), 3))
message("Median pseudotime: ", round(median(pseudotime(cds), na.rm = TRUE), 3))

# Pseudotime by cluster
pseudotime_by_cluster <- pseudotime_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    min_pseudotime = min(pseudotime, na.rm = TRUE),
    max_pseudotime = max(pseudotime, na.rm = TRUE)
  )

write.csv(pseudotime_by_cluster, file.path(output_dir, "pseudotime_by_cluster.csv"), row.names = FALSE)
message("\nPseudotime statistics by cluster:")
print(pseudotime_by_cluster)

# Create violin plot of pseudotime by cluster
message("\nGenerating pseudotime violin plot...")
p_violin <- ggplot(pseudotime_df, aes(x = cluster, y = pseudotime, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pseudotime Distribution by Cluster",
       x = "Cluster",
       y = "Pseudotime") +
  coord_flip()

safe_save_pdf(p_violin, file.path(output_dir, "pseudotime_violin_by_cluster.pdf"))

# Find genes that change as a function of pseudotime ------------------------
message("\n============================================================")
message("Finding pseudotime-dependent genes")
message("============================================================")

message("This may take a while for large datasets...")
message("Testing genes for pseudotime dependence using graph_test...")

# Use graph_test to find genes that vary along the trajectory
tryCatch({
  gene_fits <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)

  # Filter significant genes
  sig_genes <- gene_fits %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)

  message("Found ", nrow(sig_genes), " genes significantly changing along pseudotime")

  # Save results
  write.csv(gene_fits, file.path(output_dir, "pseudotime_gene_fits.csv"), row.names = FALSE)
  write.csv(sig_genes, file.path(output_dir, "pseudotime_significant_genes.csv"), row.names = FALSE)

  # Plot top genes
  if (nrow(sig_genes) > 0) {
    message("\nPlotting top pseudotime-dependent genes...")
    top_genes <- head(sig_genes$gene_short_name, 9)

    p_genes <- plot_genes_in_pseudotime(cds[top_genes,],
                                        color_cells_by = "seurat_clusters",
                                        min_expr = 0.5,
                                        ncol = 3)

    safe_save_pdf(p_genes, file.path(output_dir, "top_pseudotime_genes.pdf"), w = 18, h = 18)
  }

}, error = function(e) {
  message("Warning: Could not compute pseudotime-dependent genes: ", e$message)
})

# Final summary -------------------------------------------------------------
message("\n============================================================")
message("Pseudotime Analysis Complete!")
message("============================================================")
message("Output directory: ", output_dir)
message("Files generated:")
message("  - subset_umap.pdf: UMAP of selected clusters")
message("  - monocle3_umap_clusters.pdf: Monocle3 UMAP colored by cluster")
message("  - trajectory_by_cluster.pdf: Learned trajectory")
message("  - pseudotime.pdf: Cells colored by pseudotime")
message("  - trajectory_combined.pdf: Combined trajectory visualization")
message("  - pseudotime_values.csv: Pseudotime values for each cell")
message("  - pseudotime_by_cluster.csv: Summary statistics by cluster")
message("  - pseudotime_violin_by_cluster.pdf: Distribution of pseudotime")
message("  - monocle3_cds.rds: Saved Monocle3 CDS object")

if (file.exists(file.path(output_dir, "pseudotime_significant_genes.csv"))) {
  message("  - pseudotime_gene_fits.csv: All gene test results")
  message("  - pseudotime_significant_genes.csv: Significant genes")
  message("  - top_pseudotime_genes.pdf: Expression of top genes")
}

message("\n============================================================")
message("Completed at: ", Sys.time())
message("============================================================")

# Close log file
sink(type = "message")
sink(type = "output")
close(log_conn)
