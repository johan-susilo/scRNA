#!/usr/bin/env Rscript
# Script to regenerate QC violin plots from processed RDS files
# This uses the fixed create_qc_violin_plot function to replace blank plots

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

# Helper function for QC violin plots (fixed version)
create_qc_violin_plot <- function(seurat_obj, features, title) {
  # Get cluster identity if available
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_id <- "seurat_clusters"
  } else if (length(levels(Idents(seurat_obj))) > 1) {
    # Add current idents as a column
    seurat_obj@meta.data$temp_ident <- Idents(seurat_obj)
    cluster_id <- "temp_ident"
  } else {
    cluster_id <- NULL
  }

  # Extract metadata for specified features
  if (!is.null(cluster_id)) {
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(c(features, cluster_id))) %>%
      mutate(cell = rownames(seurat_obj@meta.data)) %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')

    # Rename cluster column for consistency
    colnames(qc_data)[colnames(qc_data) == cluster_id] <- "Identity"
    qc_data$Identity <- as.factor(qc_data$Identity)

    # Generate colors based on number of clusters
    n_clusters <- length(unique(qc_data$Identity))
    if (n_clusters <= 12) {
      colors <- colorRampPalette(RColorBrewer::brewer.pal(min(n_clusters, 12), "Paired"))(n_clusters)
    } else {
      colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
    }

    # Create violin plot grouped by cluster
    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = Identity)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      scale_fill_manual(values = colors) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray")
      ) +
      labs(title = title, x = 'Identity', y = 'Value')
  } else {
    # Fallback: no clustering information available
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(features)) %>%
      mutate(cell = rownames(seurat_obj@meta.data),
             Identity = "All") %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')

    # Create simple violin plot
    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = metric)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray")
      ) +
      labs(title = title, x = '', y = 'Value')
  }

  return(p)
}

# Directories
processed_dir <- "/home/johan/output/skin_pmh/processed"
plots_base_dir <- "/home/johan/output/skin_pmh/plots"

# Get all processed RDS files
rds_files <- list.files(processed_dir, pattern = "_processed.rds$", full.names = TRUE)

cat("Found", length(rds_files), "processed samples\n\n")

# QC features to plot
qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Process each sample
for (rds_file in rds_files) {
  sample_name <- basename(rds_file)
  sample_name <- sub("_processed.rds$", "", sample_name)

  cat("Processing:", sample_name, "\n")

  # Load processed Seurat object
  seur_obj <- readRDS(rds_file)
  cat("  Loaded:", ncol(seur_obj), "cells\n")

  # Sample plot directory
  sample_plot_dir <- file.path(plots_base_dir, sample_name)

  # Check if all QC features exist
  if (!all(qc_feats %in% colnames(seur_obj@meta.data))) {
    cat("  WARNING: Missing QC metrics, skipping\n")
    next
  }

  # Since we don't have the unfiltered/prefiltered data, we can only regenerate
  # the post-doublet and final plots from the processed object
  # We'll regenerate all 4 using the final filtered data as a placeholder

  # Note: 01 and 02 would need the original unfiltered data which we don't have
  # So we'll just regenerate 08 and 09 which use the processed (final) data

  # 08: Post-doublet removal (same as final since doublets already removed)
  cat("  Generating 08_qc_violins_post_doublet.pdf...\n")
  pdf_file <- file.path(sample_plot_dir, paste0(sample_name, "_08_qc_violins_post_doublet.pdf"))
  tryCatch({
    pdf(pdf_file, width = 15, height = 5)
    p <- create_qc_violin_plot(seur_obj, qc_feats,
                               paste0("QC Metrics (Post-Doublet Removal) - ", sample_name))
    print(p)
    dev.off()
    cat("  Saved:", pdf_file, "\n")
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    if (length(dev.list()) > 0) dev.off()
  })

  # 09: Final filtered (current processed state)
  cat("  Generating 09_qc_violins_final.pdf...\n")
  pdf_file <- file.path(sample_plot_dir, paste0(sample_name, "_09_qc_violins_final.pdf"))
  tryCatch({
    pdf(pdf_file, width = 15, height = 5)
    p <- create_qc_violin_plot(seur_obj, qc_feats,
                               paste0("QC Metrics (Final Filtered) - ", sample_name))
    print(p)
    dev.off()
    cat("  Saved:", pdf_file, "\n")
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    if (length(dev.list()) > 0) dev.off()
  })

  cat("\n")
}

cat("\n========================================\n")
cat("QC plot regeneration completed!\n")
cat("Note: Plots 01 and 02 (unfiltered/prefiltered) cannot be regenerated\n")
cat("      from processed RDS files. To regenerate those, you need to\n")
cat("      delete the processed RDS files and re-run the full pipeline.\n")
cat("========================================\n")
