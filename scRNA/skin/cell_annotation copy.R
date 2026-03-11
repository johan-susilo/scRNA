#!/usr/bin/env Rscript
# Usage: Rscript ./cell_annotation.R -r TN.combined_dim30.rds -s all
# Example: Rscript ./cell_annotation.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -s scCATCH
# Example: Rscript ./cell_annotation.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -s all -o /home/johan/output/skin_pmh/annotations
# Example: Rscript ./cell_annotation.R -r TN.combined_dim30.rds -s singleR,scCATCH -c 2
# Example: Rscript ./cell_annotation.R -r TN.combined_dim30.rds -s all --consensus

Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scCATCH)
  library(SingleR)
  library(celldex)
  library(dplyr)
  library(tidyverse)
  library(CelliD)
  library(ggpubr)
  library(parallel)
  library(future)
  library(future.apply)
  library(stringr)
})

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to RDS file (TN.combined_dim30.rds)"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Pipeline step: read_rds, singleR, markers, celliD, scCATCH, consensus, annotated_plots, combined_plots, all"),
  make_option(c("-o", "--output"), type = "character", default = "annotations",
              help = "Base output directory [default: annotations]"),
  make_option(c("-p", "--plots"), type = "character", default = NULL,
              help = "Path to plots directory (for annotated_plots step)"),
  make_option(c("--consensus"), action = "store_true", default = FALSE,
              help = "Generate consensus annotations from all methods"),
  make_option(c("--tissue"), type = "character", default = "skin",
              help = "Tissue type for scCATCH [default: skin]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Treat empty string as NULL to avoid normalizePath errors
if (!is.null(opt$plots) && identical(opt$plots, "")) {
  opt$plots <- NULL
}

# Use provided output path and create directories
output_base <- opt$output
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

resolve_plots_dir <- function(base_out, user_plots) {
  # prefer user-provided plots dir if not empty/NULL
  if (!is.null(user_plots) && nzchar(user_plots)) {
    return(user_plots)
  }
  return(file.path(base_out, "plots"))
}

plots_root <- resolve_plots_dir(output_base, opt$plots)

output_dirs <- list(
  singleR = file.path(output_base, "singleR"),
  markers = file.path(output_base, "markers"),
  celliD = file.path(output_base, "celliD"),
  scCATCH = file.path(output_base, "scCATCH"),
  consensus = file.path(output_base, "consensus"),
  logs = file.path(output_base, "logs"),
  annotated_plots = file.path(plots_root, "annotated"),
  combined_plots = file.path(plots_root, "combined_plots")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(output_dirs$logs, paste0("annotation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)  # split=TRUE means output to both console and file
sink(log_conn, type = "message")

message("Log file: ", log_file)
message("Started at: ", Sys.time())

# Set up processing mode - ALWAYS use sequential to avoid deadlocks
message("Using SEQUENTIAL processing mode to avoid deadlocks")
plan(sequential)

# Helpers and stubs ----------------------------------------------------------
# Safe reader for RDS that returns a list with expected elements.
read_rds <- function(path) {
  obj <- readRDS(path)
  # If the saved object is a single Seurat, wrap it so downstream code works.
  if (inherits(obj, "Seurat")) {
    message("Wrapping single Seurat object from RDS into a list for pipeline compatibility")
    return(list(TN.combined = obj, Joined_TN.combined = obj))
  }
  # Ensure expected names exist
  if (is.null(obj$TN.combined) && !is.null(obj$Joined_TN.combined)) {
    obj$TN.combined <- obj$Joined_TN.combined
  }
  if (is.null(obj$Joined_TN.combined) && !is.null(obj$TN.combined)) {
    obj$Joined_TN.combined <- obj$TN.combined
  }
  return(obj)
}

# Define stubs if not already present (defensive guards)
if (!exists("run_singleR", mode = "function")) {
  run_singleR <- function(seu) {
    message("[stub] run_singleR: implement SingleR labeling here. No-op for now.")
    invisible(NULL)
  }
}
if (!exists("plot_markers", mode = "function")) {
  plot_markers <- function(seu) {
    message("[stub] plot_markers: implement marker plotting here. No-op for now.")
    invisible(NULL)
  }
}
if (!exists("run_celliD", mode = "function")) {
  run_celliD <- function(seu) {
    message("[stub] run_celliD: implement CelliD annotation here. No-op for now.")
    invisible(NULL)
  }
}
if (!exists("run_scCATCH", mode = "function")) {
  run_scCATCH <- function(seu_tn, seu_joined) {
    message("[stub] run_scCATCH: implement scCATCH annotation here. No-op for now.")
    invisible(NULL)
  }
}
if (!exists("generate_annotated_plots", mode = "function")) {
  generate_annotated_plots <- function(seu) {
    message("[stub] generate_annotated_plots: implement annotated plots here. No-op for now.")
    invisible(NULL)
  }
}
if (!exists("generate_combined_plots", mode = "function")) {
  generate_combined_plots <- function(seu) {
    message("[stub] generate_combined_plots: implement combined plots here. No-op for now.")
    invisible(NULL)
  }
}

# Helper function to safely save PDF plots
safe_save_pdf <- function(plot_obj, filepath, w = 15, h = 15) {
  tryCatch({
    pdf(filepath, width = w, height = h)
    on.exit(dev.off(), add = TRUE)  # Ensure device closes even if print fails
    print(plot_obj)
    message("Saved plot: ", filepath)
  }, error = function(e) {
    message("Warning: Failed to save plot ", filepath, ": ", conditionMessage(e))
    # Make sure device is closed
    if (length(dev.list()) > 0) dev.off()
  })
}

# Normalization function for cell type names (from Cell_anno_consensus.pl)
normalize_cell_type <- function(cell_type) {
  # Normalization mappings
  normalization_map <- c(
    "Monocyte" = "Monocytes",
    "Endothelial_cells" = "Endothelial cells",
    "Endothelial cells" = "Endothelial cells",
    "Macrophage" = "Macrophages",
    "DC" = "Dendritic cells",
    "Killer Cell" = "NK cells",
    "Natural Killer Cell" = "NK cells",
    "NK cell" = "NK cells",
    "T cells" = "T cells",
    "Fibroblasts" = "Fibroblasts",
    "Keratinocytes" = "Keratinocytes",
    "Epithelial cells" = "Epithelial cells"
  )

  # Replace underscores with spaces
  cell_type <- gsub("_", " ", cell_type)

  # Apply normalization mapping if exists
  if (cell_type %in% names(normalization_map)) {
    return(normalization_map[cell_type])
  }

  return(cell_type)
}

# Helper: find common words between labels (fallback when no clear consensus)
common_words_label <- function(labels) {
  # tokenize, lowercase, keep alphanumerics
  toks <- lapply(labels, function(x) {
    str_split(str_to_lower(x), "\\W+")[[1]] %>% discard(~ .x == "")
  })
  # words appearing in at least 2 labels
  tab <- table(unlist(toks))
  shared <- names(tab)[tab >= 2]
  if (length(shared) == 0) return(NULL)
  # compose phrase from shared words in original order of SingleR label where possible
  sr_words <- toks[[1]]
  ordered <- sr_words[sr_words %in% shared]
  if (length(ordered) > 0) {
    paste(ordered, collapse = " ")
  } else {
    paste(shared, collapse = " ")
  }
}

# Compute consensus label given multiple methods for one cluster
consensus_label <- function(labels, prefer = c("SingleR", "MethodB", "MethodC")) {
  # labels is a named character vector c(SingleR="...", MethodB="...", MethodC="...")
  lbls <- na.omit(labels)
  if (length(lbls) == 0) return(NA_character_)
  # frequency vote
  freq <- sort(table(lbls), decreasing = TRUE)
  if (length(freq) == 1 || freq[1] > freq[2]) {
    # normalize winning label
    return(normalize_cell_type(names(freq)[1]))
  }
  # full tie: if all labels equal length table entries and same counts
  # prefer SingleR if available
  if (!is.null(labels["SingleR"]) && labels["SingleR"] != "") {
    return(normalize_cell_type(labels["SingleR"]))
  }
  # otherwise find common words they agree on
  cw <- common_words_label(unname(lbls))
  if (!is.null(cw) && nchar(cw) > 0) return(normalize_cell_type(cw))
  # final fallback: first non-NA label in prefer order
  for (m in prefer) {
    if (!is.null(labels[m]) && !is.na(labels[m]) && labels[m] != "") return(normalize_cell_type(labels[m]))
  }
  normalize_cell_type(unname(lbls)[1])
}

# Apply consensus to clusters of a Seurat object
apply_consensus_to_clusters <- function(seu,
                                        singleR_cluster_labels,
                                        methodB_cluster_labels = NULL,
                                        methodC_cluster_labels = NULL,
                                        cluster_key = "seurat_clusters") {
  # Build data frame of labels per cluster
  clusters <- sort(unique(seu[[cluster_key]][,1]))
  lab_df <- data.frame(
    cluster = clusters,
    SingleR = singleR_cluster_labels[as.character(clusters)],
    MethodB = if (!is.null(methodB_cluster_labels)) methodB_cluster_labels[as.character(clusters)] else NA_character_,
    MethodC = if (!is.null(methodC_cluster_labels)) methodC_cluster_labels[as.character(clusters)] else NA_character_,
    stringsAsFactors = FALSE
  )
  # Compute consensus per cluster
  lab_df$consensus <- apply(lab_df[, c("SingleR","MethodB","MethodC")], 1, function(row) {
    r <- setNames(as.character(row), c("SingleR","MethodB","MethodC"))
    consensus_label(r)
  })
  # ensure factor levels reflect unique consensus names
  lab_df$consensus <- vapply(lab_df$consensus, normalize_cell_type, character(1))
  # Map clusters to consensus labels
  map <- setNames(lab_df$consensus, lab_df$cluster)
  # Create a new factor ident with consensus names
  new_ids <- map[as.character(seu[[cluster_key]][,1])]
  seu$consensus_ident <- factor(new_ids, levels = unique(lab_df$consensus))
  return(list(seu = seu, label_table = lab_df))
}

# Plotting with consensus labels
save_consensus_plots <- function(seu, outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  # UMAP labeled by consensus
  p_umap <- DimPlot(seu, reduction = "umap", group.by = "consensus_ident", label = TRUE, repel = TRUE) +
    ggplot2::ggtitle("UMAP (Consensus Annotation)")
  ggplot2::ggsave(file.path(outdir, "consensus_umap.pdf"), p_umap, width = 9, height = 7)
  # Dot plot example (adjust features list to your marker set)
  # If you have a marker list per consensus group, provide it; else use default markers
  features <- unique(head(rownames(seu), 20)) # placeholder; replace with your markers
  p_dot <- DotPlot(seu, features = features, group.by = "consensus_ident") +
    ggplot2::ggtitle("Summary Dot Plot (Consensus)")
  ggplot2::ggsave(file.path(outdir, "consensus_summary_dot_plot.pdf"), p_dot, width = 12, height = 6)
  # Proportion plot by sample using consensus labels
  df_prop <- seu@meta.data %>%
    dplyr::count(sample = .data$orig.ident, group = .data$consensus_ident) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = n / sum(n))
  p_prop <- ggplot2::ggplot(df_prop, ggplot2::aes(x = sample, y = prop, fill = group)) +
    ggplot2::geom_col(position = "fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::ggtitle("Cell Type Proportions (Consensus)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(outdir, "consensus_proportion_plot.pdf"), p_prop, width = 9, height = 7)
  invisible(NULL)
}

# Example driver (adjust inputs to your environment)
run_consensus_annotation <- function(seu,
                                     singleR_cluster_labels,
                                     methodB_cluster_labels = NULL,
                                     methodC_cluster_labels = NULL,
                                     outdir = NULL) {
  # default to annotated subdir under plots_root
  if (is.null(outdir) || !nzchar(outdir)) {
    outdir <- output_dirs$annotated_plots
  }
  res <- apply_consensus_to_clusters(
    seu = seu,
    singleR_cluster_labels = singleR_cluster_labels,
    methodB_cluster_labels = methodB_cluster_labels,
    methodC_cluster_labels = methodC_cluster_labels,
    cluster_key = "seurat_clusters"
  )
  seu2 <- res$seu
  # Save plots to new PDF filenames with consensus names
  save_consensus_plots(seu2, outdir = outdir)
  return(res$label_table)
}

# Global object for sharing between steps
seurat_objects <- NULL

# Step execution function
execute_step <- function(step) {
  switch(step,
         read_rds = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
         },
         singleR = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_singleR(seurat_objects$Joined_TN.combined)
         },
         markers = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           plot_markers(seurat_objects$Joined_TN.combined)
         },
         celliD = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_celliD(seurat_objects$TN.combined)
         },
         scCATCH = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
         },
         consensus = {
           # Expect precomputed per-cluster labels; user can adapt these sources
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           # ...existing code to build singleR/methodB/methodC label vectors...
           # Example placeholders (replace with real mappings)
           singleR_labels <- NULL  # named vector by cluster id
           methodB_labels <- NULL  # named vector by cluster id
           methodC_labels <- NULL  # named vector by cluster id

           if (is.null(singleR_labels)) stop("Provide singleR_cluster_labels mapping for consensus step.")
           run_consensus_annotation(
             seu = seurat_objects$TN.combined,
             singleR_cluster_labels = singleR_labels,
             methodB_cluster_labels = methodB_labels,
             methodC_cluster_labels = methodC_labels,
             outdir = output_dirs$annotated_plots
           )
         },
         annotated_plots = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           generate_annotated_plots(seurat_objects$TN.combined)
         },
         combined_plots = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           generate_combined_plots(seurat_objects$TN.combined)
         },
         all = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
           run_singleR(seurat_objects$Joined_TN.combined)
           plot_markers(seurat_objects$Joined_TN.combined)
           run_celliD(seurat_objects$TN.combined)
           run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
           if (opt$consensus) {
             generate_consensus_annotation()
           }
         },
         stop("Invalid step. Valid options: read_rds, singleR, markers, celliD, scCATCH, consensus, annotated_plots, combined_plots, all")
  )
}

# Main Execution ------------------------------------------------------------
execute_step(opt$step)
message("\n============================================================")
message("Cell annotation pipeline step '", opt$step, "' completed successfully")
message("============================================================")

# Close log file
message("Completed at: ", Sys.time())
sink(type = "message")
sink(type = "output")
close(log_conn)
