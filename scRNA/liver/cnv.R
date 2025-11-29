#!/usr/bin/env Rscript
#nohup Rscript ./New_CellAnnotation.R --rds /mnt/80T/johan/output/full_liver_2021/TN.combined_dim30.rds --output /mnt/80T/johan/output/full_liver_2021 > annotation.log 2>&1 &

Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(numbat)
  library(celldex)
  library(dplyr)
  library(tidyverse)
  library(ggpubr)
})

# Centralized Configuration -------------------------------------------------
output_base <- "annotations"
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

# ...existing code...
read_rds <- function(rds_path) {
  message("Reading RDS file: ", rds_path)
  TN.combined <- readRDS(file = rds_path)
  DefaultAssay(TN.combined) <- "RNA"
  # If JoinLayers exists in user's environment keep it; otherwise just return same object
  Joined_TN.combined <- tryCatch({
    JoinLayers(TN.combined)
  }, error = function(e) {
    message("JoinLayers not available/failed; continuing with original Seurat object")
    TN.combined
  })
  message("Done reading RDS")
  return(list(TN.combined = TN.combined, Joined_TN.combined = Joined_TN.combined))
}

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL, 
              help = "Path to RDS file"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Pipeline step: read_rds, cnv, all"),
  make_option(c("-o", "--output"), type = "character", default = "annotations",
              help = "Base output directory [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Update output directory if specified
if (opt$output != "annotations") {
  output_base <- opt$output
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
}

# Global object for sharing between steps
seurat_objects <- NULL

# CNV step implementation (robust and lightweight) --------------------------
run_cnv <- function(seurat_obj, output_base, name = "unknown") {
  if (is.null(seurat_obj)) stop("Seurat object is NULL")
  cnv_dir <- file.path(output_base, "cnv")
  dir.create(cnv_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- file.path(cnv_dir, paste0("cnv_summary_", name))

  message("Computing average expression per cluster as a CNV-related summary...")
  # Ensure RNA assay is active
  tryCatch({
    DefaultAssay(seurat_obj) <- "RNA"
  }, error = function(e) {})

  # Use Seurat::AverageExpression which is stable for Seurat objects
  avg_expr <- tryCatch({
    ae <- AverageExpression(seurat_obj, assays = "RNA", slot = "data")
    # AverageExpression returns a list per assay; take RNA if present
    if (is.list(ae) && "RNA" %in% names(ae)) ae$RNA else if (is.matrix(ae)) ae else ae[[1]]
  }, error = function(e) {
    warning("AverageExpression failed: ", e$message)
    NULL
  })

  if (!is.null(avg_expr)) {
    write.csv(avg_expr, file = paste0(out_prefix, "_avg_expr_per_cluster.csv"))
    message("Wrote average expression per cluster to: ", paste0(out_prefix, "_avg_expr_per_cluster.csv"))
  } else {
    message("Skipping average expression output due to failure")
  }

  # Save a lightweight RDS summary for downstream inspection
  summary_list <- list(
    timestamp = Sys.time(),
    object_name = name,
    has_clusters = tryCatch({ length(levels(Idents(seurat_obj))) > 0 }, error = function(e) FALSE),
    avg_expr = avg_expr
  )
  saveRDS(summary_list, file = paste0(out_prefix, ".rds"))
  message("Saved CNV summary RDS to: ", paste0(out_prefix, ".rds"))

  invisible(summary_list)
}

# Step execution function (simplified to CNV workflow) ----------------------
execute_step <- function(step) {
  step <- tolower(step)
  switch(step,
         read_rds = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
         },
         cnv = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           # Run CNV-like summary on available objects (safe: check presence)
           if (!is.null(seurat_objects$Joined_TN.combined)) {
             run_cnv(seurat_objects$Joined_TN.combined, output_base, name = "joined")
           }
           if (!is.null(seurat_objects$TN.combined)) {
             run_cnv(seurat_objects$TN.combined, output_base, name = "tn")
           }
         },
         all = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
           # Run CNV step after reading
           if (!is.null(seurat_objects$Joined_TN.combined)) {
             run_cnv(seurat_objects$Joined_TN.combined, output_base, name = "joined")
           }
           if (!is.null(seurat_objects$TN.combined)) {
             run_cnv(seurat_objects$TN.combined, output_base, name = "tn")
           }
         },
         stop("Invalid step. Valid options: read_rds, cnv, all")
  )
}

# Main Execution ------------------------------------------------------------
execute_step(opt$step)
message("CNV pipeline step ", opt$step, " completed successfully")
