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

output_dirs <- list(
  singleR = file.path(output_base, "singleR"),
  markers = file.path(output_base, "markers"),
  celliD = file.path(output_base, "celliD"),
  scCATCH = file.path(output_base, "scCATCH")
)

lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Pipeline Functions --------------------------------------------------------
read_rds <- function(rds_path) {
  message("Reading RDS file: ", rds_path)
  TN.combined <- readRDS(file = rds_path)
  DefaultAssay(TN.combined) <- "RNA"
  Joined_TN.combined <- JoinLayers(TN.combined)
  message("Done reading RDS")
  return(list(TN.combined = TN.combined, Joined_TN.combined = Joined_TN.combined))
}




# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL, 
              help = "Path to RDS file"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Pipeline step: read_rds, singleR, markers, celliD, scCATCH, all"),
  make_option(c("-o", "--output"), type = "character", default = "annotations",
              help = "Base output directory [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Update output directory if specified
if (opt$output != "annotations") {
  output_base <- opt$output
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
  
  output_dirs <- list(
    singleR = file.path(output_base, "singleR"),
    markers = file.path(output_base, "markers"),
    celliD = file.path(output_base, "celliD"),
    scCATCH = file.path(output_base, "scCATCH")
  )
  
  lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
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
           run_scCATCH(seurat_objects$TN.combined)
         },
         all = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
           run_singleR(seurat_objects$Joined_TN.combined)
           plot_markers(seurat_objects$Joined_TN.combined)
           run_celliD(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
           run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
         },
         stop("Invalid step. Valid options: read_rds, singleR, markers, celliD, scCATCH, all")
  )
}

# Main Execution ------------------------------------------------------------
execute_step(opt$step)
message("Cell annotation pipeline step ", opt$step, " completed successfully")
