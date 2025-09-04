#!/usr/bin/env Rscript
Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
})

# Centralized Configuration -------------------------------------------------
output_base <- "/mnt/80T/johan/output/liver_2021"
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

output_dirs <- list(
  processed = file.path(output_base, "processed"),
  plots = file.path(output_base, "plots"),
  tables = file.path(output_base, "tables")
)

lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Pipeline Functions --------------------------------------------------------
read_samples_tsv <- function(tsv_file) {
  samples_df <- read.delim(tsv_file, header = FALSE, 
                          col.names = c("sample_names", "ident1", "ident2"),
                          stringsAsFactors = FALSE)
  saveRDS(samples_df, file.path(output_base, "samples_df.rds"))
  samples_df
}

process_sample <- function(sample_name, ident1, ident2) {
  output_rds <- file.path(output_dirs$processed, paste0(sample_name, "_processed.rds"))
  
  if(file.exists(output_rds)) {
    message("Loading preprocessed: ", sample_name)
    return(readRDS(output_rds))
  }
  
  message("\nProcessing sample: ", sample_name)
  data_dir <- file.path("/mnt/80T/johan/data/liver_R/data", sample_name)
  
  # Original processing logic preserved
  seur_obj <- CreateSeuratObject(
    counts = Read10X(data.dir = data_dir),
    project = sample_name,
    min.cells = 3,
    min.features = 10
  )
  
  # Preserved analysis pipeline
  seur_obj <- seur_obj %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  
  # DoubletFinder processing (original logic)
  sweep.res.list <- paramSweep(seur_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  # Preserved doublet removal logic
  annotations <- seur_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.08*nrow(seur_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK, 
                           nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Preserved QC and filtering
  seur_obj$orig.ident1 <- ident1
  seur_obj$orig.ident2 <- ident2
  seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")
  seur_obj <- subset(seur_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
  
  saveRDS(seur_obj, output_rds)
  seur_obj
}

integrate_samples <- function(sample_list) {
  integrated_rds <- file.path(output_base, "TN.combined_dim30.rds")
  
  if(file.exists(integrated_rds)) {
    message("Loading pre-integrated data")
    return(readRDS(integrated_rds))
  }
  
  message("Integrating samples...")
  TN.anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:30)
  TN.combined <- IntegrateData(anchorset = TN.anchors, dims = 1:30) %>%
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  
  saveRDS(TN.combined, integrated_rds)
  TN.combined
}

generate_plots <- function(TN_combined) {
  message("Generating visualizations...")
  
  # UMAP plots
  pdf(file.path(output_dirs$plots, "TNcombined_umap_labelT.pdf"), width = 15, height = 15)
  print(DimPlot(TN_combined, reduction = "umap", label = TRUE))
  dev.off()

  # Get valid cluster IDs
  valid_clusters <- levels(Idents(TN_combined))
  
  # Heatmap generation with proper formatting
  Heatmapall <- subset(TN_combined, idents = valid_clusters)
  
  Heatmapall.markers <- FindAllMarkers(
    Heatmapall,
    only.pos = TRUE,
    min.pct = 0.25,  # Fixed decimal format
    logfc.threshold = 0.25
  )
  
  top10 <- Heatmapall.markers %>% 
    group_by(cluster) %>% 
    top_n(n = 10, wt = avg_log2FC)
  
  pdf(file.path(output_dirs$plots, "heatmap_top10.pdf"), width = 25, height = 25)
  print(DoHeatmap(Heatmapall, features = top10$gene))
  dev.off()
}



# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "Input TSV file"),
  make_option(c("-s", "--step"), type = "character", 
              help = "Pipeline step: read_tsv, process, integrate, plot, all")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

execute_step <- function(step) {
  switch(step,
    read_tsv = {
      read_samples_tsv(opt$file)
    },
    process = {
      samples_df <- readRDS(file.path(output_base, "samples_df.rds"))
      lapply(1:nrow(samples_df), function(i) {
        process_sample(samples_df$sample_names[i],
                      samples_df$ident1[i],
                      samples_df$ident2[i])
      })
    },
    integrate = {
      sample_list <- list.files(output_dirs$processed, 
                               pattern = "_processed.rds$", 
                               full.names = TRUE) %>%
        lapply(readRDS)
      integrate_samples(sample_list)
    },
    plot = {
      TN_combined <- readRDS(file.path(output_base, "TN.combined_dim30.rds"))
      generate_plots(TN_combined)
    },
    all = {
      execute_step("read_tsv")
      execute_step("process")
      execute_step("integrate")
      execute_step("plot")
    },
    stop("Invalid step. Valid options: read_tsv, process, integrate, plot, all")
  )
}

# Main Execution ------------------------------------------------------------
if(is.null(opt$file)) stop("Must specify input file with -f")
if(is.null(opt$step)) stop("Must specify execution step with -s")

execute_step(opt$step)
message("Pipeline step '", opt$step, "' completed successfully")
