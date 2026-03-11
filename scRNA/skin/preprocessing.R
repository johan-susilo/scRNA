#!/usr/bin/env Rscript
# Usage: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -s all -o /home/johan/output/skin_pmh/

# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /path/to/data -s read_csv
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /path/to/data -s process -c 4
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d . -s integrate
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /home/johan/data/PMH_scRNA-seq -s plot
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -s all -o /home/johan/output/skin_pmh -c 8
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+K -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+K
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task2.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+SS -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+SS
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task3.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+PMH -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+PMH
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task4.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+HS -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+HS
#Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task5.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+DFonly1 -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+DFonly1



Sys.time()
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(ggsci)
  library(Matrix)
  library(ggpubr)
  library(cowplot)
  library(gridExtra)
  # library(clusterProfiler)  # Not used in this pipeline, commented out to avoid dependency issues
  library(gplots)
  library(ggplot2)
  library(ggnewscale)
  library(RColorBrewer)
  library(tidyr)
  library(presto)
  library(ggrepel)
  library(stringr)
  library(patchwork)
  library(scales)
  library(parallel)
  library(future)
  library(future.apply)
  library(scCATCH)
})

# --- Increase Global Memory Limit for Parallel Processing ---
options(future.globals.maxSize = 100 * 1024^3)

# --- Helper Function for QC Violin Plots (workaround for Seurat v5.3.1 VlnPlot bug) ---
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
      colors <- colorRampPalette(brewer.pal(min(n_clusters, 12), "Paired"))(n_clusters)
    } else {
      colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
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

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "Input CSV file (columns: sample_names, ident1, ident2)"),
  make_option(c("-d", "--datadir"), type = "character", default = ".",
              help = "Base directory containing sample folders [default: current directory]"),
  make_option(c("-s", "--step"), type = "character",
              help = "Pipeline step: read_csv, process, integrate, plot, all"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Base output directory [overrides default output_base]"),
  make_option(c("-c", "--cores"), type = "integer", default = NULL,
              help = "Number of cores for parallel processing [default: detect available cores]"),
  make_option(c("--doublet_rate"), type = "numeric", default = 0.08,
              help = "Expected doublet formation rate [default: 0.08 for ~10k cells]"),
  make_option(c("--min_features"), type = "integer", default = 200,
              help = "Minimum number of features per cell [default: 200]"),
  make_option(c("--max_features"), type = "integer", default = 5000,
              help = "Maximum number of features per cell [default: 5000]"),
  make_option(c("--max_mt"), type = "numeric", default = 30,
              help = "Maximum mitochondrial percentage [default: 30]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Provide a sensible default output base and directories
output_base <- "output"
output_dirs <- list(
  processed = file.path(output_base, "processed"),
  plots = file.path(output_base, "plots"),
  tables = file.path(output_base, "tables"),
  logs = file.path(output_base, "logs")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# If -o provided, override output_base and recreate directories
if (!is.null(opt$output) && nzchar(opt$output)) {
  output_base <- opt$output
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
  output_dirs <- list(
    processed = file.path(output_base, "processed"),
    plots = file.path(output_base, "plots"),
    tables = file.path(output_base, "tables"),
    logs = file.path(output_base, "logs")
  )
  lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}

# Set up logging
log_file <- file.path(output_dirs$logs, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)  # split=TRUE means output to both console and file
sink(log_conn, type = "message")

message("Log file: ", log_file)
message("Started at: ", Sys.time())

# Set up processing mode - ALWAYS use sequential to avoid deadlocks
message("Using SEQUENTIAL processing mode to avoid deadlocks")
plan(sequential)

# Define color schemes
mycolor <- c("#CCCCCC", # gray for other
  "#A6CEE3","#FF7F00","#09519C","#FFFB33","#556b2f",
  "#FF00FF","#377EB8","#bb8fce","#666666","#90EE90",
  "#ff4500","#A6761D","#E67E22","#323695","#E81E32",
  "#006837","#CBC3E3","#F1C40F","#3498DB","#34495E",
  "#FA9399","#48C9B0","#7D3C98","#ff4500","#8b4513",
  "#8a2be2","#f0e68c","#00ffff","#32CD32","#b03060")

doublet_color <- c("Doublet" = "#e35473", "Singlet" = "#54cdeb")

read_samples_csv <- function(csv_file) {
  Sys.time()
  # Read CSV file with sample information
  message("Reading sample information from ", csv_file)

  if (!file.exists(csv_file)) {
    stop("Input CSV file not found: ", csv_file, call. = FALSE)
  }
  fi <- file.info(csv_file)
  if (is.na(fi$size) || fi$size == 0) {
    stop("Input CSV exists but is empty: ", csv_file,
         "\nPlease provide a CSV with a header and at least one sample row.", call. = FALSE)
  }

  samples_df <- tryCatch({
    read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Failed to read CSV file: ", csv_file, "\nError: ", e$message, call. = FALSE)
  })

  if (nrow(samples_df) == 0) {
    stop("Input CSV has a header but no data rows. Add at least one sample row to ", csv_file, call. = FALSE)
  }

  # Verify expected columns are present
  required_cols <- c("sample_names", "ident1", "ident2")
  if (!all(required_cols %in% colnames(samples_df))) {
    stop("CSV file must contain columns: ", paste(required_cols, collapse=", "), call. = FALSE)
  }

  # Save the dataframe for future steps
  saveRDS(samples_df, file.path(output_base, "samples_df.rds"))

  message("Found ", nrow(samples_df), " samples in the CSV file")
  return(samples_df)
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

process_sample <- function(sample_name, sample_ident1, sample_ident2, base_data_dir, out_dirs) {
  # Wrap the ENTIRE function in tryCatch so one bad sample doesn't kill the parallel run
  Sys.time()
  tryCatch({
    output_rds <- file.path(out_dirs$processed, paste0(sample_name, "_processed.rds"))
    sample_id <- sample_name

    if(file.exists(output_rds)) {
      message("Loading preprocessed: ", sample_name)
      return(output_rds)  # Return path, not object
    }

    message("\n==================== Processing sample: ", sample_name, " ====================")
    data_dir <- file.path(base_data_dir, sample_name)

    # Create sample-specific plot directories
    sample_plot_dir <- file.path(out_dirs$plots, sample_name)
    dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)

    # Create a Seurat object
    seur_obj <- CreateSeuratObject(
      counts = Read10X(data.dir = data_dir),
      project = sample_name,
      min.cells = 3,
      min.features = 10
    )
    message("Created Seurat object for ", sample_name, " with dimensions: ", dim(seur_obj)[1], " features and ", dim(seur_obj)[2], " cells")

    # ============ QC METRICS (BEFORE ANY FILTERING) ============
    message("Calculating QC metrics...")
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")

    # QC violin plots BEFORE filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_01_qc_violins_unfiltered.pdf"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Unfiltered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (unfiltered) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ PRE-FILTERING (BEFORE DOUBLET FINDER) ============
    # This significantly speeds up DoubletFinder by removing obvious empty droplets
    message("Pre-filtering to remove empty droplets...")

    # Apply basic pre-filtering to remove obvious junk
    cells_before_prefilter <- ncol(seur_obj)
    seur_obj <- subset(seur_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 50)
    cells_after_prefilter <- ncol(seur_obj)

    message("Pre-filtering removed ", cells_before_prefilter - cells_after_prefilter,
            " empty droplets (", round((cells_before_prefilter - cells_after_prefilter) / cells_before_prefilter * 100, 2),
            "%), keeping ", cells_after_prefilter, " cells for DoubletFinder")

    # QC violin plots AFTER pre-filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_02_qc_violins_prefiltered.pdf"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Pre-filtered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (pre-filtered) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ DOUBLET FINDER WORKFLOW  ============
    message("Starting DoubletFinder workflow...")

    # Normalize, find variable features, scale data, run PCA
    seur_obj <- seur_obj %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA()

    # Plot elbow plot
    safe_save_pdf(ElbowPlot(seur_obj),
                  file.path(sample_plot_dir, paste0(sample_id, "_03_elbow.pdf")))

    # Find neighbors and clusters (needed for DoubletFinder)
    seur_obj <- seur_obj %>%
      FindNeighbors(dims = 1:50) %>%
      FindClusters() %>%
      RunUMAP(dims = 1:30)

    # Plot initial UMAP
    safe_save_pdf(DimPlot(seur_obj, reduction = "umap", label = TRUE),
                  file.path(sample_plot_dir, paste0(sample_id, "_04_umap_initial.pdf")))

    # pK Identification (no ground-truth) - optimal parameter for DoubletFinder
    message("Identifying optimal pK parameter...")
    sweep.res.list <- paramSweep(seur_obj, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    safe_save_pdf(ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
                    geom_point() + geom_line() +
                    ggtitle(paste0("pK Identification - ", sample_name)) +
                    theme_bw(),
                  file.path(sample_plot_dir, paste0(sample_id, "_05_doubletfinder_pk.pdf")))

    # Select optimal pK
    pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
    pK <- as.numeric(as.character(pK[[1]]))
    message("Optimal pK: ", pK)

    # Homotypic Doublet Proportion Estimate
    annotations <- seur_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(opt$doublet_rate * nrow(seur_obj@meta.data))  # Use configurable doublet rate
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    message("Expected doublets: ", nExp_poi, " (adjusted: ", nExp_poi.adj, ")")

    # Identify doublet cells
    message("Running DoubletFinder...")
    seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK,
                              nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

    # Dynamic column detection (critical fix from previous version)
    DF_cols <- grep("DF.classifications", colnames(seur_obj@meta.data), value = TRUE)

    if (length(DF_cols) == 0) {
      stop("DoubletFinder failed to add classification column for sample: ", sample_name)
    }
    DF_classification <- DF_cols[1]
    message("DoubletFinder classification column: ", DF_classification)

    # Plot unfiltered doublet results
    safe_save_pdf(DimPlot(seur_obj, reduction = 'umap', group.by = DF_classification, cols = doublet_color) +
                    ggtitle(paste0("Doublets (Before Removal) - ", sample_name)),
                  file.path(sample_plot_dir, paste0(sample_id, "_06_doublets_unfiltered.pdf")))

    # Check doublet statistics
    doublet_table <- table(seur_obj@meta.data[[DF_classification]])
    message("Doublet statistics:")
    print(doublet_table)

    # Remove doublet cells (keep only Singlets)
    singlet_indices <- which(seur_obj@meta.data[[DF_classification]] == "Singlet")
    seur_obj <- seur_obj[, singlet_indices]
    message("Cells after doublet removal: ", ncol(seur_obj))

    # Plot filtered doublet results
    safe_save_pdf(DimPlot(seur_obj, reduction = "umap", group.by = DF_classification, cols = doublet_color) +
                    ggtitle(paste0("Doublets (After Removal) - ", sample_name)),
                  file.path(sample_plot_dir, paste0(sample_id, "_07_doublets_filtered.pdf")))

    # QC violin plots AFTER doublet removal
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_08_qc_violins_post_doublet.pdf"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Post-Doublet Removal) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (post-doublet) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ QUALITY CONTROL  ============
    message("Performing final quality control filtering...")

    # Add metadata identifiers
    seur_obj$orig.ident1 <- sample_ident1
    seur_obj$orig.ident2 <- sample_ident2

    # Quality filtering using configurable parameters
    message("Filtering cells with criteria: nFeature_RNA > ", opt$min_features,
            " & < ", opt$max_features, ", percent.mt < ", opt$max_mt)
    cells_before_final_qc <- ncol(seur_obj)
    seur_obj <- subset(seur_obj, subset = nFeature_RNA > opt$min_features &
                                          nFeature_RNA < opt$max_features &
                                          percent.mt < opt$max_mt)
    cells_after_final_qc <- ncol(seur_obj)
    message("Final QC filtering removed ", cells_before_final_qc - cells_after_final_qc,
            " cells, keeping ", cells_after_final_qc, " high-quality cells")

    # QC violin plots AFTER final filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_09_qc_violins_final.pdf"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Final Filtered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (final) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # Save processed object
    saveRDS(seur_obj, output_rds)
    message("Processed object saved: ", output_rds)
    message("==================== Completed: ", sample_name, " ====================\n")

    # Return the file path instead of the object to avoid serialization issues in parallel
    return(output_rds)

  }, error = function(e) {
    # If anything fails inside this function, catch it here so the other parallel workers don't die
    message("\nCRITICAL ERROR processing sample ", sample_name, ": ", e$message)
    return(NULL) # Return NULL to indicate failure
  })
}


integrate_samples <- function(sample_list) {
  Sys.time()
  # Remove NULLs from failed samples
  sample_list <- Filter(Negate(is.null), sample_list)

  if (length(sample_list) == 0) {
    stop("No valid samples available for integration.")
  }

  message("\n==================== Starting Integration ====================")

  # Cell cycle scoring and normalization 
  message("Starting normalization and cell cycle scoring for ", length(sample_list), " samples...")
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  TN.list <- lapply(X = sample_list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x, features = rownames(x), verbose = TRUE)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    return(x)
  })

  # Integration workflow
  message("Selecting integration features...")
  features <- SelectIntegrationFeatures(object.list = TN.list)

  # Check if integration anchors already exist
  anchors_rds <- file.path(output_base, "TN.anchorsdim30.rds")

  if (file.exists(anchors_rds)) {
    message("Loading existing integration anchors from: ", anchors_rds)
    TN.anchors <- readRDS(anchors_rds)
    message("Integration anchors loaded successfully")
  } else {
    message("Finding integration anchors (this may take a while)...")
    TN.anchors <- FindIntegrationAnchors(object.list = TN.list, dims = 1:30, verbose = TRUE)
    saveRDS(TN.anchors, anchors_rds)
    message("Integration anchors saved to: ", anchors_rds)
  }

  message("Integrating data...")
  TN.combined <- IntegrateData(anchorset = TN.anchors, dims = 1:30, verbose = TRUE)
  DefaultAssay(TN.combined) <- "integrated"
  message("Data integration completed")

  # Scaling with regression (from Yan's approach - regress out cell cycle and MT%)
  message("Scaling integrated data with regression...")
  TN.combined <- ScaleData(TN.combined,
                          vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                          features = rownames(TN.combined), verbose = TRUE)

  message("Running PCA...")
  TN.combined <- RunPCA(TN.combined, npcs = 50, verbose = TRUE)

  # Save elbow plot
  safe_save_pdf(ElbowPlot(TN.combined, ndims = 50),
                file.path(output_dirs$plots, "TNcombined_elbow.pdf"))

  message("Finding neighbors...")
  TN.combined <- FindNeighbors(TN.combined, reduction = "pca", dims = 1:30)

  message("Finding clusters...")
  TN.combined <- FindClusters(TN.combined, resolution = 0.5)

  message("Running UMAP...")
  TN.combined <- RunUMAP(TN.combined, reduction = "pca", dims = 1:30)

  # Save tables
  message("Saving cell count tables...")
  CellNumber <- table(Idents(TN.combined), TN.combined$orig.ident1)
  write.csv(CellNumber, file = file.path(output_dirs$tables, "CellNumber_bygroup.csv"))

  # Save integrated object
  message("Saving integrated object...")
  saveRDS(TN.combined, file = file.path(output_base, "TN.combined_dim30.rds"))
  message("==================== Integration complete! ====================\n")

  return(TN.combined)
}


create_summary_dot_plot <- function(seurat_object, output_dir) {
  Sys.time()
  message("Generating summary marker gene dot plot...")

  # Define the exact list of marker genes
  summary_marker_genes <- c(
    # Immune - T-cells / NK
    "PTPRC", "CD3D", "CD3E", "CD8A", "GZMA", "GZMB", "KLRD1", "NKG7",
    # Immune - B-cells / Plasma
    "CD79A", "MS4A1", "IGHG1", "JCHAIN",
    # Immune - Myeloid
    "LYZ", "CD68", "CD163", "CD14", "FCGR3A", "MS4A7",
    # Epithelial Cells
    "EPCAM", "KRT8", "KRT18",
    # Endothelial Cells
    "PECAM1", "VWF",
    # Fibroblasts
    "DCN", "LUM", "COL1A1",
    # Smooth Muscle Cells / Pericytes
    "ACTA2", "TAGLN", "RGS5",
    # Melanocytes
    "MLANA", "PMEL",
    # Mast Cells
    "TPSAB1", "CPA3"
  )

  # Ensure the features are present in the object to avoid errors
  genes_in_object <- rownames(seurat_object)
  plot_genes <- intersect(summary_marker_genes, genes_in_object)

  if (length(plot_genes) == 0) {
    warning("None of the specified marker genes were found in the Seurat object.")
    return(NULL)
  }

  # Create the DotPlot
  p <- DotPlot(seurat_object, features = plot_genes) +
    RotatedAxis() +
    ggtitle("Summary Marker Gene Expression") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))

  # Save the plot to a PDF file
  ggsave(file.path(output_dir, "summary_dot_plot.pdf"), plot = p, width = 18, height = 10)

  message("Summary dot plot saved to: ", file.path(output_dir, "summary_dot_plot.pdf"))

  print(p)
}

generate_plots <- function(TN.combined) {
  Sys.time()
  message("\n==================== Generating Plots ====================")

  # --- 1. Attempt to Apply Annotations ---
  # Look for the consensus file in the standard location relative to output_base
  consensus_file <- file.path(output_base, "annotations","consensus", "consensus_annotation.tsv")
  
  # If output_base was customized in opt, check there too
  if (!file.exists(consensus_file) && !is.null(opt$output)) {
     consensus_file <- file.path(opt$output, "annotations", "consensus", "consensus_annotation.tsv")
  }

  TN.plotting <- TN.combined
  is_annotated <- FALSE

  if (file.exists(consensus_file)) {
    message("Found consensus annotation file. Applying labels...")
    consensus_data <- read.delim(consensus_file, sep = "\t", stringsAsFactors = FALSE)
    
    # Ensure Cluster format matches (remove X)
    consensus_data$Cluster <- as.character(gsub("^X", "", consensus_data$Cluster))
    
    current_ids <- levels(TN.combined)
    new_names <- character(length(current_ids))
    names(new_names) <- current_ids
    
    for (id in current_ids) {
      clean_id <- as.character(gsub("^X", "", id))
      match_row <- consensus_data[consensus_data$Cluster == clean_id, ]
      
      if (nrow(match_row) > 0) {
        # Format: "0: Cell Type" (Preserves unique colors for all 20 clusters)
        new_names[id] <- paste0(clean_id, ": ", match_row$Cell_Type[1])
      } else {
        new_names[id] <- paste0(clean_id, ": Unknown")
      }
    }
    
    TN.plotting <- RenameIdents(TN.combined, new_names)
    is_annotated <- TRUE
  } else {
    message("Consensus file not found. Using numeric cluster labels.")
  }

  # Define colors 
  # If annotated, we might have more than 30 levels if something went wrong, 
  # but usually it's the same count. Safe color generation:
  n_groups <- length(levels(TN.plotting))
  if (n_groups <= length(mycolor)) {
    plot_colors <- mycolor[1:n_groups]
  } else {
    plot_colors <- colorRampPalette(mycolor)(n_groups)
  }

  # --- 2. UMAP Plots ---
  # Note: label.size set to 3 to accommodate longer text labels
  
  safe_save_pdf(DimPlot(TN.plotting, reduction = "umap", label = FALSE, pt.size = 0.8, cols = plot_colors),
                file.path(output_dirs$plots, "TNcombined_umap_labelF.pdf"))

  safe_save_pdf(DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.8, cols = plot_colors),
                file.path(output_dirs$plots, "TNcombined_umap_labelT.pdf"))

  safe_save_pdf(DimPlot(TN.plotting, group.by = "orig.ident2", pt.size = 0.8, cols = plot_colors),
                file.path(output_dirs$plots, "TNcombined_umap_groupbyorigident.pdf"))

  # Split by Condition (ident1) - Matches your uploaded image
  safe_save_pdf(DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, split.by = "orig.ident1",
                        pt.size = 0.8, ncol = 2, cols = plot_colors),
                file.path(output_dirs$plots, "TNcombined_umap_labelT_splitorigident1.pdf"))

  # Split by Sample (ident2)
  safe_save_pdf(DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, split.by = "orig.ident2",
                        pt.size = 0.8, ncol = 2, cols = plot_colors),
                file.path(output_dirs$plots, "TNcombined_umap_labelT_splitsample.pdf"))

  # --- 3. Heatmap ---
  # Using numeric clusters for Heatmap is usually cleaner for calculation, 
  # but we use the plotting object to match the UMAP labels.
  CellNumber <- table(Idents(TN.plotting), TN.plotting$orig.ident1)
  cluster_count <- nrow(CellNumber)

  message("Finding markers for heatmap...")
  # We subset using indices to get all clusters
  Heatmapall <- subset(TN.plotting, idents = levels(TN.plotting))
  Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  write.csv(Heatmapall.markers, file = file.path(output_dirs$tables, "Findallmarkers.csv"), row.names = FALSE)

  top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

  safe_save_pdf(DoHeatmap(Heatmapall, features = top10$gene, group.colors = plot_colors),
                file.path(output_dirs$plots, "heatmap_top10.pdf"),
                w = 25, h = 25)

  # --- 4. Proportion Plots ---
  # Matches your uploaded proportion_plot.jpg
  
  # By Condition (ident1)
  Cellproportion <- table(Idents(TN.plotting), TN.plotting$orig.ident1)
  Cellproportion <- round(sweep(Cellproportion, MARGIN = 2, STATS = colSums(Cellproportion), FUN = "/") * 100, 2)
  write.csv(Cellproportion, file = file.path(output_dirs$tables, "Cellproportion.csv"), row.names = TRUE)

  Cellproportion_df <- as.data.frame(Cellproportion)
  
  proportion_plot <- ggplot(Cellproportion_df, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.6) +
    xlab("Sample") +
    ylab("Proportion") +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 8)) + # Smaller text for longer labels
    scale_fill_manual(values = plot_colors)

  safe_save_pdf(proportion_plot, file.path(output_dirs$plots, "proportion_plot.pdf"))

  # By Sample (ident2)
  Sampleproportion <- table(Idents(TN.plotting), TN.plotting$orig.ident2)
  Sampleproportion <- round(sweep(Sampleproportion, MARGIN = 2, STATS = colSums(Sampleproportion), FUN = "/") * 100, 2)
  write.csv(Sampleproportion, file = file.path(output_dirs$tables, "Sampleproportion.csv"), row.names = TRUE)

  Sampleproportion_df <- as.data.frame(Sampleproportion)
  colnames(Sampleproportion_df) <- c("Cluster", "Sample", "Freq")

  sample_proportion_plot <- ggplot(Sampleproportion_df, aes(x = Sample, y = Freq, fill = Cluster)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.6) +
    xlab("Sample") +
    ylab("Proportion") +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_fill_manual(values = plot_colors)

  safe_save_pdf(sample_proportion_plot, file.path(output_dirs$plots, "sample_proportion_plot.pdf"))

  message("==================== Plots completed! ====================\n")
}

run_and_plot_pca <- function(seurat_object, output_plots_dir) {
  Sys.time()
  message("Generating sample similarity plot using pseudo-bulk PCA...")

  # Create "pseudo-bulk" profiles by averaging expression for each sample (using orig.ident2)
  avg_expr <- AverageExpression(
    seurat_object,
    group.by = "orig.ident2",
    assays = "RNA",
    layer = "data"
  )

  # The result is a matrix of genes x samples. Transpose it for PCA.
  avg_expr_matrix <- t(avg_expr$RNA)

  # Run PCA on the pseudo-bulk profiles
  pca_results <- prcomp(avg_expr_matrix, scale. = TRUE)

  # Prepare the results for plotting
  pca_data <- as.data.frame(pca_results$x)
  pca_data$sample <- rownames(pca_data)

  # Calculate the percentage of variance explained by each PC
  pca_variance <- pca_results$sdev^2 / sum(pca_results$sdev^2)

  # Create the Plot
  similarity_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 6, alpha = 0.9) +
    geom_text_repel(aes(label = sample), size = 4, box.padding = 0.5) +
    theme_bw() +
    guides(color = "none") +
    ggtitle("PCA of Sample Similarity (Pseudo-bulk)") +
    labs(
      x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

  # Save the plot
  ggsave(file.path(output_plots_dir, "pca_sample_similarity_plot.png"), plot = similarity_plot, width = 10, height = 8)

  message("Sample similarity PCA plot saved.")

  print(similarity_plot)
}

execute_step <- function(step) {
  switch(step,
    read_csv = {
      read_samples_csv(opt$file)
    },
    process = {
      samples_df <- readRDS(file.path(output_base, "samples_df.rds"))
      message("\n==================== Checking Sample Processing Status ====================")

      # Check which samples already have processed RDS files
      existing_samples <- character()
      new_samples <- character()

      for (i in 1:nrow(samples_df)) {
        sample_name <- samples_df$sample_names[i]
        output_rds <- file.path(output_dirs$processed, paste0(sample_name, "_processed.rds"))

        if (file.exists(output_rds)) {
          existing_samples <- c(existing_samples, sample_name)
        } else {
          new_samples <- c(new_samples, sample_name)
        }
      }

      if (length(existing_samples) > 0) {
        message("Found ", length(existing_samples), " already processed samples:")
        message("  ", paste(existing_samples, collapse=", "))
      }

      if (length(new_samples) > 0) {
        message("Will process ", length(new_samples), " new samples:")
        message("  ", paste(new_samples, collapse=", "))
      } else {
        message("All samples already processed! Skipping processing step.")
        return(invisible(NULL))
      }

      message("\n==================== Starting Sample Processing ====================")

      # Process samples sequentially using lapply
      results <- lapply(1:nrow(samples_df), function(i) {
        process_sample(samples_df$sample_names[i],
                      samples_df$ident1[i],
                      samples_df$ident2[i],
                      opt$datadir,
                      output_dirs)
      })

      # Check for failures
      failed_indices <- which(sapply(results, is.null))
      if (length(failed_indices) > 0) {
        message("\nWARNING: The following samples FAILED to process:")
        message(paste(samples_df$sample_names[failed_indices], collapse=", "))
        message("Check the logs for details.")
      } else {
        message("All samples processed successfully")
      }
      invisible(results)
    },
    integrate = {
      # Check if integrated object already exists
      integrated_rds <- file.path(output_base, "TN.combined_dim30.rds")

      if (file.exists(integrated_rds)) {
        message("Loading existing integrated object from: ", integrated_rds)
        TN_combined <- readRDS(integrated_rds)
        message("Integrated object loaded successfully with ", ncol(TN_combined), " cells")
        return(TN_combined)
      }

      # If not, proceed with integration
      sample_files <- list.files(output_dirs$processed,
                               pattern = "_processed.rds$",
                               full.names = TRUE)

      if (length(sample_files) == 0) stop("No processed sample files found in ", output_dirs$processed)

      message("Found ", length(sample_files), " processed samples for integration")
      sample_list <- lapply(sample_files, readRDS)
      integrate_samples(sample_list)
    },
    plot = {
      TN_combined <- readRDS(file.path(output_base, "TN.combined_dim30.rds"))
      generate_plots(TN_combined)
      run_and_plot_pca(TN_combined, output_dirs$plots)
      create_summary_dot_plot(TN_combined, output_dirs$plots)
    },
    all = {
      execute_step("read_csv")
      execute_step("process")
      execute_step("integrate")
      execute_step("plot")
    },
    stop("Invalid step. Valid options: read_csv, process, integrate, plot, all")
  )
}

# Main Execution
if(is.null(opt$file)) stop("Must specify input file with -f")
if(is.null(opt$step)) stop("Must specify execution step with -s")

execute_step(opt$step)
message("Pipeline step '", opt$step, "' completed successfully")

# Close log file
message("Completed at: ", Sys.time())
sink(type = "message")
sink(type = "output")
close(log_conn)
