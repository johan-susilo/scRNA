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
})

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to RDS file (TN.combined_dim30.rds)"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Pipeline step: read_rds, singleR, markers, celliD, scCATCH, consensus, all"),
  make_option(c("-o", "--output"), type = "character", default = "annotations",
              help = "Base output directory [default: annotations]"),
  make_option(c("--consensus"), action = "store_true", default = FALSE,
              help = "Generate consensus annotations from all methods"),
  make_option(c("--tissue"), type = "character", default = "skin",
              help = "Tissue type for scCATCH [default: skin]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Use provided output path and create directories
output_base <- opt$output
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

output_dirs <- list(
  singleR = file.path(output_base, "singleR"),
  markers = file.path(output_base, "markers"),
  celliD = file.path(output_base, "celliD"),
  scCATCH = file.path(output_base, "scCATCH"),
  consensus = file.path(output_base, "consensus"),
  logs = file.path(output_base, "logs")
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

# Pipeline Functions --------------------------------------------------------
read_rds <- function(rds_path) {
  message("============================================================")
  message("Reading RDS file: ", rds_path)
  message("============================================================")
  TN.combined <- readRDS(file = rds_path)
  DefaultAssay(TN.combined) <- "RNA"
  Joined_TN.combined <- JoinLayers(TN.combined)
  message("Done reading RDS")
  message("Cells: ", ncol(Joined_TN.combined), ", Features: ", nrow(Joined_TN.combined))
  message("Clusters: ", length(unique(Idents(Joined_TN.combined))))
  return(list(TN.combined = TN.combined, Joined_TN.combined = Joined_TN.combined))
}

#-----------------------------------------------SingleR------------------------------------------------------------------
run_singleR <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Starting SingleR Annotation")
  message("============================================================")
  counts <- GetAssayData(Joined_TN.combined)

  # hpca database
  message("Running SingleR with HumanPrimaryCellAtlas database...")
  hpca.se <- tryCatch({
    HumanPrimaryCellAtlasData()
  }, error = function(e) {
    message("Error loading HumanPrimaryCellAtlas: ", e$message)
    message("Attempting workaround...")
    # Clear cache and retry
    eh <- ExperimentHub::ExperimentHub()
    celldex::HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = "nonna")
  })

  pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_hpca, file = file.path(output_dirs$singleR, "SingleR_hpca.csv"), col.names = TRUE)

  # bpe database
  message("Running SingleR with BlueprintEncode database...")
  bpe.se <- tryCatch({
    BlueprintEncodeData()
  }, error = function(e) {
    message("Error loading BlueprintEncode: ", e$message)
    message("Attempting workaround...")
    # Clear cache and retry
    eh <- ExperimentHub::ExperimentHub()
    celldex::BlueprintEncodeData(ensembl = FALSE, cell.ont = "nonna")
  })

  pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1, labels = bpe.se$label.main)
  clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_bpe, file = file.path(output_dirs$singleR, "SingleR_bpe.csv"), col.names = TRUE)

  # Process and summarize results
  message("Processing SingleR results...")
  clustering.table_hpca <- read.csv(file.path(output_dirs$singleR, "SingleR_hpca.csv"))
  rownames(clustering.table_hpca) <- clustering.table_hpca[,1]
  clustering.table_hpca <- clustering.table_hpca[,-1]
  clustering.table_hpca["annotation",] <- rownames(clustering.table_hpca)[apply(clustering.table_hpca,2,which.max)]
  write.table(clustering.table_hpca, file = file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"),
              col.names = TRUE, sep= "\t", row.names = TRUE, quote = FALSE)

  clustering.table_bpe <- read.csv(file.path(output_dirs$singleR, "SingleR_bpe.csv"))
  rownames(clustering.table_bpe) <- clustering.table_bpe[,1]
  clustering.table_bpe <- clustering.table_bpe[,-1]
  clustering.table_bpe["annotation",] <- rownames(clustering.table_bpe)[apply(clustering.table_bpe,2,which.max)]
  write.table(clustering.table_bpe, file = file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"),
              sep= "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  message("SingleR annotation completed!")
  message("============================================================\n")
  return(list(hpca = clustering.table_hpca, bpe = clustering.table_bpe))
}
#-------------------------------------------SingleR end-----------------------------------------------------------

#-------------------------------------------Classical Markers------------------------------------------------------
plot_markers <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Generating Classical Marker Plots")
  message("============================================================")

  # Define marker sets for different cell types 
  marker_sets <- list(
    "Epithelial" = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP"),
    "Sweat_gland" = c("MUCL1","PIP","AQP5"),
    "SMC" = c("MCAM","ACTA2","MYL9","TAGLN","MYH11"),
    "Pericyte" = c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11"),
    "Fibroblasts" = c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248"),
    "Vascular_EC" = c("PECAM1","VWF"),
    "Lymphatic_EC" = c("PROX1","LYVE1"),
    "T_cells" = c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7"),
    "NK_cells" = c("GNLY","NKG7"),
    "B_cells" = c("MS4A1","CD79A","SEC11C","CD79B"),
    "Plasma_cells" = c("IGJ","MZB1","XBP1","CD79A","CD79B"),
    "Monocytes" = c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ"),
    "Dendritic_cells" = c("IRF7","HLA-DRA","LYZ","S100B","CD1C"),
    "Neutrophils" = c("ITGAX","ITGAM","FCGR2A","ANPEP"),
    "Mast_cells" = c("ADCYAP1","CPA3","TPSAB1","VWA5A"),
    "Melanocytes" = c("DCT","MLANA"),
    "Neuronal_cells" = c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR"),
    "Schwann_cells" = c("NRXN1","CCN3","MPZ","PTN","S100B")
  )

  message("Processing ", length(marker_sets), " marker sets...")

  # Process each marker set sequentially
  results <- lapply(names(marker_sets), function(cell_type) {
    markers <- marker_sets[[cell_type]]
    plot_title <- gsub("_", " ", cell_type)

    p <- DotPlot(Joined_TN.combined, features = markers, cols = c("white", "darkred"), dot.scale = 8) +
      RotatedAxis() + labs(title=plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size=24))

    plot_file <- file.path(output_dirs$markers, paste0("Classical_markers_", cell_type, ".pdf"))
    tryCatch({
      pdf(plot_file, width = 15, height = 15)
      print(p)
      message("Saved marker plot: ", plot_file)
    }, error = function(e) {
      message("Warning: Failed to save marker plot for ", cell_type, ": ", conditionMessage(e))
    }, finally = {
      if (length(dev.list()) > 0) dev.off()
    })

    return(cell_type)
  })

  message("Marker plots generated for ", length(results), " cell types")
  message("============================================================\n")
}


#-------------------------------------------CelliD----------------------------------------------------------------
run_celliD <- function(seurat_object) {
  message("\n============================================================")
  message("Starting CelliD Annotation")
  message("============================================================")

  # Step 1: Downsample the object to a manageable size to prevent memory errors.
  if (ncol(seurat_object) > 90000) {
    message("Object has > 90,000 cells. Downsampling to 90,000 for CelliD analysis.")
    seurat_object_subset <- subset(seurat_object, cells = sample(Cells(seurat_object), 90000))
  } else {
    seurat_object_subset <- seurat_object
  }

  # Step 2: Join layers ON THE SMALL, SUBSETTED OBJECT.
  message("Joining layers on the downsampled object...")
  seurat_subset_joined <- JoinLayers(seurat_object_subset)

  # Step 3: CelliD dimensionality reduction through MCA
  message("Running MCA (Multiple Correspondence Analysis)...")
  Baron <- RunMCA(seurat_subset_joined, slot = "RNA")

  # Download all cell-type gene signatures from PanglaoDB
  message("Downloading PanglaoDB gene signatures...")
  panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", show_col_types = FALSE)

  # Filter to get human specific genes
  panglao_all <- panglao %>% filter(str_detect(species,"Hs"))

  # Convert dataframes to a list of named vectors which is the format for CelliD input
  panglao_all <- panglao_all %>%
    group_by(`cell type`) %>%
    summarise(geneset = list(`official gene symbol`))

  all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)

  # Remove very short signatures
  all_gs <- all_gs[sapply(all_gs, length) >= 10]

  message("Running RunCellHGT for cell type prediction...")
  HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
  all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
  Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")

  # Safe plotting with device cleanup
  safe_save_pdf(
    DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",
            label = TRUE, label.size = 3, repel = TRUE) +
      theme(legend.text = element_text(size = 7), aspect.ratio = 1),
    file.path(output_dirs$celliD, "Baron_dimplot.pdf")
  )

  # Output data
  message("Summarizing CelliD results...")
  clustering.table_CelliD <- table(Baron@meta.data[["all_gs_prediction_signif"]], Baron@active.ident)
  write.csv(clustering.table_CelliD, file = file.path(output_dirs$celliD, "CelliD_PanglaoDB.csv"))

  table_for_summary <- clustering.table_CelliD
  if ("unassigned" %in% rownames(table_for_summary)) {
    table_for_summary <- table_for_summary[rownames(table_for_summary) != "unassigned", , drop = FALSE]
  }

  get_annotation <- function(col) {
    if (all(col == 0)) {
      return("unassigned")
    } else {
      return(rownames(table_for_summary)[which.max(col)])
    }
  }

  annotation_row <- apply(table_for_summary, 2, get_annotation)
  summary_df <- as.data.frame.matrix(clustering.table_CelliD)
  summary_df["annotation",] <- annotation_row

  write.table(summary_df, file = file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"),
              col.names = TRUE, sep= "\t", row.names =TRUE, quote = FALSE)

  message("CelliD annotation completed!")
  message("============================================================\n")
  return(Baron)
}

#-------------------------------------------scCATCH----------------------------------------------------------------
run_scCATCH <- function(TN.combined, Joined_TN.combined) {
  message("\n============================================================")
  message("Starting scCATCH Annotation")
  message("============================================================")

  # Get normalized data matrix
  data.input <- GetAssayData(Joined_TN.combined, assay = "RNA", layer = "data")

  # Revise gene symbols
  message("Revising gene symbols...")
  data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)

  # Create scCATCH object
  labels <- Idents(TN.combined)
  meta <- data.frame(group = labels, row.names = names(labels))
  obj <- createscCATCH(data = data.input, cluster = as.character(meta$group))

  # Define tissue types based on input
  tissue_list <- if (opt$tissue == "skin") {
    c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage','Subcutaneous adipose tissue',
      'Hair follicle','Lung','Muscle','Skin','Dermis','Lymph node','Lymphoid tissue',
      'Pluripotent stem cell','Skeletal muscle','Umbilical cord blood','Plasma',
      'Umbilical cord','Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell','Kidney',
      'Pancreas','Pancreatic islet','Pyloric gland','Pancreatic acinar tissue')
  } else {
    c('Blood','Peripheral blood','Lymph node','Lymphoid tissue','Bone marrow','Spleen')
  }

  message("Finding marker genes with tissue filter: ", opt$tissue)
  obj <- findmarkergene(object = obj,
                        species = "Human",
                        marker = cellmatch,
                        tissue = tissue_list,
                        use_method = "1")

  # Find cell types
  message("Finding cell types...")
  obj <- findcelltype(object = obj)

  # Save results
  write.csv(obj@celltype, file = file.path(output_dirs$scCATCH, "scCATCH.csv"), col.names = TRUE)

  # Process scCATCH output for consensus (similar to Cell_anno_consensus.pl)
  message("Processing scCATCH results for consensus...")
  data <- read.csv(file.path(output_dirs$scCATCH, "scCATCH.csv"), header = TRUE, stringsAsFactors = FALSE)

  if (colnames(data)[1] == "") {
    colnames(data)[1] <- "Empty"
  }
  colnames(data) <- c("Empty", "Cluster", "Marker", "Cell_Type", "Score", "Related_Marker", "PMID")

  # Normalize cell types
  data$Cell_Type <- sapply(data$Cell_Type, normalize_cell_type)

  # Create summary with cluster prefix
  unique_mapping <- data %>%
    mutate(Cluster = ifelse(is.na(Cluster) | Cluster == "", "", paste0("X", Cluster))) %>%
    select(Cluster, Cell_Type) %>%
    distinct() %>%
    arrange(Cluster)

  wide_data <- unique_mapping %>%
    pivot_wider(names_from = Cluster, values_from = Cell_Type)

  write.table(wide_data, file = file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("scCATCH annotation completed!")
  message("============================================================\n")
  return(obj)
}
#-------------------------------------------scCATCH end-----------------------------------------------------------

#-------------------------------------------Consensus Annotation---------------------------------------------------
generate_consensus_annotation <- function() {
  message("\n============================================================")
  message("Generating Consensus Annotations")
  message("============================================================")

  # Initialize storage for annotations by cluster
  cell_annotation <- list()

  # Read SingleR HPCA results
  if (file.exists(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"))) {
    message("Reading SingleR HPCA annotations...")
    hpca_data <- read.delim(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"),
                           header = TRUE, stringsAsFactors = FALSE, sep = "\t")

    if ("annotation" %in% rownames(hpca_data)) {
      annotation_row <- hpca_data[rownames(hpca_data) == "annotation", ]
      for (col_name in colnames(annotation_row)) {
        if (col_name != "" && !is.na(annotation_row[[col_name]]) && annotation_row[[col_name]] != "") {
          cluster_key <- paste0("X", col_name)
          normalized_type <- normalize_cell_type(annotation_row[[col_name]])
          if (!cluster_key %in% names(cell_annotation)) {
            cell_annotation[[cluster_key]] <- list()
          }
          if (!normalized_type %in% names(cell_annotation[[cluster_key]])) {
            cell_annotation[[cluster_key]][[normalized_type]] <- 0
          }
          cell_annotation[[cluster_key]][[normalized_type]] <- cell_annotation[[cluster_key]][[normalized_type]] + 1
        }
      }
    }
  }

  # Read SingleR BPE results
  if (file.exists(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"))) {
    message("Reading SingleR BPE annotations...")
    bpe_data <- read.delim(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"),
                          header = TRUE, stringsAsFactors = FALSE, sep = "\t")

    if ("annotation" %in% rownames(bpe_data)) {
      annotation_row <- bpe_data[rownames(bpe_data) == "annotation", ]
      for (col_name in colnames(annotation_row)) {
        if (col_name != "" && !is.na(annotation_row[[col_name]]) && annotation_row[[col_name]] != "") {
          cluster_key <- paste0("X", col_name)
          normalized_type <- normalize_cell_type(annotation_row[[col_name]])
          if (!cluster_key %in% names(cell_annotation)) {
            cell_annotation[[cluster_key]] <- list()
          }
          if (!normalized_type %in% names(cell_annotation[[cluster_key]])) {
            cell_annotation[[cluster_key]][[normalized_type]] <- 0
          }
          cell_annotation[[cluster_key]][[normalized_type]] <- cell_annotation[[cluster_key]][[normalized_type]] + 1
        }
      }
    }
  }

  # Read CelliD results
  if (file.exists(file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"))) {
    message("Reading CelliD annotations...")
    cellid_data <- read.delim(file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"),
                             header = TRUE, stringsAsFactors = FALSE, sep = "\t")

    if ("annotation" %in% rownames(cellid_data)) {
      annotation_row <- cellid_data[rownames(cellid_data) == "annotation", ]
      for (col_name in colnames(annotation_row)) {
        if (col_name != "" && col_name != "X" && !is.na(annotation_row[[col_name]]) && annotation_row[[col_name]] != "unassigned") {
          cluster_key <- paste0("X", col_name)
          normalized_type <- normalize_cell_type(annotation_row[[col_name]])
          if (!cluster_key %in% names(cell_annotation)) {
            cell_annotation[[cluster_key]] <- list()
          }
          if (!normalized_type %in% names(cell_annotation[[cluster_key]])) {
            cell_annotation[[cluster_key]][[normalized_type]] <- 0
          }
          cell_annotation[[cluster_key]][[normalized_type]] <- cell_annotation[[cluster_key]][[normalized_type]] + 1
        }
      }
    }
  }

  # Read scCATCH results
  if (file.exists(file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"))) {
    message("Reading scCATCH annotations...")
    sccatch_data <- read.delim(file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"),
                              header = TRUE, stringsAsFactors = FALSE, sep = "\t")

    for (row_idx in 1:nrow(sccatch_data)) {
      for (col_name in colnames(sccatch_data)) {
        cell_type_str <- sccatch_data[row_idx, col_name]
        if (!is.na(cell_type_str) && cell_type_str != "") {
          # Split multi-type annotations
          individual_types <- strsplit(cell_type_str, ",")[[1]]
          for (cell_type in individual_types) {
            cell_type <- trimws(cell_type)
            cluster_key <- col_name
            normalized_type <- normalize_cell_type(cell_type)
            if (!cluster_key %in% names(cell_annotation)) {
              cell_annotation[[cluster_key]] <- list()
            }
            if (!normalized_type %in% names(cell_annotation[[cluster_key]])) {
              cell_annotation[[cluster_key]][[normalized_type]] <- 0
            }
            cell_annotation[[cluster_key]][[normalized_type]] <- cell_annotation[[cluster_key]][[normalized_type]] + 1
          }
        }
      }
    }
  }

  # Generate consensus by finding max count for each cluster
  message("Computing consensus annotations...")
  consensus_results <- data.frame(
    Cluster = character(),
    Cell_Type = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )

  for (cluster in sort(names(cell_annotation))) {
    max_count <- 0
    max_cell_types <- c()

    for (cell_type in names(cell_annotation[[cluster]])) {
      count <- cell_annotation[[cluster]][[cell_type]]
      if (count > max_count) {
        max_count <- count
        max_cell_types <- c(cell_type)
      } else if (count == max_count) {
        max_cell_types <- c(max_cell_types, cell_type)
      }
    }

    cell_types_str <- paste(max_cell_types, collapse = ", ")
    consensus_results <- rbind(consensus_results,
                              data.frame(Cluster = cluster,
                                       Cell_Type = cell_types_str,
                                       Count = max_count,
                                       stringsAsFactors = FALSE))
  }

  # Save consensus results
  write.table(consensus_results,
              file = file.path(output_dirs$consensus, "consensus_annotation.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("Consensus annotation completed!")
  message("Results saved to: ", file.path(output_dirs$consensus, "consensus_annotation.tsv"))
  message("============================================================\n")

  return(consensus_results)
}
#-------------------------------------------Consensus Annotation end-----------------------------------------------


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
           generate_consensus_annotation()
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
         stop("Invalid step. Valid options: read_rds, singleR, markers, celliD, scCATCH, consensus, all")
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
