#!/usr/bin/env Rscript
#nohup Rscript ./New_CellAnnotation.R --rds /mnt/80T/johan/output/full_liver_2021/TN.combined_dim30.rds --output /mnt/80T/johan/output/full_liver_2021 > annotation.log 2>&1 &

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

run_singleR <- function(Joined_TN.combined) {
  message("Starting SingleR annotation")
  counts <- GetAssayData(Joined_TN.combined)
  
  # hpca database
  hpca.se <- HumanPrimaryCellAtlasData()
  pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_hpca, file = file.path(output_dirs$singleR, "SingleR_hpca.csv"), col.names = TRUE)
  
  # bpe database
  bpe.se <- BlueprintEncodeData()
  pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1, labels = bpe.se$label.main)
  clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_bpe, file = file.path(output_dirs$singleR, "SingleR_bpe.csv"), col.names = TRUE)
  
  # Process and summarize results
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
  
  message("SingleR annotation completed")
  return(list(hpca = clustering.table_hpca, bpe = clustering.table_bpe))
}

plot_markers <- function(Joined_TN.combined) {
  message("Generating marker plots")
  
  # Define marker sets
  marker_sets <- list(
    "Epi" = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP"),
    "Sweat_gland" = c("MUCL1","PIP","AQP5"),
    "SMC" = c("MCAM","ACTA2","MYL9","TAGLN","MYH11"),
    "peicyte" = c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11"),
    "fibroblasts" = c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248"),
    "vEC" = c("PECAM1","VWF"),
    "lEC" = c("PROX1","LYVE1"),
    "TC" = c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7"),
    "NK" = c("GNLY","NKG7"),
    "BC" = c("MS4A1","CD79A","SEC11C","CD79B"),
    "Plasma" = c("IGJ","MZB1","XBP1","CD79A","CD79B"),
    "monocyte" = c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ"),
    "dendritic" = c("IRF7","HLA-DRA","LYZ","S100B","CD1C"),
    "Neutrophils" = c("ITGAX","ITGAM","FCGR2A","ANPEP"),
    "Mast" = c("ADCYAP1","CPA3","TPSAB1","VWA5A"),
    "Melanocytes" = c("DCT","MLANA"),
    "Neuronal" = c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR"),
    "Schwann" = c("NRXN1","CCN3","MPZ","PTN","S100B")
  )
  
  # Process each marker set
  for(cell_type in names(marker_sets)) {
    markers <- marker_sets[[cell_type]]
    plot_title <- gsub("_", " ", cell_type)
    
    p <- DotPlot(Joined_TN.combined, features = markers, cols = c("white", "darkred"), dot.scale = 8) +
      RotatedAxis() + labs(title=plot_title) + 
      theme(plot.title = element_text(hjust = 0.5, size=24))
    
    pdf(file.path(output_dirs$markers, paste0("Classical_markers_", cell_type, ".pdf")), width = 15, height = 15)
    print(p)
    dev.off()
  }
  
  message("Marker plots generated")
}


run_celliD <- function(seurat_object) {
  message("Starting CelliD annotation")
  
  # Step 1: Downsample the object to a manageable size to prevent memory errors.
  if (ncol(seurat_object) > 90000) {
    message("Object has > 90,000 cells. Downsampling to 90,000 for CelliD analysis.")
    seurat_object_subset <- subset(seurat_object, cells = sample(Cells(seurat_object), 90000))
  } else {
    seurat_object_subset <- seurat_object
  }
  
  # Step 2: Join layers ON THE SMALL, SUBSETTED OBJECT.
  # This creates the simplified data structure that older functions expect.
  message("Joining layers on the downsampled object...")
  seurat_subset_joined <- JoinLayers(seurat_object_subset)
  
  # Step 3: Run MCA on the joined subset.
  # We pass the 'RNA' assay name to the 'slot' argument, which directly addresses the error message.
  message("Running MCA on the joined subset...")
  Baron <- RunMCA(seurat_subset_joined, slot = "RNA")
  
  
  # --- The rest of the function proceeds as before ---
  message("Downloading PanglaoDB gene signatures...")
  panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", show_col_types = FALSE)
  panglao_all <- panglao %>% filter(str_detect(species,"Hs"))
  
  panglao_all <- panglao_all %>%
    group_by(`cell type`) %>%
    summarise(geneset = list(`official gene symbol`))
  
  all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)
  all_gs <- all_gs[sapply(all_gs, length) >= 10]
  
  message("Running RunCellHGT for cell type prediction...")
  HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
  all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
  Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")
  
  pdf(file.path(output_dirs$celliD, "Baron_dimplot.pdf"), width = 15, height = 15)
  print(DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",
                label = TRUE, label.size = 3, repel = TRUE) + 
          theme(legend.text = element_text(size = 7), aspect.ratio = 1))
  dev.off()
  
  message("Summarizing CelliD results...")
  clustering.table_CelliD <- table(Baron@meta.data[["all_gs_prediction_signif"]], Baron@active.ident)
  write.csv(clustering.table_CelliD, file = file.path(output_dirs$celliD, "CelliD_PanglaoDB_full_table.csv"))
  
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
  
  message("CelliD annotation completed")
  return(Baron)
}

run_scCATCH <- function(seurat_object) {
  message("Starting scCATCH annotation")
  
  # --- Step 1: Downsample to prevent memory errors ---
  # If the object is large, we create a smaller, representative subset.
  if (ncol(seurat_object) > 100000) {
    message("Object has > 100,000 cells. Downsampling to 100,000 for scCATCH analysis.")
    seurat_subset <- subset(seurat_object, cells = sample(Cells(seurat_object), 100000))
  } else {
    seurat_subset <- seurat_object
  }
  
  # --- Step 2: Join Layers on the downsampled subset ---
  # This correctly combines data from all samples into a single layer.
  message("Joining layers on the subset...")
  seurat_subset_joined <- JoinLayers(seurat_subset)
  
  # --- Step 3: Extract the data and matching metadata ---
  # We use GetAssayData, which works on the now-joined object.
  data.input <- GetAssayData(seurat_subset_joined, assay = "RNA", layer = "data")
  
  # CRITICAL: We must get the cluster labels from the SAME subset object
  # to ensure the number of cells matches the data matrix.
  labels <- Idents(seurat_subset_joined)
  meta <- data.frame(group = labels, row.names = names(labels))
  
  # --- The rest of the function now runs correctly ---
  data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)
  
  obj <- createscCATCH(data = data.input, cluster = as.character(meta$group))
  
  tissues <- c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage','Subcutaneous adipose tissue',
               'Hair follicle','Lung','Muscle','Skin','Dermis','Lymph node','Lymphoid tissue',
               'Pluripotent stem cell','Skeletal muscle','Umbilical cord blood','Plasma',
               'Umbilical cord','Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell','Kidney',
               'Pancreas','Pancreatic islet','Pyloric gland','Pancreatic acinar tissue')
  
  obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch,
                       tissue = tissues, use_method = "1")
  
  obj <- findcelltype(object = obj)
  write.csv(obj@celltype, file = file.path(output_dirs$scCATCH, "scCATCH.csv"), col.names = TRUE)
  
  data <- read.csv(file.path(output_dirs$scCATCH, "scCATCH.csv"), header = TRUE, stringsAsFactors = FALSE)
  
  if (colnames(data)[1] == "") {
    colnames(data)[1] <- "Empty"
  }
  
  colnames(data) <- c("Empty", "Cluster", "Marker", "Cell_Type", "Score", "Related_Marker", "PMID")
  data$Cell_Type <- gsub("_", " ", data$Cell_Type)
  
  unique_mapping <- data %>%
    mutate(Cluster = ifelse(is.na(Cluster) | Cluster == "", "", paste0("X", Cluster))) %>%
    select(Cluster, Cell_Type) %>%
    distinct() %>%
    arrange(Cluster)
  
  wide_data <- unique_mapping %>%
    pivot_wider(names_from = Cluster, values_from = Cell_Type)
  
  write.table(wide_data, file = file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("scCATCH annotation completed")
  return(obj)
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
