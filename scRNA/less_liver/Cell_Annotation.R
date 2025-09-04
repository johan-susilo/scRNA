#!/usr/bin/env Rscript
#usage: ./CellAnnotation.R -rds $PATH_to_RDS &> log
#Progressive plasticity during colorectal cancer metastasis
#A single-cell atlas of liver metastases of colorectal cancer reveals reprogramming of the tumor microenvironment in response to preoperative chemotherapy
Sys.time()

library(Seurat)
library(scCATCH)
library(SingleR)
library(celldex)
library(dplyr)
library(tidyverse)
library(CelliD)
library(ggpubr) #library for plotting


args <- commandArgs(trailingOnly = TRUE)

#create dir for output
dir.create("annotations")



# Check if any arguments are provided
if (length(args) < 1) {
  stop("Please provide sample names and identifiers as arguments.")
}

# Initialize variables and parse arguments
rds <- NULL


if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "-rds") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        rds <- strsplit(args[i + 1], ",")[[1]]
        options(rds = rds)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list of rds file-path after -rds")
      }
    } else if (args[i] == "-ident1") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        ident1 <- strsplit(args[i + 1], ",")[[1]]
        options(ident1 = ident1)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list for -ident1")
      }
    } else if (args[i] == "-ident2") {
      # Check if the next argument exists and is not another flag
      if (i + 1 <= length(args) && substr(args[i + 1], 1, 1) != "-") {
        ident2 <- strsplit(args[i + 1], ",")[[1]]
        options(ident2 = ident2)
        i <- i + 1
      } else {
        stop("Please provide a comma-separated list for -ident2")
      }
    }
    i <- i + 1
  }
}

# Use getOption to retrieve the options set via command-line
rds <- getOption("rds")
ident1 <- getOption("ident1")
ident2 <- getOption("ident2")

# Check if sample names are set, then process each
if (!is.null(rds)) {
  # Loop over each rds
  for(i in seq_along(rds)) {
  rds <- rds[i]
  rds_path <- paste0("./",rds)
  print(rds_path)
  ##load the merged Seurat object
  TN.combined <- readRDS(file = rds_path)
  DefaultAssay(TN.combined) <- "RNA"  
  Joined_TN.combined <- JoinLayers(TN.combined)
  print("------------done reading rds---------------")

  #-----------------------------------------------SingleR------------------------------------------------------------------

  print("------------start SingleR---------------")
  counts <- GetAssayData(Joined_TN.combined)

  #hpca database
  hpca.se  <- HumanPrimaryCellAtlasData()
  hpca.se
  pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1,labels = hpca.se$label.main)
  pred.hpca
  table(pred.hpca$pruned.labels)
  clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  clustering.table_hpca
  write.csv(clustering.table_hpca , file = "./annotations/SingleR_hpca.csv", col.names = TRUE)

  #bpe database
  bpe.se <- BlueprintEncodeData()
  bpe.se
  pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1,labels = bpe.se$label.main)
  pred.bpe
  table(pred.bpe$pruned.labels)
  clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  clustering.table_bpe
  write.csv(clustering.table_bpe , file = "./annotations/SingleR_bpe.csv", col.names = TRUE)

  #load hpca analysis
  clustering.table_hpca <- read.csv("./annotations/SingleR_hpca.csv")
  rownames(clustering.table_hpca) <- clustering.table_hpca[,1]
  clustering.table_hpca <- clustering.table_hpca[,-1]
  clustering.table_hpca["annotation",] <- rownames(clustering.table_hpca)[apply(clustering.table_hpca,2,which.max)]
  # Add an empty tab in front of the first row
  #colnames(clustering.table_hpca) <- c("", colnames(clustering.table_hpca))  # Add empty column at the beginning
  write.table(clustering.table_hpca , file = "./annotations/SingleR_hpca_summary.tsv", col.names = TRUE, sep= "\t", row.names = TRUE, quote = FALSE)

  #load bpe analysis
  clustering.table_bpe <- read.csv("./annotations/SingleR_bpe.csv")
  rownames(clustering.table_bpe) <- clustering.table_bpe[,1]
  clustering.table_bpe <- clustering.table_bpe[,-1]
  clustering.table_bpe["annotation",] <- rownames(clustering.table_bpe)[apply(clustering.table_bpe,2,which.max)]
  # Add an empty tab in front of the first row
  #colnames(clustering.table_hpca) <- c("", colnames(clustering.table_hpca))  # Add empty column at the beginning  
  write.table(clustering.table_bpe , file = "./annotations/SingleR_bpe_summary.tsv", sep= "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  DefaultAssay(TN.combined) <- "RNA"

  ##Classical markers
  #Epithelial cells
  markers.to.plot <- c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP")
  Classical_markers_Epi <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Epithelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Epi.pdf", width = 15, height = 15)
  print(Classical_markers_Epi)
  dev.off()

  #Sweat gland cell
  markers.to.plot <- c("MUCL1","PIP","AQP5")
  Classical_markers_Sweat_gland <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Sweat gland cell") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_Sweat_gland.pdf", width = 15, height = 15)
  Classical_markers_Sweat_gland
  dev.off()

  #SMC and pili muscle
  markers.to.plot <- c("MCAM","ACTA2","MYL9","TAGLN","MYH11")
  Classical_markers_SMC <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="SMC and pili muscle") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_SMC.pdf", width = 15, height = 15)
  print(Classical_markers_SMC)
  dev.off()

  #peicyte
  markers.to.plot <- c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11")
  Classical_markers_peicyte <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Peicyte") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_peicyte.pdf", width = 15, height = 15)
  print(Classical_markers_peicyte)
  dev.off()

  #fibroblasts
  markers.to.plot <- c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248")
  Classical_markers_fibroblasts <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Fibroblasts") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_fibroblasts.pdf", width = 15, height = 15)
  print(Classical_markers_fibroblasts)
  dev.off()

  #Vascular endothelial cells
  markers.to.plot <- c("PECAM1","VWF")
  Classical_markers_vEC <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Vascular endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_vEC.pdf", width = 15, height = 15)
  print(Classical_markers_vEC)
  dev.off()

  #lymphatic endothelial cells
  markers.to.plot <- c("PROX1","LYVE1")
  Classical_markers_lEC <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Lymphatic endothelial cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_lEC.pdf", width = 15, height = 15)
  print(Classical_markers_lEC)
  dev.off()

  #T cells
  markers.to.plot <- c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7")
  Classical_markers_TC <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="T cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_TC.pdf", width = 15, height = 15)
  print(Classical_markers_TC)
  dev.off()

  #NK cells
  markers.to.plot <- c("GNLY","NKG7")
  Classical_markers_NK <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="NK cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plots
  pdf("./annotations/Classical_markers_NK.pdf", width = 15, height = 15)
  print(Classical_markers_NK)
  dev.off()

  #B cells
  markers.to.plot <- c("MS4A1","CD79A","SEC11C","CD79B")
  Classical_markers_BC <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="B cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_BC.pdf", width = 15, height = 15)
  print(Classical_markers_BC)
  dev.off()

  #Plasma cells
  markers.to.plot <- c("IGJ","MZB1","XBP1","CD79A","CD79B")
  Classical_markers_Plasma <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Plasma cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Plasma.pdf", width = 15, height = 15)
  print(Classical_markers_Plasma)
  dev.off()

  #Myeloid cells (monocyte)
  markers.to.plot <- c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ")
  Classical_markers_monocyte <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Monocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_monocyte.pdf", width = 15, height = 15)
  print(Classical_markers_monocyte)
  dev.off()

  #Myeloid cells (dendritic cell)
  markers.to.plot <- c("IRF7","HLA-DRA","LYZ","S100B","CD1C")
  Classical_markers_dendritic <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Dendritic cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Dendritic cells.pdf", width = 15, height = 15)
  print(Classical_markers_dendritic)
  dev.off()

  #Granulocytes&Neutrophils
  markers.to.plot <- c("ITGAX","ITGAM","FCGR2A","ANPEP")
  Classical_markers_Neutrophils <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Granulocytes and Neutrophils") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Neutrophils.pdf", width = 15, height = 15)
  print(Classical_markers_Neutrophils)
  dev.off()

  #Mast cells
  markers.to.plot <- c("ADCYAP1","CPA3","TPSAB1","VWA5A")
  Classical_markers_Mast <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Mast cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Mast.pdf", width = 15, height = 15)
  print(Classical_markers_Mast)
  dev.off()

  #Melanocytes
  markers.to.plot <- c("DCT","MLANA")
  Classical_markers_Melanocytes <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Melanocytes") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Melanocytes.pdf", width = 15, height = 15)
  print(Classical_markers_Melanocytes)
  dev.off()

  #Neuronal cells
  markers.to.plot <- c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR")
  Classical_markers_Neuronal <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Neuronal cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Neuronal.pdf", width = 15, height = 15)
  print(Classical_markers_Neuronal)
  dev.off()

  #Schwann cells
  markers.to.plot <- c("NRXN1","CCN3","MPZ","PTN","S100B")
  Classical_markers_Schwann <- DotPlot(Joined_TN.combined, features = markers.to.plot, cols = c("white", "darkred"), dot.scale = 8) +
  RotatedAxis() + labs(title="Schwann cells") + theme(plot.title = element_text(hjust = 0.5, size=24))
  #output marker plot
  pdf("./annotations/Classical_markers_Schwann.pdf", width = 15, height = 15)
  print(Classical_markers_Schwann)
  dev.off()
  print("--------------------------SingleR done!--------------------------------\n")

  #-------------------------------------------SingleR end-----------------------------------------------------------

  #-------------------------------------------CelliD----------------------------------------------------------------
  print("--------------------------Start CelliD!--------------------------------\n")
  #CelliD dimensionality reduction through MCA
  DefaultAssay(TN.combined) <- "RNA"
  Joined_TN.combined = JoinLayers(TN.combined)
  Baron <- RunMCA(Joined_TN.combined,slot="RNA")
  #download all cell-type gene signatures from panglaoDB
  panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
  #filter to get human specific genes
  panglao_all <- panglao %>%  filter(str_detect(species,"Hs"))
  # convert dataframes to a list of named vectors which is the format for CelliD input
  panglao_all <- panglao_all %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))
  all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)
  #remove very short signatures
  all_gs <- all_gs[sapply(all_gs, length) >= 10]
  #RunCellHGT
  HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
  all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
  Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")
  #Output data
  pdf("./annotations/Baron_dimplot.pdf", width = 15, height = 15)
  print(DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",label = TRUE, label.size = 3,repel = TRUE)+theme(legend.text = element_text(size = 7), aspect.ratio = 1))
  clustering.table_CelliD <- table(Baron@meta.data[["all_gs_prediction_signif"]], Baron@active.ident)
  clustering.table_CelliD
  write.csv(clustering.table_CelliD , file = "./annotations/CelliD_PanglaoDB.csv")

  #load CelliD analysis
  clustering.table_CelliD <- read.csv("./annotations/CelliD_PanglaoDB.csv")
  rownames(clustering.table_CelliD) <- clustering.table_CelliD[,1]
  clustering.table_CelliD <- clustering.table_CelliD[-nrow(clustering.table_CelliD),-1]
  clustering.table_CelliD["annotation",] <- rownames(clustering.table_CelliD)[apply(clustering.table_CelliD,2,which.max)]
  # Add an empty tab in front of the first row
  #colnames(clustering.table_hpca) <- c("", colnames(clustering.table_hpca))  # Add empty column at the beginning
  write.table(clustering.table_CelliD , file = "./annotations/CelliD_PanglaoDB_summary_default.tsv", col.names = TRUE,sep= "\t", row.names =TRUE,quote = FALSE)
  
  print("--------------------------CelliD Done!--------------------------------\n")
  print("--------------------------Start scCatch!--------------------------------\n")
  # revise gene symbols
  data.input <- GetAssayData(Joined_TN.combined, assay = "RNA", layer = "data") # normalized data matrix
  data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)
  #create scCATCH object with createscCATCH(). Users need to provide the normalized data and the cluster for each cell.
  labels <- Idents(TN.combined)
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  obj <- createscCATCH(data = data.input, cluster = as.character(meta$group))
  # demo_geneinfo
  demo_marker()
  # The most strict condition to identify marker genes
  #for the tissue part might need to modify due to the sample type
  obj <- findmarkergene(object = obj,  species = "Human",    marker = cellmatch,
                      tissue = c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage','Subcutaneous adipose tissue',
                                 'Hair follicle','Lung','Muscle','Skin','Dermis','Lymph node','Lymphoid tissue',
                                 'Pluripotent stem cell','Skeletal muscle','Umbilical cord blood','Plasma',
                                 'Umbilical cord','Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell','Kidney',
                                 'Pancreas','Pancreatic islet','Pyloric gland','Pancreatic acinar tissue' ),  use_method = "1")
  obj <- findcelltype(object = obj)
  obj@celltype
  write.csv(obj@celltype, file = "./annotations/scCATCH.csv", col.names = TRUE)

  input_file <- "./annotations/scCATCH.csv"  # Replace with your actual file name
  data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
# Check and correct column names, assuming the first column is empty
  if (colnames(data)[1] == "") {
     colnames(data)[1] <- "Empty"  # Temporarily name the empty column
     }
     # Rename columns for easier reference
     colnames(data) <- c("Empty", "Cluster", "Marker", "Cell_Type", "Score", "Related_Marker", "PMID")
     # Normalize the Cell_Type column (replace underscores with spaces)
     data$Cell_Type <- gsub("_", " ", data$Cell_Type)

     # Update Cluster values: Add "X" prefix only to non-empty and non-NA values
     unique_mapping <- data %>%
     mutate(Cluster = ifelse(is.na(Cluster) | Cluster == "", "", paste0("X", Cluster))) %>%  # Prefix 'X' only if valid
     select(Cluster, Cell_Type) %>%
     distinct() %>%
     arrange(Cluster)  # Sort by Cluster
 
 # Transform the data into a wide format
 wide_data <- unique_mapping %>%
 pivot_wider(names_from = Cluster, values_from = Cell_Type)

 # Write the transformed data to a file
 output_file <- "./annotations/scCATCH_summary.tsv"
 write.table(wide_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

 # Print success message
 print(paste("File transformed and saved to", output_file))

 print("--------------------------scCatch Done!--------------------------------\n")

}
  } else {
    print("No valid rds given. Please use -rds to specify rds file.")
}   
