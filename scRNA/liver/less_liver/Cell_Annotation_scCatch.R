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
