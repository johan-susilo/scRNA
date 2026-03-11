library(Seurat)
library(dplyr)
library(ggplot2)

# Load your full dataset
# Replace 'your_pmh_data.rds' with your actual filename
pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

# quick check to ensure it loaded
head(pmh_obj@meta.data)

table(Idents(pmh_obj))

colnames(pmh_obj@meta.data)
