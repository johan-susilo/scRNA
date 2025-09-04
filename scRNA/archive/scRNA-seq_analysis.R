# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         

setwd("/mnt/80T/johan/liver_R")
# load libraries
library(Seurat)
library(tidyverse)


# Load the dataset file
if (file.exists("data_seurat_object.rds")) {
    data.seurat <- readRDS("data_seurat_object.rds") #when you want to load the data from the file
    print("Seurat object loaded successfully")
} else {
    data <- Read10X(data.dir = "data/")
    data.seurat <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 10)
    saveRDS(data.seurat, file = "data_seurat_object.rds")
}


data.seurat 
#25122 features across 140281 samples


# 1. QC -------

# % MT reads
data.seurat[["percent.mt"]] <- PercentageFeatureSet(data.seurat, pattern = "^MT-")
View(data.seurat@meta.data)



# 2. Filtering -----------------
data.seurat <- subset(data.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                          percent.mt < 5)

# make violin plot
pdf("violin_plot.pdf", width = 15, height = 15)
VlnPlot(data.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() #close the pdf file

# make scatter plot
pdf("featureScatter_plot.pdf", width = 15, height = 15)
FeatureScatter(data.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
dev.off() #close the pdf file


# 3. Normalize data ----------
data.seurat <- NormalizeData(data.seurat)
#str(data.seurat)


# 4. Identify highly variable features --------------
data.seurat <- FindVariableFeatures(data.seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data.seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data.seurat)
pdf("high_variable_feature.pdf", width = 15, height = 15)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off() #close the pdf file

# 5. Scaling -------------
all.genes <- rownames(data.seurat)
data.seurat <- ScaleData(data.seurat, features = all.genes)

str(data.seurat)

# 6. Perform Linear dimensionality reduction --------------
data.seurat <- RunPCA(data.seurat, features = VariableFeatures(object = data.seurat))


# visualize PCA results
print(data.seurat[["pca"]], dims = 1:5, nfeatures = 5)
pdf("pca_plot.pdf", width = 15, height = 15)
DimHeatmap(data.seurat, dims = 1, cells = 500, balanced = TRUE)
dev.off() #close the pdf file

# determine dimensionality of the data
pdf("elbow_plot.pdf", width = 15, height = 15)
ElbowPlot(data.seurat)
dev.off() #close the pdf file

# 7. Clustering ------------
data.seurat <- FindNeighbors(data.seurat, dims = 1:15)

# understanding resolution
data.seurat <- FindClusters(data.seurat, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(data.seurat@meta.data)

DimPlot(data.seurat, group.by = "RNA_snn_res.0.1", label = TRUE)

# setting identity of clusters
Idents(data.seurat)
Idents(data.seurat) <- "RNA_snn_res.0.1"
Idents(data.seurat)

# non-linear dimensionality reduction --------------
data.seurat <- RunUMAP(data.seurat, dims = 1:15)


pdf("umap_plot.pdf", width = 15, height = 15)
DimPlot(data.seurat, reduction = "umap")
dev.off()

saveRDS(data.seurat, file = "data_seurat_object_complete.rds")

