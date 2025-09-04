#packages install
#remotes::install_version("SeuratObject", version= "5.0.1")
#remotes::install_version("Seurat", version= "5.0.2")
# The easiest way to get dplyr is to install the whole tidyverse:
#remotes::install_version("matrixStats", version= "1.1.0")
#install.packages("tidyverse")
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#install.packages("Matrix")
#install.packages("ggpubr")
#install.packages("gridExtra")
#BiocManager::install("clusterProfiler")
#install.packages('gplots')
#install.packages("ggnewscale")
#install.packages("RColorBrewer")
#

#Library
library(Seurat)
library(dplyr)
library(ggsci)
library(Matrix)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(clusterProfiler)
library(gplots)
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)
library(tidyr) 
library(DoubletFinder)
#setwd("C:/Users/user/Desktop/202403/")

#Because the need of integration, need to input all data and preprocess at once.

##Load all the dataset 
#HC131#Healthy skin
HC131.data <- Read10X(data.dir = "TN131_HTY/")
HC131 <- CreateSeuratObject(counts = HC131.data, project = "131", min.cells = 3, min.features = 10)

#HC244#Healthy skin
HC244.data <- Read10X(data.dir = "TN244_HTY/")
HC244 <- CreateSeuratObject(counts = HC244.data, project = "244", min.cells = 3, min.features = 10)

#UA258#Unaffected
UA258.data <- Read10X(data.dir = "TN258_Unaffected/")
UA258 <- CreateSeuratObject(counts = UA258.data, project = "258", min.cells = 3, min.features = 10)

#AC259#Acute skin
AC259.data <- Read10X(data.dir = "TN259_acute/")
AC259 <- CreateSeuratObject(counts = AC259.data, project = "259", min.cells = 3, min.features = 10)

#CH260#Chronic skin
CH260.data <- Read10X(data.dir = "TN260_chronic/")
CH260 <- CreateSeuratObject(counts = CH260.data, project = "260", min.cells = 3, min.features = 10)

##Re-check
#HC131 Healthy 
HC131 #16770 genes x 2457 samples
#HC244 Healthy
HC244 #21577 genes x 9753 samples
#UA258 Unaffected
UA258 #16694 gene x 2863 samples
#AC259 Acute
AC259 19818 genes x 7076 samples
#CH260 Chronic
CH260 #19188 genes x 6519 samples

##Remove doublet cells by DoubletFinder package
#HC131
#Run general flow of scRNA-seq by Seurat package
HC131 <- NormalizeData(object = HC131)
HC131 <- FindVariableFeatures(object = HC131)
HC131 <- ScaleData(object = HC131)
HC131 <- RunPCA(object = HC131)
pdf("./202408/HC131/HC131_elbow1.pdf", width = 15, height = 15)
ElbowPlot(HC131)
dev.off()
HC131 <- FindNeighbors(object = HC131, dims = 1:50)
HC131 <- FindClusters(object = HC131)
HC131 <- RunUMAP(object = HC131, dims = 1:30)
pdf("./202408/HC131/HC131_umap1.pdf", width = 15, height = 15)
DimPlot(HC131, reduction = "umap", label = T)
dev.off()
#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(HC131, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pdf("./202408/HC1331/HC131_pkplot.pdf", width = 15, height = 15)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- HC131@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(HC131@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
HC131 <- doubletFinder(HC131, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(HC131)

#group.by will changed due to different sample. ##Try to make it as a variable.
pdf("./202408/HC131/HC131_unfiltered_doublet.pdf", width = 15, height = 15)
DimPlot(HC131, reduction = 'umap', group.by = "DF.classifications_0.25_0.02_176")
dev.off()
table(HC131@meta.data$DF.classifications_0.25_0.02_176)

#Remove doublet cells
HC131 <- subset(HC131, subset = DF.classifications_0.25_0.02_176 == "Singlet")
pdf("./202408/HC131/HC131_filtered_doublet.pdf", width = 15, height = 15)
DimPlot(HC131, reduction = "umap", group.by="DF.classifications_0.25_0.02_176")
dev.off()
HC131

#HC244
#Run general flow of scRNA-seq by Seurat package
HC244 <- NormalizeData(object = HC244)
HC244 <- FindVariableFeatures(object = HC244)
HC244 <- ScaleData(object = HC244)
HC244 <- RunPCA(object = HC244)
pdf("./202408/HC244/HC244_elbow1.pdf", width = 15, height = 15)
ElbowPlot(HC244)
dev.off()
HC244 <- FindNeighbors(object = HC244, dims = 1:50)
HC244 <- FindClusters(object = HC244)
HC244 <- RunUMAP(object = HC244, dims = 1:30)
pdf("./202408/HC244/HC244_umap1.pdf", width = 15, height = 15)
DimPlot(HC244, reduction = "umap", label = T)
dev.off()

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(HC244, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pdf("./202408/HC244/HC244_pk.pdf", width = 15, height = 15)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- HC244@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(HC244@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
HC244 <- doubletFinder(HC244, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(HC244)

#group.by will changed due to different sample. ##Try to make it as a variable.
pdf("./202408/HC244/HC244_unfiltered_doublet.pdf", width = 15, height = 15)
DimPlot(HC244, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_710")
dev.off()
table(HC244@meta.data$DF.classifications_0.25_0.005_710)

#Remove doublet cells
HC244 <- subset(HC244, subset = DF.classifications_0.25_0.005_710 == "Singlet")
pdf("./202408/HC244/HC244_filtered_doublet.pdf", width = 15, height = 15)
DimPlot(HC244, reduction = "umap", group.by="DF.classifications_0.25_0.005_710")
dev.off()
HC244

#UA258
#Run general flow of scRNA-seq by Seurat package
UA258 <- NormalizeData(object = UA258)
UA258 <- FindVariableFeatures(object = UA258)
UA258 <- ScaleData(object = UA258)
UA258 <- RunPCA(object = UA258)
pdf("./202408/UA258/UA258_elbow1.pdf", width = 15, height = 15)
ElbowPlot(UA258)
dev.off()
UA258 <- FindNeighbors(object = UA258, dims = 1:50)
UA258 <- FindClusters(object = UA258)
UA258 <- RunUMAP(object = UA258, dims = 1:30)
pdf("./202408/UA258/UA258_umap1.pdf", width = 15, height = 15)
DimPlot(UA258, reduction = "umap", label = T)
dev.off()

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(UA258, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pdf("./202408/UA258/UA258_pk.pdf", width = 15, height = 15)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- UA258@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(UA258@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
UA258 <- doubletFinder(UA258, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(UA258)

#group.by will changed due to different sample. ##Try to make it as a variable.
pdf("./202408/UA258/UA258_unfiltered_doublet.pdf", width = 15, height = 15)
DimPlot(UA258, reduction = 'umap', group.by = "DF.classifications_0.25_0.24_210")
dev.off()
table(UA258@meta.data$DF.classifications_0.25_0.24_210)

#Remove doublet cells
UA258 <- subset(UA258, subset = DF.classifications_0.25_0.24_210 == "Singlet")
pdf("./202408/UA258/UA258_filtered_doublet.pdf", width = 15, height = 15)
DimPlot(UA258, reduction = "umap", group.by="DF.classifications_0.25_0.24_210")
dev.off()
UA258

#AC259
#Run general flow of scRNA-seq by Seurat package
AC259 <- NormalizeData(object = AC259)
AC259 <- FindVariableFeatures(object = AC259)
AC259 <- ScaleData(object = AC259)
AC259 <- RunPCA(object = AC259)
pdf("./202408/AC259/AC259_elbow1.pdf", width = 15, height = 15)
ElbowPlot(AC259)
dev.off()
AC259 <- FindNeighbors(object = AC259, dims = 1:50)
AC259 <- FindClusters(object = AC259)
AC259 <- RunUMAP(object = AC259, dims = 1:30)
pdf("./202408/AC259/AC259_umap1.pdf", width = 15, height = 15)
DimPlot(AC259, reduction = "umap", label = T)
dev.off()

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(AC259, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pdf("./202408/AC259/AC259_pk.pdf", width = 15, height = 15)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- AC259@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(AC259@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
AC259 <- doubletFinder(AC259, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(AC259)

#group.by will changed due to different sample. ##Try to make it as a variable.
pdf("./202408/AC259/AC259_unfiltered_doublet.pdf", width = 15, height = 15)
DimPlot(AC259, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_519")
dev.off()
table(AC259@meta.data$DF.classifications_0.25_0.3_519)

#Remove doublet cells
AC259 <- subset(AC259, subset = DF.classifications_0.25_0.3_519 == "Singlet")
pdf("./202408/AC259/AC259_filtered_doublet.pdf", width = 15, height = 15)
DimPlot(AC259, reduction = "umap", group.by="DF.classifications_0.25_0.3_519")
dev.off()
AC259

#CH260
#Run general flow of scRNA-seq by Seurat package
CH260 <- NormalizeData(object = CH260)
CH260 <- FindVariableFeatures(object = CH260)
CH260 <- ScaleData(object = CH260)
CH260 <- RunPCA(object = CH260)
pdf("./202408/CH260/CH260_elbow1.pdf", width = 15, height = 15)
ElbowPlot(CH260)
dev.off()
CH260 <- FindNeighbors(object = CH260, dims = 1:50)
CH260 <- FindClusters(object = CH260)
CH260 <- RunUMAP(object = CH260, dims = 1:30)
pdf("./202408/CH260/CH260_umap1.pdf", width = 15, height = 15)
DimPlot(CH260, reduction = "umap", label = T)
dev.off()

#pK Identification (no ground-truth)
sweep.res.list <- paramSweep(CH260, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pdf("./202408/CH260/CH260_pk.pdf", width = 15, height = 15)
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line()
dev.off()
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))
head(pK)

#Homotypic Doublet Proportion Estimate
annotations <- CH260@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.08*nrow(CH260@meta.data)) ## 8% doublet formation rate from 10X Genomics form if 10,000 cells are obtained.
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
head(nExp_poi.adj)

#Identify doublet cells
CH260 <- doubletFinder(CH260, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
head(CH260)

#group.by will changed due to different sample. ##Try to make it as a variable.
pdf("./202408/CH260/CH260_unfiltered_doublet.pdf", width = 15, height = 15)
DimPlot(CH260, reduction = 'umap', group.by = "DF.classifications_0.25_0.3_464")
dev.off()
table(CH260@meta.data$DF.classifications_0.25_0.3_464)

#Remove doublet cells
CH260 <- subset(CH260, subset = DF.classifications_0.25_0.3_464 == "Singlet")
pdf("./202408/CH260/CH260_filtered_doublet.pdf", width = 15, height = 15)
DimPlot(CH260, reduction = "umap", group.by="DF.classifications_0.25_0.3_464")
dev.off()
CH260

##Quality control
#HC131
#The $ operator can add columns to object metadata.
HC131$orig.ident1 <- "1_Healthy_skin"
HC131$orig.ident2 <- "HC131"
sample <- "Healthy_skin_131"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
HC131[["percent.mt"]] <- PercentageFeatureSet(HC131, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(HC131@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(HC131, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("./202408/HC131/HC131_violin.pdf", width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
HC131 <- subset(HC131, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#HC244
#The $ operator can add columns to object metadata.
HC244$orig.ident1 <- "1_Healthy_skin"
HC244$orig.ident2 <- "HC244"
sample <- "Healthy_skin_HC244"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
HC244[["percent.mt"]] <- PercentageFeatureSet(HC244, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(HC244@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(HC244, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("./202408/HC244/HC244_violin.pdf", width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
HC244 <- subset(HC244, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#UA258
#The $ operator can add columns to object metadata.
UA258$orig.ident1 <- "2_Unaffected_skin"
UA258$orig.ident2 <- "UA258"
sample <- "Unaffected_skin_258"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
UA258[["percent.mt"]] <- PercentageFeatureSet(UA258, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(UA258@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(UA258, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("./202408/UA258/UA258_violin.pdf", width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
UA258 <- subset(UA258, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#AC259
#The $ operator can add columns to object metadata.
AC259$orig.ident1 <- "3_Acute_skin"
AC259$orig.ident2 <- "AC259"
sample <- "Acute_skin_259"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AC259[["percent.mt"]] <- PercentageFeatureSet(AC259, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(AC259@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(AC259, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("./202408/AC259/AC259_violin.pdf", width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
AC259 <- subset(AC259, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

#CH260
#The $ operator can add columns to object metadata.
CH260$orig.ident1 <- "4_Chronic_skin"
CH260$orig.ident2 <- "CH260"
sample <- "Chronic_skin_260"
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
CH260[["percent.mt"]] <- PercentageFeatureSet(CH260, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(CH260@meta.data, 5)
#Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(CH260, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pdf("./202408/CH260/CH260_violin.pdf", width = 15, height = 15)
violin_plot
dev.off()
#Filter out for Feature and mitochondria pct.
CH260 <- subset(CH260, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

##Re-check (after DoubletFinder)
HC131  ##
HC244  ##
UA258  ##
AC259  ##
CH260  ##


##scRNA-seq sample integration
#Set merge file
TN.C <- merge (HC131, y = c(HC244, UA258, AC259, AC260), 
               add.cell.ids = c("HC131","HC244","UA258","AC259","CH260"), project = "groups")
head(TN.C[[]])
table(TN.C$orig.ident)
#Set splitobject file
TN.list<- SplitObject(TN.C, split.by = "orig.ident")
TN.list

# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# NormalizeData,FindVariableFeatures, ScaleData, CellCycleScoring
TN.list <- lapply(X = TN.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = rownames(x), verbose = TRUE)
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
})

#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = TN.list)
#Perform integration 
#reference.list <- TN.list[c("HC131","HC244","UA258","AC259","CH260")]
TN.anchors <- FindIntegrationAnchors(object.list = TN.list, dims = 1:30)
saveRDS(TN.anchors, file = "./202408/TN.anchorsdim30.rds")

#this command creates an 'integrated' data assay
TN.combined <- IntegrateData(anchorset = TN.anchors, dims = 1:30)
#Perform an integrated analysis
#specify that we will perform downstream analysis on the corrected data note that the
#original unmodified data still resides in the 'RNA' assay
DefaultAssay(TN.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
#TN.combined <- ScaleData(TN.combined, verbose = TRUE)
TN.combined <- ScaleData(TN.combined, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(TN.combined), verbose = TRUE)
TN.combined <- RunPCA(TN.combined, npcs = 50, verbose = TRUE)
VizDimLoadings(TN.combined, dims = 1:2)
ElbowPlot(TN.combined, ndims= 50)
# Cluster the cells (1:30)
TN.combined <- FindNeighbors(TN.combined, reduction = "pca", dims = 1:30)
TN.combined <- FindClusters(TN.combined, resolution = 0.4)
table(TN.combined@active.ident)
table(Idents(TN.combined), TN.combined$orig.ident1)

#CellNumber
CellNumber <- table(Idents(TN.combined), TN.combined$orig.ident1)
write.csv(CellNumber , file = "./202408/CellNumber.csv") 

#TN.combined <- RunUMAP(TN.combined, reduction="pca", dims = 1:30)
#DimPlot(TN.combined,reduction = "umap",group.by = "orig.ident")
##UMAP (1:30)
TN.combined <- RunUMAP(TN.combined, reduction = "pca", dims = 1:30)
#output these plots
pdf("./202408/#diffumap_plot.pdf", width = 15, height = 15)
DimPlot(TN.combined, reduction = "umap",label = FALSE, pt.size = 0.8)
DimPlot(TN.combined, reduction = "umap",label = TRUE, pt.size = 0.8)
DimPlot(TN.combined, group.by = "orig.ident", pt.size = 0.8)
DimPlot(TN.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2)
DimPlot(TN.combined, reduction = "umap",label = TRUE, split.by = "TN.merge", pt.size = 0.8, ncol = 3)

##All cells cluster heatmap
Heatmapall <- subset(TN.combined, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34))
Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(Heatmapall.markers)
write.csv(Heatmapall.markers , file = "./202408/Findallmarkers.csv", col.names = TRUE)
top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Heatmapall, features = top10$gene)

##Save the merged Seurat object
saveRDS(TN.combined, file = "./202408/TN.combined_dim30.rds")


#Cellproportion
Cellproportion <- table(Idents(TN.combined), TN.combined$orig.ident1)
Cellproportion <- round(sweep(Cellproportion,MARGIN=2, STATS=colSums(Cellproportion), FUN = "/")*100,2)
write.csv(Cellproportion, file = "C:/Users/user/Desktop/202403/Cellproportion.csv", row.names = T)
Cellproportion <- as.data.frame(Cellproportion)

nb.cols <- nrow(CellNumber)
mycolors <- colorRampPalette(brewer.pal(nb.cols, "Paired"))(nb.cols)

proportion_plot <- ggplot(Cellproportion, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) + 
  geom_col(position = "fill", width = 0.6) + xlab("Sample") + ylab("Proportion") + 
  theme(legend.title = element_blank())+ scale_fill_manual(values = mycolors)

pdf("./202408/proportion_plot.pdf", width = 15, height = 15)
proportion_plot
dev.off()

###RUN CELL ANNOTATION TOOLS HERE and do annotation manually


TN.combined <- readRDS(file = "./202408/TN.combined_dim30.rds")
##Rename cluster (for example_Fibroblast, Macrophage......)
new.cluster.ids <- c("C0_FB","C1_KC","C2_Mac","C3_vEC","C4_KC","C5_Tcell", 
                     "C6_FB","C7_FB","C8_FB","C9_KC","C10_KC","C11_SMC",
                     "C12_KC","C13_KC","C14_Mo","C15_MC","C16_MLA","C17_lEC",
                     "C18_Mo","C19_KC","C20_KC","C21_DC")
names(new.cluster.ids) <- levels(TN.combined)
TNname.combined <- RenameIdents(TN.combined, new.cluster.ids)
#DimPlot(TNname.combined, reduction = "umap",label = FALSE, pt.size = 0.8, cols = mycolors)
#DimPlot(TNname.combined, reduction = "umap",label = TRUE, pt.size = 0.8, repel = T, cols = mycolors)
#DimPlot(TNname.combined, group.by = "orig.ident1", pt.size = 0.8, cols = c("#00BA38","#00B9E3","#F8766D","#DB72FB"))
#DimPlot(TNname.combined, reduction = "umap",label = FALSE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, cols = mycolors)
UMAP_cluster_plot_all<- DimPlot(TNname.combined, reduction = "umap",label = TRUE, pt.size = 0.8, label.size=7, repel = T, cols = mycolors)
UMAP_plot_subset<- DimPlot(TNname.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, repel = T, cols = mycolors)

pdf("C:/Users/user/Desktop/202403/UMAP_cluster_plot_all.pdf", width = 15, height = 15)
UMAP_cluster_plot_all
dev.off()

pdf("C:/Users/user/Desktop/202403/UMAP_cluster_plot_subset.pdf", width = 15, height = 15)
UMAP_plot_subset
dev.off()

#Cell type propotion 

# rename identify from cluster to cell type
celltype.ids <- c("FB","KC","Mac","vEC","KC","Tcell", 
                  "FB","FB","FB","KC","KC","SMC",
                  "KC","KC","Mo","MC","MLA","lEC",
                  "Mo","KC","KC","DC")

names(celltype.ids) <- levels(TN.combined)
TNtype.combined <- RenameIdents(TN.combined, celltype.ids)

nb.cols <- length(unique(celltype.ids))
typecolors <- colorRampPalette(brewer.pal(nb.cols, "Paired"))(nb.cols)

# plot UMAP 
UMAP_type_plot_all<- DimPlot(TNtype.combined, reduction = "umap",label = TRUE, pt.size = 0.8, label.size=7, repel = T, cols = typecolors)
UMAP_type_plot_subset<- DimPlot(TNtype.combined, reduction = "umap",label = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, repel = T, cols = typecolors)

pdf("C:/Users/user/Desktop/202403/UMAP_celltype_plot_all.pdf", width = 15, height = 15)
UMAP_type_plot_all
dev.off()

pdf("C:/Users/user/Desktop/202403/UMAP_celltype_plot_subset.pdf", width = 15, height = 15)
UMAP_type_plot_subset
dev.off()

#Cell type proportion
celltypeproportion <- table(Idents(TNtype.combined), TN.combined$orig.ident1)
celltypeproportion <- round(sweep(celltypeproportion,MARGIN=2, STATS=colSums(celltypeproportion), FUN = "/")*100,2)
write.csv(celltypeproportion, file = "C:/Users/user/Desktop/202403/Cell_type_proportion.csv", row.names = T)
celltypeproportion <- as.data.frame(celltypeproportion)
type_proportion_plot <- ggplot(celltypeproportion, aes(x = Var2, y = Freq, fill = Var1)) + theme_bw(base_size = 15) + 
  geom_col(position = "fill", width = 0.6) + xlab("Sample") + ylab("Proportion") + labs(title="Cell type propotion") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size=24))+ scale_fill_manual(values = typecolors)

pdf("C:/Users/user/Desktop/202403/proportion_plot.pdf", width = 15, height = 15)
type_proportion_plot
dev.off()



##All cells cluster heatmap
Heatmapall <- subset(TN.combined, idents = seq(0,21)) 
Heatmapall.markers <- FindAllMarkers(Heatmapall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
head(Heatmapall.markers)
write.csv(Heatmapall.markers , file = "C:/Users/user/Desktop/202403/Findallmarkers_2.csv", col.names = TRUE) 
top10 <- Heatmapall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
marker_heatmap <- DoHeatmap(Heatmapall, features = top10$gene)

pdf("C:/Users/user/Desktop/202403/marker_heatmap.pdf", width = 15, height = 15)
marker_heatmap
dev.off()

##Save the merged Seurat object
saveRDS(TN.combined, file = "C:/Users/user/Desktop/202403/TN.combined_20240315.rds")


