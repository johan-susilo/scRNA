library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tidyr)
library(dplyr)
library(ggrepel)

out_dir <- "/home/johan/output/skin_pmh/dge/macrophage"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

pmh_obj <- readRDS("/home/johan/output/skin_pmh/TN.combined_dim30.rds")

pmh_obj$condition <- ifelse(
  grepl("HTY|UA", pmh_obj$orig.ident2, ignore.case = TRUE),
  "Healthy",
  "PMH"
)

# verify
table(pmh_obj$orig.ident2, pmh_obj$condition)


# === STEP 1: Fresh AggregateExpression ===
pb_list <- AggregateExpression(
  pmh_obj, 
  assays = "RNA", 
  slot = "counts",
  group.by = c("seurat_clusters", "orig.ident", "condition"),
  return.seurat = TRUE
)

print("pb_list contents:")
print(names(pb_list))
print("pb_list length:")
print(length(pb_list))

# === STEP 2: Extract RNA assay ===
mac_pb_full <- pb_list[["RNA"]]
print(colnames(mac_pb_full)[1:10]) #g0 represent cluster 0 aggregate

# === STEP 3: Select cluster 6 ===
cluster6_cols <- grep("^g6_", colnames(mac_pb_full), value = TRUE)

# 2. Subset WITH drop = FALSE to preserve colnames even if only 1 match
mac_pb_subset <- mac_pb_full[, cluster6_cols, drop = FALSE]

# 4. Create metadata
mac_meta <- data.frame(pseudobulk_id = cluster6_cols) %>%
  # remove = FALSE keeps the original column so we can use it for rownames later
  tidyr::separate(pseudobulk_id, into = c("cluster", "sample", "condition"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(sample_id = sample) %>%
  dplyr::select(pseudobulk_id, sample_id, condition)   # Keep pseudobulk_id in the select statement

# Set rownames
rownames(mac_meta) <- mac_meta$pseudobulk_id
print(head(mac_meta)) # Verify it worked

#extract count matrix for DESeq2
counts_matrix <- as.matrix(GetAssayData(mac_pb_full, slot = "counts"))
counts_cluster6 <- counts_matrix[, rownames(mac_meta), drop = FALSE]

# --- Verify the order matches perfectly ---
counts_cluster6 <- counts_cluster6[, rownames(mac_meta)]
all(colnames(counts_cluster6) == rownames(mac_meta)) # This should print TRUE

# --- Create the DESeq2 object ---
dds <- DESeqDataSetFromMatrix(
  countData = counts_cluster6,
  colData = mac_meta,
  design = ~ condition
)

# run deseq pipeline
dds <- DESeq(dds)
#extract result
res <- results(dds, contrast = c("condition", "PMH", "Healthy"))

res_ordered <- res[order(res$padj), ]
print(head(res_ordered))

#create a data frame for plotting
res_df <- as.data.frame(res)
res_df <- res_df %>% filter(!is.na(padj)) #filter out na
res_df$gene <- rownames(res_df) #move row names (gene) into column

res_df <- res_df %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  select(gene, everything())

write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)

# 1. Create a subset of the top 10 UP and top 10 DOWN genes to label
top_genes <- res_df %>%
  filter(significance != "Not Significant") %>%  # Keep only the significant ones
  group_by(significance) %>%                     # Split them into the "Up" and "Down" groups
  arrange(padj) %>%                              # Sort by the most significant (lowest padj)
  slice_head(n = 10) %>%                         # Grab the top 10 rows from EACH group
  ungroup()                                      # Ungroup to return to a normal dataframe

#ma plot
# Open a PDF file
pdf("/home/johan/output/skin_pmh/dge/macrophage/MA_plot_cluster6.pdf", width = 6, height = 5)

# Generate the MA plot directly from the DESeq2 results object (res)
# ylim shrinks the y-axis to keep extreme outliers from zooming the plot out too far
DESeq2::plotMA(res, main = "MA Plot: PMH vs Healthy", ylim = c(-8, 8))

# Close and save the file
dev.off()

#volcano plot

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_genes, 
                  aes(label = gene), 
                  color = "black",   # Makes the text easy to read
                  box.padding = 0.5, # Adds a little physical space around the text
                  max.overlaps = Inf) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Volcano Plot: PMH vs Healthy")

ggsave(filename = file.path(out_dir, "volcano_plot_cluster6.pdf"),
       plot = volcano_plot,
       width = 6,
       height = 5)

