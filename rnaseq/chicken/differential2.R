#!/usr/bin/env Rscript

# Simple DESeq2 differential expression analysis
# Usage: Rscript differential.R --counts counts.tsv --metadata metadata.csv --outdir results/ [options...]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)  # Add tibble for column_to_rownames
  library(rlang)   # Add this for %||% operator
  library(ggrepel) # <- for label styling
})

# Small helper (used by the added style block)
ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- gsub("^--", "", args[i])
      if (i < length(args) && !startsWith(args[i+1], "--")) {
        params[[key]] <- args[i+1]
        i <- i + 2
      } else {
        params[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

# Set defaults - fix the %||% usage
counts_file <- if(is.null(params$counts)) stop("--counts required") else params$counts
metadata_file <- if(is.null(params$metadata)) stop("--metadata required") else params$metadata
outdir <- if(is.null(params$outdir)) "." else params$outdir
sample_col <- if(is.null(params$sampleColumn)) "sample" else params$sampleColumn
contrast_factor <- if(is.null(params$contrastFactor)) "egg_production" else params$contrastFactor
numerator <- if(is.null(params$contrastNumerator)) "high_egg" else params$contrastNumerator
denominator <- if(is.null(params$contrastDenominator)) "low_egg" else params$contrastDenominator
batch_col <- params$batchColumn
padj_thresh <- if(is.null(params$padj)) 0.05 else as.numeric(params$padj)
lfc_thresh <- if(is.null(params$lfc)) 1 else as.numeric(params$lfc)

cat("Loading data...\n")
# Load counts and metadata
counts <- read_tsv(counts_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

# Prepare count matrix
count_matrix <- counts %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

# Ensure metadata matches count columns
metadata <- metadata %>%
  filter(.data[[sample_col]] %in% colnames(count_matrix)) %>%
  arrange(match(.data[[sample_col]], colnames(count_matrix)))

count_matrix <- count_matrix[, metadata[[sample_col]]]

cat(sprintf("Counts: %d genes x %d samples\n", nrow(count_matrix), ncol(count_matrix)))
cat(sprintf("Metadata: %d samples\n", nrow(metadata)))

# Create DESeq2 object
cat("Creating DESeq2 object...\n")
if (!is.null(batch_col) && batch_col %in% colnames(metadata)) {
  design_formula <- as.formula(paste("~", batch_col, "+", contrast_factor))
  cat(sprintf("Design: ~ %s + %s\n", batch_col, contrast_factor))
} else {
  design_formula <- as.formula(paste("~", contrast_factor))
  cat(sprintf("Design: ~ %s\n", contrast_factor))
}

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = design_formula
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("After filtering: %d genes\n", nrow(dds)))

# Run DESeq2
cat("Running DESeq2...\n")
dds <- DESeq(dds)

# Get results
res <- results(dds, 
               contrast = c(contrast_factor, numerator, denominator),
               alpha = padj_thresh)

cat(sprintf("Contrast: %s %s vs %s\n", contrast_factor, numerator, denominator))
cat(sprintf("Significant genes (padj < %g, |LFC| > %g): %d\n", 
            padj_thresh, lfc_thresh, 
            sum(res$padj < padj_thresh & abs(res$log2FoldChange) > lfc_thresh, na.rm = TRUE)))

# Save results
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
results_file <- file.path(outdir, "deseq2_results.csv")
res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  arrange(padj, desc(abs(log2FoldChange)))

write_csv(res_df, results_file)
cat(sprintf("Results saved: %s\n", results_file))

# PCA plot (your original)
cat("Generating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c(contrast_factor), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = .data[[contrast_factor]])) +
  geom_point(size = 3) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    title = "PCA Plot",
    color = contrast_factor
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(outdir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# ---------------------------- ADDED: styled PCA + scree + exports ----------------------------
cat("Generating styled PCA + scree (additional outputs)…\n")

# Build matrix for PCA from your VST
mat_for_pca <- assay(vsd)

# Select top variable genes (safe default)
top_var <- min(5000L, nrow(mat_for_pca))
if (top_var > 1L && top_var < nrow(mat_for_pca)) {
  v <- apply(mat_for_pca, 1, var, na.rm = TRUE)
  sel <- order(v, decreasing = TRUE)[seq_len(top_var)]
  mat_pca <- mat_for_pca[sel, , drop = FALSE]
} else {
  mat_pca <- mat_for_pca
}

# Run PCA on samples
pca <- prcomp(t(mat_pca), center = TRUE, scale. = FALSE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

# Bind metadata for plotting (pc_df)
pc_meta <- as.data.frame(metadata)
rownames(pc_meta) <- pc_meta[[sample_col]]
pc_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  pc_meta[rownames(pca$x), , drop = FALSE],
  check.names = FALSE
)

# Styling controls (use your existing columns)
color_by   <- contrast_factor
shape_by   <- NULL                  # set e.g. "tissue" if you want shapes
label_by   <- sample_col
ellipse_by <- color_by
title_txt  <- "PCA (DESeq2 VST)"

# ---------------- plotting (unchanged visualization) ----------------
xlab_txt <- sprintf("PC1 (%.1f%%)", 100*var_expl[1])
ylab_txt <- sprintf("PC2 (%.1f%%)", 100*var_expl[2])

aes_base <- aes(x = PC1, y = PC2)
if (!is.null(color_by) && color_by %in% colnames(pc_df)) aes_base <- modifyList(aes_base, aes(color = .data[[color_by]]))
if (!is.null(shape_by) && shape_by %in% colnames(pc_df)) aes_base <- modifyList(aes_base, aes(shape = .data[[shape_by]]))

label_layer <- NULL
if (!is.null(label_by) && label_by %in% colnames(pc_df)) {
  label_layer <- geom_text_repel(aes(label = .data[[label_by]]), max.overlaps = 100, size = 3)
}
ellipse_layer <- NULL
if (!is.null(ellipse_by) && ellipse_by %in% colnames(pc_df)) {
  gdf <- pc_df[!is.na(pc_df[[ellipse_by]]), , drop = FALSE]
  ok_groups <- names(which(table(gdf[[ellipse_by]]) >= 3))
  if (length(ok_groups) > 0) {
    ellipse_layer <- stat_ellipse(aes(group = .data[[ellipse_by]]), type = "norm", level = 0.95,
                                  linewidth = 0.5, linetype = 2, alpha = 0.2)
  }
}

p <- ggplot(pc_df, aes_base) +
  geom_point(size = 3, alpha = 0.9) +
  label_layer +
  ellipse_layer +
  theme_bw(base_size = 12) +
  labs(title = title_txt, x = xlab_txt, y = ylab_txt,
       color = if (!is.null(color_by) && color_by %in% colnames(pc_df)) color_by else NULL,
       shape = if (!is.null(shape_by) && shape_by %in% colnames(pc_df)) shape_by else NULL) +
  theme(plot.title = element_text(face = "bold"))

# Scree (kept)
scree_df <- data.frame(PC = paste0("PC", seq_along(var_expl)), Var = 100*var_expl)
pscree <- ggplot(scree_df[1:min(10, nrow(scree_df)), ], aes(x = PC, y = Var)) +
  geom_col() + geom_text(aes(label = sprintf("%.1f%%", Var)), vjust = -0.4, size = 3) +
  labs(title = "Scree plot (variance explained)", x = NULL, y = "% variance") +
  theme_bw(base_size = 12)

# ---------------- outputs (unchanged filenames adapted to outdir) ----------------
fig_dir <- file.path(outdir, "figs"); ensure_dir(fig_dir)

# Matrices (based on raw counts in dds)
libsizes <- colSums(counts(dds))
cpm      <- sweep(counts(dds), 2, libsizes, "/") * 1e6
logcpm   <- log2(cpm + 1)

cpm_fp      <- file.path(outdir, "cpm_counts.tsv")
logcpm_fp   <- file.path(outdir, "logcpm_counts.tsv")
pcacoord_fp <- file.path(outdir, "pca_coordinates.tsv")
pdf_fp      <- file.path(fig_dir, "pca_samples.pdf")
png_fp      <- file.path(fig_dir, "pca_samples.png")
scree_fp    <- file.path(fig_dir, "pca_scree.png")
vst_fp      <- file.path(outdir, "vst_counts.tsv")

write.table(round(cpm,    6), file = cpm_fp,    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(round(logcpm, 6), file = logcpm_fp, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(round(mat_for_pca, 6), file = vst_fp, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# PCA coords + metadata
pc_out <- pc_df; rownames(pc_out) <- NULL
write.table(pc_out, file = pcacoord_fp, sep = "\t", quote = FALSE, row.names = FALSE)

# figures
ggsave(pdf_fp,   plot = p,      width = 8, height = 6, units = "in")
ggsave(png_fp,   plot = p,      width = 8, height = 6, units = "in", dpi = 300)
ggsave(scree_fp, plot = pscree, width = 7, height = 4, units = "in", dpi = 300)

cat(sprintf("Wrote: %s\n", cpm_fp))
cat(sprintf("Wrote: %s\n", logcpm_fp))
cat(sprintf("Wrote: %s\n", pcacoord_fp))
cat(sprintf("Wrote: %s\n", pdf_fp))
cat(sprintf("Wrote: %s\n", png_fp))
cat(sprintf("Wrote: %s\n", scree_fp))

# -------------------------------------------------------------------------------------------

# Volcano plot
cat("Generating volcano plot...\n")
volcano_data <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    significant = padj < padj_thresh & abs(log2FoldChange) > lfc_thresh,
    direction = case_when(
      log2FoldChange > lfc_thresh & padj < padj_thresh ~ "Up",
      log2FoldChange < -lfc_thresh & padj < padj_thresh ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_plot <- ggplot(volcano_data, aes(log2FoldChange, -log10(padj), color = direction)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", alpha = 0.5) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    title = sprintf("Volcano Plot: %s vs %s", numerator, denominator),
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(outdir, "volcano_plot.png"), volcano_plot, width = 8, height = 6, dpi = 300)

# Summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Total genes analyzed: %d\n", nrow(res_df)))
cat(sprintf("Upregulated (LFC > %g, padj < %g): %d\n", 
            lfc_thresh, padj_thresh,
            sum(volcano_data$direction == "Up")))
cat(sprintf("Downregulated (LFC < -%g, padj < %g): %d\n", 
            lfc_thresh, padj_thresh,
            sum(volcano_data$direction == "Down")))
cat(sprintf("Results saved to: %s\n", outdir))
