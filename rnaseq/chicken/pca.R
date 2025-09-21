#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  ok <- c(
    requireNamespace("ggplot2", quietly = TRUE),
    requireNamespace("ggrepel", quietly = TRUE)
  )
  if (!all(ok)) {
    stop("Please install: install.packages(c('ggplot2','ggrepel'))")
  }
})

library(ggplot2)
library(ggrepel)

# --------- arg parsing ---------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

# I/O and options (defaults compatible with your pipeline)
output_base <- get_arg("output_base", "/home/johan/johan/output/chicken")
metadata_fp <- get_arg("metadata",    "/home/johan/pipeline/rnaseq/chicken/sample.csv")
combined_fp <- get_arg("combined",    NA)    # e.g. /mnt/data/combined.HTseq_report
pattern     <- get_arg("pattern",     ".htseq.no.counts.tsv")
min_count   <- as.integer(get_arg("min_count", "10"))

# PCA/transform options
pseudocount <- as.numeric(get_arg("pseudocount", "1"))  # for log2(CPM + p)
top_var     <- as.integer(get_arg("top_var", "5000"))   # keep top-N variable genes (after logCPM); set 0 to keep all
center_genes <- tolower(get_arg("center", "true")) %in% c("true","t","1","yes","y")
scale_genes  <- tolower(get_arg("scale",  "false")) %in% c("true","t","1","yes","y")

# Plot aesthetics (metadata columns)
color_by   <- get_arg("color_by",   "tissue")
shape_by   <- get_arg("shape_by",   "egg_production")
label_by   <- get_arg("label_by",   "sample")
ellipse_by <- get_arg("ellipse_by", color_by)

# --------- helpers ---------
msg <- function(...) message(sprintf(...))
ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

read_combined <- function(fp) {
  df <- read.table(fp, header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  first_col <- colnames(df)[1]
  # be flexible on first column name
  if (!grepl("gene|id", tolower(first_col))) {
    stop("Combined file must have first column as gene id (e.g., GeneID). Found: ", first_col)
  }
  # guard against duplicated 'GeneID' columns
  if (any(duplicated(colnames(df)))) {
    dup <- unique(colnames(df)[duplicated(colnames(df))])
    stop("Duplicate column names in combined matrix (e.g., ", paste(dup, collapse=", "),
         "). Please fix headers.")
  }
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  # coerce to integer-like numeric
  for (j in seq_len(ncol(df))) df[[j]] <- suppressWarnings(as.numeric(as.character(df[[j]])))
  df[is.na(df)] <- 0
  as.data.frame(df, check.names = FALSE)
}

read_counts_from_output <- function(base_dir, patt) {
  files <- Sys.glob(file.path(base_dir, "*", "count", paste0("*", patt)))
  if (length(files) == 0) stop("No htseq count files found with pattern '", patt, "' under ", base_dir)
  msg("Found %d count files", length(files))
  mats <- list()
  for (f in files) {
    sm <- sub("\\.htseq.*$", "", basename(f))
    dat <- tryCatch(
      read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("GeneID", sm), quote = "", comment.char = ""),
      error = function(e) stop("Failed reading ", f, ": ", e$message)
    )
    dat <- dat[!grepl("^__", dat$GeneID), , drop = FALSE]
    mats[[sm]] <- dat
  }
  merged <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), mats)
  rownames(merged) <- merged$GeneID
  merged$GeneID <- NULL
  merged[is.na(merged)] <- 0
  for (j in seq_len(ncol(merged))) merged[[j]] <- suppressWarnings(as.numeric(as.character(merged[[j]])))
  as.data.frame(merged, check.names = FALSE)
}

# --------- load data ---------
stopifnot(dir.exists(output_base))
if (!file.exists(metadata_fp)) stop("metadata file not found: ", metadata_fp)

counts <- if (!is.na(combined_fp) && file.exists(combined_fp)) {
  msg("Reading combined matrix: %s", combined_fp)
  read_combined(combined_fp)
} else {
  msg("Merging per-sample HTSeq result files from: %s", output_base)
  read_counts_from_output(output_base, pattern)
}

if (ncol(counts) < 2) stop("At least 2 samples are required for PCA; found: ", ncol(counts))

# Basic gene filter by total counts
keep <- rowSums(counts) >= min_count
counts <- counts[keep, , drop = FALSE]
msg("Kept %d genes after filtering (rowSums >= %d)", nrow(counts), min_count)

# metadata
meta <- read.csv(metadata_fp, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(meta)) stop("metadata must contain a 'sample' column")

# align to counts (order metadata to match columns)
meta <- meta[match(colnames(counts), meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample

# fill missing rows minimally
if (any(is.na(meta$sample))) {
  miss <- colnames(counts)[is.na(meta$sample)]
  warning("Samples missing in metadata: ", paste(miss, collapse = ", "),
          ". Creating minimal metadata rows.")
  add <- data.frame(sample = miss, stringsAsFactors = FALSE)
  meta[is.na(meta$sample), "sample"] <- miss[is.na(meta$sample)]
  rownames(meta)[is.na(rownames(meta))] <- miss[is.na(rownames(meta))]
}

# ensure final order
meta <- meta[colnames(counts), , drop = FALSE]

# --------- PURE normalization & transform (no DESeq2) ---------
# CPM
libsizes <- colSums(counts)
if (any(libsizes == 0)) {
  bad <- names(libsizes)[libsizes == 0]
  stop("Zero library size for samples: ", paste(bad, collapse=", "))
}
cpm <- sweep(counts, 2, libsizes, "/") * 1e6

# log2(CPM + pseudocount)
logcpm <- log2(cpm + pseudocount)

# Optionally keep top-N variable genes to reduce noise/high-dim effects
if (!is.na(top_var) && top_var > 0 && top_var < nrow(logcpm)) {
  v <- apply(logcpm, 1, var, na.rm = TRUE)
  sel <- order(v, decreasing = TRUE)[seq_len(top_var)]
  logcpm_pca <- logcpm[sel, , drop = FALSE]
  msg("Selected top %d variable genes for PCA", length(sel))
} else {
  logcpm_pca <- logcpm
}

# --------- PCA ---------
pca <- prcomp(t(logcpm_pca), center = center_genes, scale. = scale_genes)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

# Avoid duplicate 'sample' column when binding metadata
pc_meta <- meta[rownames(pca$x), , drop = FALSE]
if ("sample" %in% colnames(pc_meta)) pc_meta$sample <- NULL
pc_df <- data.frame(sample = rownames(pca$x),
                    PC1 = pca$x[, 1],
                    PC2 = pca$x[, 2],
                    pc_meta,
                    check.names = FALSE)

# --------- plotting ---------
title_txt <- sprintf("PCA (log2(CPM + %g))%s%s",
                     pseudocount,
                     if (center_genes) ", centered" else "",
                     if (scale_genes) ", scaled" else "")
xlab_txt  <- sprintf("PC1 (%.1f%%)", 100*var_expl[1])
ylab_txt  <- sprintf("PC2 (%.1f%%)", 100*var_expl[2])

aes_base <- aes(x = PC1, y = PC2)
if (!is.null(color_by)  && color_by  %in% colnames(pc_df)) aes_base <- modifyList(aes_base, aes(color = .data[[color_by]]))
if (!is.null(shape_by)  && shape_by  %in% colnames(pc_df)) aes_base <- modifyList(aes_base, aes(shape = .data[[shape_by]]))

# Build optional layers first to avoid inline if/else parse issues
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

# Scree plot
scree_df <- data.frame(PC = paste0("PC", seq_along(var_expl)),
                       Var = 100*var_expl)
pscree <- ggplot(scree_df[1:min(10, nrow(scree_df)), ], aes(x = PC, y = Var)) +
  geom_col() + geom_text(aes(label = sprintf("%.1f%%", Var)), vjust = -0.4, size = 3) +
  labs(title = "Scree plot (variance explained)", x = NULL, y = "% variance") +
  theme_bw(base_size = 12)

# --------- outputs ---------
fig_dir <- file.path(output_base, "figs")
out_dir <- file.path(output_base)
ensure_dir(fig_dir)

# matrices
cpm_fp     <- file.path(out_dir, "cpm_counts.tsv")
logcpm_fp  <- file.path(out_dir, "logcpm_counts.tsv")
pcacoord_fp<- file.path(out_dir, "pca_coordinates.tsv")

# figures
pdf_fp     <- file.path(fig_dir, "pca_samples.pdf")
png_fp     <- file.path(fig_dir, "pca_samples.png")
scree_fp   <- file.path(fig_dir, "pca_scree.png")

# save matrices
write.table(round(cpm, 6),   file = cpm_fp,    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(round(logcpm, 6),file = logcpm_fp, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# PCA coords + metadata
pc_out <- pc_df
rownames(pc_out) <- NULL
write.table(pc_out, file = pcacoord_fp, sep = "\t", quote = FALSE, row.names = FALSE)

# save plots
ggsave(pdf_fp, plot = p, width = 8, height = 6, units = "in")
ggsave(png_fp, plot = p, width = 8, height = 6, units = "in", dpi = 300)
ggsave(scree_fp, plot = pscree, width = 7, height = 4, units = "in", dpi = 300)

msg("Wrote: %s", cpm_fp)
msg("Wrote: %s", logcpm_fp)
msg("Wrote: %s", pcacoord_fp)
msg("Wrote: %s", pdf_fp)
msg("Wrote: %s", png_fp)
msg("Wrote: %s", scree_fp)
