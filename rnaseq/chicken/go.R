#!/usr/bin/env Rscript

.libPaths("/home/johan/R/x86_64-pc-linux-gnu-library/4.2")


suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(org.Gg.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  
})

# Inputs
input_file <- "/home/johan/johan/output/chicken/DGE/DEseq2_normalized.csv"
gaf_file   <- "/home/johan/johan/reference/chicken/go/goa_chicken.gaf"
output_dir <- "/home/johan/johan/output/chicken/GO"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Load gene list (expect Ensembl gene IDs in rownames)
gene_df <- read.csv(input_file, row.names = 1, check.names = FALSE)
genes_ensem <- rownames(gene_df)
if (length(genes_ensem) == 0) stop("Gene table has no rows -> cannot perform GO enrichment.")

# 2) Prepare ID mappings for two strategies
# 2a) Ensembl -> ENTREZ (for enrichGO with OrgDb)
map_entrez <- AnnotationDbi::mapIds(
  org.Gg.eg.db,
  keys = genes_ensem,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
map_entrez <- map_entrez[!is.na(map_entrez)]
genes_entrez <- unique(unname(map_entrez))

# 2b) Ensembl -> UniProt (to match GOA GAF accessions if we use GAF)
suppressWarnings({
  sel <- AnnotationDbi::select(org.Gg.eg.db, keys = genes_ensem,
                               keytype = "ENSEMBL", columns = c("UNIPROT"))
})
sel <- sel[!is.na(sel$UNIPROT) & nzchar(sel$UNIPROT), , drop = FALSE]
genes_uniprot <- unique(sel$UNIPROT)

cat("Mapped to ENTREZ:", length(genes_entrez), "| to UniProt:", length(genes_uniprot),
    "| from Ensembl:", length(genes_ensem), "\n")

# Helper to safely save plots
save_plot <- function(plot_obj, filepath, width = 9, height = 7) {
  tryCatch({
    ggsave(filename = filepath, plot = plot_obj, width = width, height = height, limitsize = FALSE)
    cat("Saved:", filepath, "\n")
  }, error = function(e) {
    cat("Failed to save", basename(filepath), ":", conditionMessage(e), "\n")
  })
}

# Build plots from an enrichResult
generate_plots <- function(eresult, label) {
  if (is.null(eresult) || nrow(eresult) == 0) {
    cat("No enriched GO terms found for:", label, "\n")
    return(invisible(NULL))
  }

  cat("Top terms for", label, ":\n"); print(utils::head(eresult))

  # Compute similarity matrix where possible
  es <- tryCatch({ pairwise_termsim(eresult) }, error = function(e) eresult)

  # Build plots one by one to be robust against missing deps
  # Dotplot
  p <- tryCatch({ dotplot(eresult, showCategory = min(20, nrow(eresult))) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(output_dir, paste0("GO_", label, "_dotplot.pdf")))

  # Barplot
  p <- tryCatch({ barplot(eresult, showCategory = min(20, nrow(eresult))) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(output_dir, paste0("GO_", label, "_barplot.pdf")))

  # Cnetplot
  p <- tryCatch({ cnetplot(es, showCategory = min(10, nrow(eresult)), circular = FALSE, colorEdge = TRUE) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(output_dir, paste0("GO_", label, "_cnetplot.pdf")))

  # Emapplot
  p <- tryCatch({ emapplot(es, layout = "kk") }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(output_dir, paste0("GO_", label, "_emapplot.pdf")))

  # Treeplot
  p <- tryCatch({ treeplot(eresult) }, error = function(e) NULL)
  if (!is.null(p)) save_plot(p, file.path(output_dir, paste0("GO_", label, "_treeplot.pdf")))
}

# Strategy A: Use GAF if available (build custom TERM2GENE and use enricher)
use_gaf <- file.exists(gaf_file)
if (use_gaf && length(genes_uniprot) > 0) {
  cat("\nUsing GAF:", gaf_file, "\n")
  gaf <- tryCatch(
    read.delim(gaf_file, header = FALSE, comment.char = "!", quote = "", stringsAsFactors = FALSE, fill = TRUE),
    error = function(e) NULL
  )
  if (!is.null(gaf) && ncol(gaf) >= 9) {
    # Columns (GAF 2.x): V2 = DB Object ID (UniProt), V5 = GO ID, V9 = Aspect (P/F/C)
    gaf <- gaf[!is.na(gaf$V2) & nzchar(gaf$V2) & !is.na(gaf$V5) & nzchar(gaf$V5), ]
    term2gene <- unique(gaf[, c(5, 2)])
    colnames(term2gene) <- c("GO", "GENE")

    # Optional: GO term names if GO.db installed
    term2name <- NULL
    if (requireNamespace("GO.db", quietly = TRUE)) {
      suppressPackageStartupMessages(library(GO.db))
      term2name <- AnnotationDbi::select(GO.db, keys = unique(term2gene$GO), columns = "TERM", keytype = "GOID")
      colnames(term2name) <- c("GO", "NAME")
      term2name <- term2name[!is.na(term2name$NAME), , drop = FALSE]
    }

    # Per-ontology mapping using Aspect (P/BP, F/MF, C/CC)
    gaf$V9[is.na(gaf$V9)] <- ""
    onto_map <- list(
      BP = unique(gaf[gaf$V9 == "P", c(5, 2)]),
      MF = unique(gaf[gaf$V9 == "F", c(5, 2)]),
      CC = unique(gaf[gaf$V9 == "C", c(5, 2)])
    )
    for (nm in names(onto_map)) {
      df <- onto_map[[nm]]
      colnames(df) <- c("GO", "GENE")
      if (nrow(df) == 0) {
        cat("No mapping rows in GAF for", nm, "\n"); next
      }
      cat("\nRunning GAF-based enrichment for", nm, "(genes:", length(genes_uniprot), ")\n")
      eres <- tryCatch({
        enricher(gene = genes_uniprot, TERM2GENE = df, TERM2NAME = term2name, pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      }, error = function(e) NULL)
      generate_plots(eres, paste0(nm, "_GAF"))
    }
  } else {
    cat("Failed to read/parse GAF. Falling back to OrgDb enrichGO.\n")
    use_gaf <- FALSE
  }
} else if (use_gaf && length(genes_uniprot) == 0) {
  cat("GAF found but no UniProt IDs mapped from your Ensembl list. Falling back to OrgDb enrichGO.\n")
  use_gaf <- FALSE
}

# Strategy B: Fallback to standard enrichGO with OrgDb
if (!use_gaf) {
  if (length(genes_entrez) == 0) stop("No genes mapped to ENTREZ -> cannot perform GO enrichment.")
  for (ont in c("BP", "CC", "MF")) {
    cat("\nRunning enrichGO for", ont, "(ENTREZ IDs:", length(genes_entrez), ")...\n")
    ego <- enrichGO(
      gene = genes_entrez,
      OrgDb = org.Gg.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    generate_plots(ego, ont)
  }
}

cat("\nGO enrichment workflow complete. Plots saved to:\n", output_dir, "\n")
