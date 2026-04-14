suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db) # Human database for gene translation
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# ==============================================================================
# 1. Define the clusterProfiler Analysis Function
# ==============================================================================
run_pathway_analysis <- function(csv_path, out_base_dir) {
  
  # Extract names for titling and saving
  file_name <- basename(csv_path)
  clean_name <- str_remove(file_name, "^deseq2_") %>% str_remove("\\.csv$")
  cell_type <- str_extract(clean_name, "^[^_]+")
  
  message(paste("\n-> Running clusterProfiler Analysis for:", clean_name))
  
  out_dir <- file.path(out_base_dir, cell_type, "pathways")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  
  # Load DESeq2 DGE data
  dge_data <- read.csv(csv_path)
  
  # --- STEP 1: Gene ID Translation (Symbol to Entrez) ---
  # clusterProfiler requires Entrez IDs. We must map your gene symbols.
  dge_data <- dge_data %>% filter(!is.na(padj))
  
  mapped_genes <- bitr(dge_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Merge the Entrez IDs back into your DGE results
  dge_data <- left_join(dge_data, mapped_genes, by = c("gene" = "SYMBOL")) %>% 
    filter(!is.na(ENTREZID))
  
  # Define the custom background (Universe)
  universe_entrez <- dge_data$ENTREZID
  
  # --- STEP 2: Extract Significant Upregulated and Downregulated Genes ---
  genes_up <- dge_data %>%
    filter(log2FoldChange > 1 & padj < 0.05) %>%
    pull(ENTREZID)
  
  genes_down <- dge_data %>%
    filter(log2FoldChange < -1 & padj < 0.05) %>%
    pull(ENTREZID)
  
  # Create a named list for clusterProfiler's compareCluster
  query_list <- list()
  if(length(genes_up) > 0) query_list[["Upregulated_in_PMH"]] <- genes_up
  if(length(genes_down) > 0) query_list[["Downregulated_in_PMH"]] <- genes_down
  
  if (length(query_list) == 0) {
    message("   - SKIPPING: No significant DE genes found to analyze.")
    return(NULL)
  }
  
  # ============================================================================
  # 3A. Run GO: Biological Process (GO:BP)
  # ============================================================================
  message("   - Running GO:BP...")
  go_res <- compareCluster(
    geneCluster = query_list,
    fun = "enrichGO",
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE # Automatically translates Entrez back to human-readable Symbols for the plot!
  )
  
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    
    # --- FIXED CSV SAVING: Replace slashes with commas ---
    go_df <- as.data.frame(go_res)
    go_df$geneID <- str_replace_all(go_df$geneID, "/", ", ")
    write.csv(go_df, file.path(out_dir, paste0(clean_name, "_GOBP_results.csv")), row.names = FALSE)
    # -----------------------------------------------------
    
    # Generate native clusterProfiler DotPlot
    p_go <- dotplot(go_res, showCategory = 10) +
      ggtitle(paste("Biological Processes (GO:BP) -", tools::toTitleCase(str_replace_all(clean_name, "_", " ")))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(filename = file.path(out_dir, paste0(clean_name, "_GOBP_bubbleplot.png")), plot = p_go, width = 10, height = 7, dpi = 300)
  } else {
    message("   - No significant GO:BP terms found.")
  }

  # ============================================================================
  # 3B. Run KEGG Pathways
  # ============================================================================
  message("   - Running KEGG Pathways...")
  kegg_res <- compareCluster(
    geneCluster = query_list,
    fun = "enrichKEGG",
    universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    # Translate KEGG Entrez IDs back to Symbols for readable CSVs and plots
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    # --- FIXED CSV SAVING: Replace slashes with commas ---
    kegg_df <- as.data.frame(kegg_res)
    kegg_df$geneID <- str_replace_all(kegg_df$geneID, "/", ", ")
    write.csv(kegg_df, file.path(out_dir, paste0(clean_name, "_KEGG_results.csv")), row.names = FALSE)
    # -----------------------------------------------------
    
    # Generate native clusterProfiler DotPlot
    p_kegg <- dotplot(kegg_res, showCategory = 10) +
      ggtitle(paste("KEGG Pathways -", tools::toTitleCase(str_replace_all(clean_name, "_", " ")))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(filename = file.path(out_dir, paste0(clean_name, "_KEGG_bubbleplot.png")), plot = p_kegg, width = 10, height = 7, dpi = 300)
  } else {
    message("   - No significant KEGG pathways found.")
  }
  
  message(paste("   - Success! Saved separate GO and KEGG plots to:", out_dir))
}
# ==============================================================================
# 2. Automated Execution Loop with Automated Logging
# ==============================================================================

base_subset_dir <- "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster"

# --- 1. SET UP THE LOG FILE ---
# Create a logs folder inside the subset_cluster directory if it doesn't exist
log_dir <- file.path(base_subset_dir, "logs")
if (!dir.exists(log_dir)) { dir.create(log_dir, recursive = TRUE) }

# Create a unique filename using the current date and time
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file_path <- file.path(log_dir, paste0("pathway_analysis_", timestamp, ".log"))

# Open the connection and start recording
log_conn <- file(log_file_path, open = "wt")
sink(log_conn, type = "output", split = TRUE) # split = TRUE prints to BOTH the screen and the file
sink(log_conn, type = "message")

message("==================================================================")
message(paste("Log file created at:", log_file_path))
message("=== Starting clusterProfiler Analysis (Separated KEGG/GO:BP) ===")
message("==================================================================")

# --- 2. RUN THE LOOP ---
cell_type_folders <- list.dirs(base_subset_dir, recursive = FALSE, full.names = FALSE)
# Remove the "logs" folder from the list of cell types to process!
cell_type_folders <- cell_type_folders[cell_type_folders != "logs"]

for (cell_type in cell_type_folders) {
  
  dge_dir <- file.path(base_subset_dir, cell_type, "dge_pseudobulk")
  if (!dir.exists(dge_dir)) { next }
  
  dge_files <- list.files(dge_dir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE, recursive = TRUE)
  
  if (length(dge_files) == 0) {
    message(paste("No DGE CSV files found for", cell_type))
    next
  }
  
  for (csv_file in dge_files) {
    # If the analysis crashes on one file, tryCatch ensures the loop continues to the next!
    tryCatch({
      run_pathway_analysis(csv_file, base_subset_dir)
    }, error = function(e) {
      message(paste("   [ERROR] Failed to process", basename(csv_file), ":", conditionMessage(e)))
    })
  }
}

message("\n==================================================================")
message("=== All Pathway Analyses Completed Successfully! ===")
message(paste("Run finished at:", Sys.time()))
message("==================================================================")

# --- 3. CLOSE THE LOG FILE ---
# This is critical! If you don't close the sink, R will keep logging forever.
sink(type = "message")
sink(type = "output")
close(log_conn)