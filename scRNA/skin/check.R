library(Seurat)

# Function to process files from a directory
process_directory <- function(data_dir, stage_name) {
  # Get all RDS files in the directory
  rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)

  # Check if any files were found
  if (length(rds_files) == 0) {
    cat(paste0("No RDS files found in directory: ", data_dir, "\n\n"))
    return(NULL)
  }

  # Create a data frame to store results
  results <- data.frame(
    Stage = character(),
    Sample = character(),
    N_Cells = numeric(),
    N_Genes = numeric(),
    Mean_Counts_Per_Cell = numeric(),
    Total_Counts = numeric(),
    stringsAsFactors = FALSE
  )

  return(list(files = rds_files, results = results, stage = stage_name))
}

# Directories to check
preprocessed_dir <- "/home/johan/data/PMH_scRNA-seq"
processed_dir <- "/home/johan/johan/johan/output/skin_pmh/processed"

# Initialize combined results
all_results <- data.frame(
  Stage = character(),
  Sample = character(),
  N_Cells = numeric(),
  N_Genes = numeric(),
  Mean_Counts_Per_Cell = numeric(),
  Total_Counts = numeric(),
  stringsAsFactors = FALSE
)

# Process both directories
dirs_to_process <- list(
  list(path = preprocessed_dir, stage = "Before Processing"),
  list(path = processed_dir, stage = "After Processing")
)

# Loop through each directory
for (dir_info in dirs_to_process) {
  cat(paste0("\n========== ", dir_info$stage, " ==========\n"))
  cat(paste0("Directory: ", dir_info$path, "\n\n"))

  if (dir_info$stage == "Before Processing") {
    # For raw data, look for 10X format directories
    sample_dirs <- list.dirs(dir_info$path, recursive = FALSE, full.names = TRUE)

    # Filter for directories containing 10X files
    sample_dirs <- sample_dirs[sapply(sample_dirs, function(d) {
      file.exists(file.path(d, "matrix.mtx.gz")) || file.exists(file.path(d, "matrix.mtx"))
    })]

    if (length(sample_dirs) == 0) {
      cat("No 10X format sample directories found.\n")
      next
    }

    # Loop through each sample directory
    for (sample_dir in sample_dirs) {
      sample_name <- basename(sample_dir)
      cat(paste0("Reading 10X data: ", sample_name, "\n"))

      tryCatch({
        # Read 10X format data
        seurat_obj <- Read10X(data.dir = sample_dir)

        # If it's a list (multiple assays), take the first one
        if (is.list(seurat_obj)) {
          counts_matrix <- seurat_obj[[1]]
        } else {
          counts_matrix <- seurat_obj
        }

        # Calculate statistics
        n_cells <- ncol(counts_matrix)
        n_genes <- nrow(counts_matrix)
        total_counts <- sum(counts_matrix)
        mean_counts <- mean(Matrix::colSums(counts_matrix))

        # Add to results
        all_results <<- rbind(all_results, data.frame(
          Stage = dir_info$stage,
          Sample = sample_name,
          N_Cells = n_cells,
          N_Genes = n_genes,
          Mean_Counts_Per_Cell = round(mean_counts, 2),
          Total_Counts = total_counts,
          stringsAsFactors = FALSE
        ))

        cat(paste0("  SUCCESS: ", n_cells, " cells, ", n_genes, " genes\n"))
      }, error = function(e) {
        cat(paste0("  ERROR: ", e$message, "\n"))
      })
    }
  } else {
    # For processed data, read RDS files
    rds_files <- list.files(dir_info$path, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)

    if (length(rds_files) == 0) {
      cat("No RDS files found in this directory.\n")
      next
    }

    # Loop through each RDS file
    for (file in rds_files) {
      # Extract sample name from filename
      sample_name <- basename(file)
      sample_name <- sub("\\.rds$", "", sample_name)
      sample_name <- sub("_processed$", "", sample_name)  # Remove _processed suffix for comparison

      # Read the Seurat object
      cat(paste0("Reading: ", sample_name, "\n"))

      tryCatch({
        seurat_obj <- readRDS(file)

        # Check if it's a Seurat object
        if (!inherits(seurat_obj, "Seurat")) {
          cat(paste0("  SKIPPED: Not a Seurat object (class: ", class(seurat_obj)[1], ")\n"))
          return(NULL)
        }

        # Check if RNA assay exists
        if (!"RNA" %in% names(seurat_obj@assays)) {
          cat(paste0("  SKIPPED: No RNA assay found\n"))
          return(NULL)
        }

        # Extract counts matrix (handle both Seurat v4 and v5)
        if (inherits(seurat_obj@assays$RNA, "Assay5")) {
          # Seurat v5: use LayerData or layers
          counts_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
        } else {
          # Seurat v4 and earlier
          counts_matrix <- seurat_obj@assays$RNA@counts
        }

        # Calculate statistics
        n_cells <- ncol(counts_matrix)
        n_genes <- nrow(counts_matrix)
        total_counts <- sum(counts_matrix)
        mean_counts <- mean(Matrix::colSums(counts_matrix))

        # Add to results
        all_results <<- rbind(all_results, data.frame(
          Stage = dir_info$stage,
          Sample = sample_name,
          N_Cells = n_cells,
          N_Genes = n_genes,
          Mean_Counts_Per_Cell = round(mean_counts, 2),
          Total_Counts = total_counts,
          stringsAsFactors = FALSE
        ))

        cat(paste0("  SUCCESS: ", n_cells, " cells, ", n_genes, " genes\n"))
      }, error = function(e) {
        cat(paste0("  ERROR: ", e$message, "\n"))
      })
    }
  }
}

# Print summary
cat("\n\n========== SUMMARY OF ALL SAMPLES ==========\n\n")
print(all_results, row.names = FALSE)

# Print statistics by stage
cat("\n========== STATISTICS BY STAGE ==========\n")
for (stage in unique(all_results$Stage)) {
  stage_data <- all_results[all_results$Stage == stage, ]
  cat(paste0("\n", stage, ":\n"))
  cat(paste0("  Number of samples: ", nrow(stage_data), "\n"))
  cat(paste0("  Total cells: ", sum(stage_data$N_Cells), "\n"))
  cat(paste0("  Average cells per sample: ", round(mean(stage_data$N_Cells), 2), "\n"))
  cat(paste0("  Average genes per sample: ", round(mean(stage_data$N_Genes), 2), "\n"))
}

# Print comparison if both stages exist
if (length(unique(all_results$Stage)) > 1) {
  cat("\n========== COMPARISON ==========\n")
  before <- all_results[all_results$Stage == "Before Processing", ]
  after <- all_results[all_results$Stage == "After Processing", ]

  if (nrow(before) > 0 && nrow(after) > 0) {
    cat(paste0("Cell retention rate: ", round(sum(after$N_Cells) / sum(before$N_Cells) * 100, 2), "%\n"))
    cat(paste0("Gene retention rate: ", round(mean(after$N_Genes) / mean(before$N_Genes) * 100, 2), "%\n"))
  }
}

