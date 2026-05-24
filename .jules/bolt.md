## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-05-24 - Avoiding `subset()` inside loops in Seurat/R scripts
**Learning:** Calling `subset()` inside a loop over unique sample IDs (or similar metadata columns) on a large Seurat object leads to redundant object parsing, repeated validation overhead, and large memory allocations on every iteration. This pattern slows down R pipelines significantly, especially when generating per-sample plots.
**Action:** Instead of subsetting inside the loop, use `SplitObject(seurat_obj, split.by = "column_name")` once outside the loop to partition the object into a list. Then iterate over the list items. This performs the subsetting optimally and cleanly in a single pass.
