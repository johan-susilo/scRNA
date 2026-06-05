## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-03-05 - Avoid repetitive Seurat subsetting inside loops
**Learning:** Repetitively calling `subset(seurat_object, condition)` inside a loop (like iterating through patients or samples) results in massive O(n) CPU overhead because the Seurat object evaluates the entire dataset on every loop iteration.
**Action:** Replace `subset()` inside loops with an unconditional `SplitObject()` before the loop. Crucially, to prevent Out-Of-Memory (OOM) failures from holding duplicate Seurat lists in memory, you MUST explicitly remove the list object immediately after the loop via `rm(); gc()`.
