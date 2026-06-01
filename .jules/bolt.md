## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2026-06-01 - Seurat SplitObject() memory overhead
**Learning:** Replacing `subset()` calls inside a loop with a single `SplitObject()` call outside the loop is computationally faster, but inherently stores the entirely duplicated dataset in memory as a list. Failing to remove this list immediately causes catastrophic Out-Of-Memory (OOM) failures in large scRNA-seq datasets.
**Action:** When using `SplitObject()` to optimize loops, always explicitly remove the list object immediately after the loop using `rm(my_split_list); gc()`.
