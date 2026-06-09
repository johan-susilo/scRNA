## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-03-08 - [Optimize Seurat Subsetting inside Loops]
**Learning:** Performing `subset()` iteratively inside a `for` loop in Seurat is extremely inefficient (O(N*K) time complexity). `SplitObject()` correctly splits the data into a list in O(N) time but creates a fully duplicated dataset in memory, which risks OOM errors.
**Action:** Use `SplitObject()` prior to loops across distinct groups (like sample IDs), iterate over the list elements, and immediately apply `rm(my_split_list); gc()` when the loop completes to reclaim memory space.
