## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-03-05 - Replacing iterative subsetting with SplitObject in Seurat
**Learning:** For performance optimization in Seurat/R scripts, repetitive `subset()` calls inside loops incur heavy O(n) subsetting CPU overhead. While unconditionally using `SplitObject()` before each loop prevents this overhead, it stores the entirely duplicated dataset in memory as a list, which can lead to catastrophic Out-Of-Memory (OOM) failures in large scRNA-seq datasets.
**Action:** Use `SplitObject()` before the loop to split the Seurat object. Iterate over the resulting list. To prevent OOM failures, explicitly remove the list object immediately after the loop using `rm()` and run `gc()`. Do not conditionally cache the split object across script sections to avoid stale data bugs.
