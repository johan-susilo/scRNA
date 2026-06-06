## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-03-05 - Avoid O(n) CPU overhead and catastrophic OOM with SplitObject in Seurat loops
**Learning:** Repetitively calling `subset()` inside loops in Seurat/R scripts forces O(n) CPU overhead per iteration, making scripts dreadfully slow on large datasets. Using `SplitObject()` unconditionally before loops prevents this, but caches the fully split list object in memory, which causes Out-Of-Memory (OOM) failures for massive datasets downstream.
**Action:** Replace `subset()` inside loops with `SplitObject()` before the loops for a massive CPU boost. However, immediately drop the cached list by using `rm(sample_objects); gc()` right after the loop finishes to avoid catastrophic memory consumption.
