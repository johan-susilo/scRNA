## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-03-05 - Avoid O(N) repetitive subsetting with SplitObject, but manage memory carefully
**Learning:** Repetitive use of Seurat's `subset()` within a `for` loop on a large object creates a massive O(N) CPU overhead due to constant data filtering. While using `SplitObject()` *before* the loop is a great performance optimization, it duplicates the entire dataset in memory as a list of chunks, leading to catastrophic Out-Of-Memory (OOM) failures on large scRNA-seq datasets if left unchecked.
**Action:** Replace inside-loop `subset()` calls with a single `split_list <- SplitObject(...)` outside the loop to optimize CPU performance. Crucially, explicitly and immediately clear the duplicated split list from memory right after the loop using `rm(split_list); gc()` to avoid OOM crashes.
