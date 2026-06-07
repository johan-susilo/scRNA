## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.
## 2026-06-07 - Replacing iterative Seurat subset calls with SplitObject
**Learning:** In Seurat analysis scripts, repeatedly calling `subset(seurat_obj, ...)` inside a loop evaluating factors (like grouping by samples) forces Re-evaluation and memory recopying of massive objects causing severe O(N) CPU and memory delays.
**Action:** Use `SplitObject()` before the loop. It splits the Seurat object efficiently into a list partitioned by the variable in one pass. To prevent Out-Of-Memory (OOM) failures from doubling memory state, immediately drop the object via `rm()` and trigger garbage collection with `gc()` post loop.
