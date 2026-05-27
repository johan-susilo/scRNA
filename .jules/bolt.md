## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.
## 2025-05-27 - Replacing repetitive Seurat subset() inside loops with SplitObject()
**Learning:** In Seurat/R scripts, iterating over `unique(object$metadata_col)` and repeatedly calling `subset(object, metadata_col == x)` inside the loop forces an O(N) subset operation on the *entire* object per iteration, leading to O(N * k) time complexity (where k is the number of unique groups).
**Action:** Unconditionally use `SplitObject(object, split.by = "metadata_col")` immediately before the loop, which partitions the object once. Then iterate over the resulting list of objects. (Avoid conditionally caching the split object across different script sections, as intermediate metadata mutations to the parent object will not be reflected).
