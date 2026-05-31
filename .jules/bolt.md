## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2025-10-24 - Avoiding repetitive subset() calls inside loops in Seurat/R
**Learning:** Subsetting a Seurat object repeatedly inside a `for` loop (e.g., iterating through sample IDs or clusters) causes the entire object to be scanned and sliced on every iteration. This leads to O(N*M) time complexity, where N is cells and M is subsets, becoming a major performance bottleneck for large datasets.
**Action:** Always unconditionally use `SplitObject(seurat_obj, split.by = "column_name")` immediately before the loop to partition the object into a list in a single pass. Then, iterate over the names of the list. Do not conditionally cache the split object, as this can cause stale data bugs if intermediate metadata mutations occur on the parent object.
