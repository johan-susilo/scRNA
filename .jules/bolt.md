## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2026-05-29 - Replacing iterative Seurat subsetting with SplitObject()
**Learning:** Repeatedly calling `subset()` on a large Seurat object inside a loop triggers costly deep copies and memory allocations per iteration, resulting in O(N * S) time complexity and severe performance degradation. However, conditional caching (e.g. `if (!exists(...))`) across different sections leads to stale data bugs.
**Action:** Unconditionally partition the Seurat object once immediately before each loop using `SplitObject(obj, split.by = 'var')` and iterate over the resulting list instead.
