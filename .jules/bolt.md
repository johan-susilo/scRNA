## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2026-06-02 - SplitObject over loop subsetting in Seurat/R
**Learning:** In Seurat/R scRNA-seq scripts, unconditionally calling `subset()` inside a loop causes O(n) CPU overhead and severely degrades performance. `SplitObject()` is much faster because it does the subsetting once upfront, but since it returns a large list, it can easily lead to catastrophic Out-Of-Memory (OOM) failures for large datasets if the list is not immediately deleted.
**Action:** Use `SplitObject()` before looping, iterate over the resulting list, and explicitly clear memory with `rm(); gc()` immediately after the loop to prevent OOM errors.
