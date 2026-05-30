## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2024-05-30 - [Seurat SplitObject Pattern for Skin scRNA]
**Learning:** In `scRNA/skin/detail_annotation.R`, repeatedly invoking `subset(obj, orig.ident2 == id)` inside `for` loops across multiple samples incurs significant O(K*N) performance overhead for large Seurat objects.
**Action:** Replace `subset()` within `for` loops by calling Seurat's `SplitObject()` exactly once directly prior to the loop. Note that `SplitObject()` must be re-invoked before later loops if metadata (e.g. `Macro_Lineage`) was mutated in between, preventing stale data mismatches.
