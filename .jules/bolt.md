## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.

## 2026-05-26 - Avoiding iterative `subset()` calls in Seurat loops
**Learning:** Calling `subset()` on a Seurat object inside a loop (e.g., iterating through sample IDs to generate per-sample plots) is extremely inefficient. It performs deep, state-mutating operations on the large Seurat object redundantly, scaling O(n * S) where S is the cost of a subset operation.
**Action:** Always use `SplitObject(seurat_obj, split.by = "column")` *before* the loop. This partitions the data once (O(1) split cost relatively speaking) into a list of subsets, allowing you to iterate over the pre-split list O(n). This yields massive speedups in per-sample marker visualization and UMAP generation.
