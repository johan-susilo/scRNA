## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.
## 2026-05-28 - Replaced iterative Seurat subsetting with SplitObject
**Learning:** In Seurat workflows, calling `subset(seurat_obj, metadata == condition)` inside loops that iterate over numerous groups (like sample IDs or clusters) severely impacts performance. `subset()` scans the entire object to filter cells and instantiates a new subsetted object on every single iteration, which scales poorly.
**Action:** Replaced repetitive `subset()` calls inside loops with a single `SplitObject(seurat_obj, split.by = "grouping_variable")` operation prior to the loop. Iterating over the resulting pre-split list is incredibly faster, turning an O(n^2) scaling issue into O(1) object creation.
