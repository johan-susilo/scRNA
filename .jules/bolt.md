## BOLT'S JOURNAL

## 2025-03-05 - Replacing iterative rbind with dplyr::bind_rows in R
**Learning:** In R, calling `rbind` inside a loop (especially one that iterates many times, like appending single rows to a data frame for every cell annotation vote) causes a complete copy of the data frame in memory on every iteration. This leads to O(n^2) time complexity and severe performance/memory bottlenecks.
**Action:** Always collect rows into a list inside loops (e.g., `my_list[[i]] <- my_data_frame`) and use `dplyr::bind_rows(my_list)` (or `data.table::rbindlist(my_list)`) after the loop finishes. This allocates memory once and is exceptionally faster.
## 2025-03-05 - Seurat Object Subsetting Loop Optimization
**Learning:** Repetitively calling `subset()` on a large Seurat object inside a `for` loop forces R to deep copy and evaluate the entire object on every iteration. This is a severe O(N) performance anti-pattern in single-cell R workflows.
**Action:** Replace `for (id in group) { obj_sub <- subset(obj, ident == id) }` with `split_objs <- SplitObject(obj, split.by = "ident")` followed by iterating over the list elements. This splits the object once, vastly improving memory usage and execution speed.
