library(dplyr)
library(fs)
library(here)

here::i_am("R/aggregate_random_forest_bias_tuned.R")

in_dir <- path("replication_data", "rf_bias_tuned_grid")
paths <- dir_ls(in_dir, glob = "*.rds", type = "file")
paths <- paths[grepl("^rf_[0-9]+\\.rds$", path_file(paths))]
if (length(paths) != 720) {
  stop("Expected 720 RF shards; found ", length(paths), ".")
}

result <- bind_rows(lapply(paths, readRDS)) %>%
  arrange(task_id)
if (!setequal(unique(result$task_id), seq_len(720))) {
  stop("RF shard task IDs are incomplete.")
}

saveRDS(result, path("replication_data", "modeling_rf.rds"))
message("Wrote replication_data/modeling_rf.rds")
