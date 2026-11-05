library(dplyr)
library(fs)
library(here)

here::i_am("R/aggregate_random_forest_tuning_by_bias.R")

in_dir <- path("replication_data", "rf_tuning_by_bias")
paths <- dir_ls(in_dir, glob = "*.rds", type = "file")
paths <- paths[grepl("^rf_tune_[0-9]+\\.rds$", path_file(paths))]
if (length(paths) != 36) {
  stop("Expected 36 RF tuning shards; found ", length(paths), ".")
}

selected <- bind_rows(lapply(paths, function(p) readRDS(p)$selected)) %>%
  arrange(task_id)
if (!setequal(selected$task_id, seq_len(36))) {
  stop("RF tuning shard IDs are incomplete.")
}

saveRDS(selected, path("replication_data", "rf_hyperparameters.rds"))
message("Wrote replication_data/rf_hyperparameters.rds")
