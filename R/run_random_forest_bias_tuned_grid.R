library(dplyr)
library(fs)
library(here)

here::i_am("R/run_random_forest_bias_tuned_grid.R")

Sys.setenv(RUBRICS_RF_SOURCE_ONLY = "true")
source(path("R", "run_random_forest_grid.R"))
Sys.unsetenv("RUBRICS_RF_SOURCE_ONLY")

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) y else x
}
task_id <- as.integer((args[1] %||% Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))[[1]])
out_dir <- (args[2] %||% path("replication_data", "rf_bias_tuned_grid"))[[1]]
if (is.na(task_id) || task_id < 1 || task_id > 720) {
  stop("Task id must be 1 through 720.")
}

grid <- rf_grid() %>%
  rename(
    feats_fit = feats,
    variant_fit = variant,
    b_fit = bias,
    iter_ran = iter
  )
job <- slice(grid, .env$task_id)

hyperparameters <- readRDS(
  path("replication_data", "rf_hyperparameters.rds")
) %>%
  filter(
    feats_fit == job$feats_fit[[1]],
    variant_fit == job$variant_fit[[1]],
    b_fit == job$b_fit[[1]]
  )
if (nrow(hyperparameters) != 1) {
  stop("Expected one tuned RF specification for task ", task_id)
}

hyperparameters <- as.list(
  hyperparameters[1, c("mtry", "min_node_size", "sample_fraction")]
)

dir_create(out_dir)
out_file <- path(out_dir, sprintf("rf_%03d.rds", task_id))
set.seed(1000L + task_id)
result <- analyze_rf(
  b = job$b_fit[[1]],
  feats = job$feats_fit[[1]],
  variant = job$variant_fit[[1]],
  iter = job$iter_ran[[1]],
  rf_data = load_rf_data(),
  hyperparameters = hyperparameters
) %>%
  mutate(
    task_id = task_id,
    tuned_mtry = hyperparameters$mtry,
    tuned_min_node_size = hyperparameters$min_node_size,
    tuned_sample_fraction = hyperparameters$sample_fraction,
    .before = 1
  )

write_rds(result, out_file)
message("Wrote ", out_file)
