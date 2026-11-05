library(dplyr)
library(tidyr)
library(fs)
library(here)

here::i_am("R/tune_random_forest_by_bias.R")

Sys.setenv(RUBRICS_RF_SOURCE_ONLY = "true")
source(path("R", "run_random_forest_grid.R"))
Sys.unsetenv("RUBRICS_RF_SOURCE_ONLY")

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) y else x
}
task_id <- as.integer((args[1] %||% Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))[[1]])
out_dir <- (args[2] %||% path("replication_data", "rf_tuning_by_bias"))[[1]]
if (is.na(task_id) || task_id < 1 || task_id > 36) {
  stop("Task id must be 1 through 36.")
}

specs <- expand_grid(
  feats_fit = c("feat_ks", "feat_rb", "feat_gp"),
  variant_fit = c("orig", "neut"),
  b_fit = c(0, 0.5, 1, 1.5, 2, 2.5)
) %>%
  mutate(task_id = row_number(), .before = 1)

job <- slice(specs, .env$task_id)
set.seed(9000L + task_id)
rf_data <- load_rf_data()
df <- rf_data[[paste0("df_", job$variant_fit[[1]])]] %>%
  mutate(biased_score = score +
    rnorm(n(), job$b_fit[[1]], 1) * (2 * demo_male - 1))
features <- rf_data[[job$feats_fit[[1]]]]

candidate_grid <- expand_grid(
  mtry_fraction = c(0.10, 0.33, 0.67),
  min_node_size = c(1L, 5L, 10L),
  sample_fraction = c(0.632, 0.8)
) %>%
  mutate(
    mtry = pmax(1L, pmin(length(features),
      round(mtry_fraction * length(features))))
  )

folds <- sample(rep(seq_len(5), length.out = nrow(df)))
score_candidate <- function(mtry, min_node_size, sample_fraction) {
  fold_mse <- vapply(seq_len(5), function(k) {
    forest <- ranger::ranger(
      reformulate(features, "biased_score"),
      data = df[folds != k, ],
      num.trees = 500,
      mtry = mtry,
      min.node.size = min_node_size,
      sample.fraction = sample_fraction,
      num.threads = 4
    )
    pred <- predict(forest, data = df[folds == k, ])$predictions
    mean((pred - df$biased_score[folds == k]) ^ 2)
  }, numeric(1))
  mean(fold_mse)
}

scored <- candidate_grid %>%
  rowwise() %>%
  mutate(cv_mse = score_candidate(mtry, min_node_size, sample_fraction)) %>%
  ungroup() %>%
  arrange(cv_mse, mtry, min_node_size, sample_fraction)

selected <- slice(scored, 1) %>%
  mutate(
    feats_fit = job$feats_fit[[1]],
    variant_fit = job$variant_fit[[1]],
    b_fit = job$b_fit[[1]],
    task_id = task_id
  )

dir_create(out_dir)
saveRDS(
  list(selected = selected, candidates = scored),
  path(out_dir, sprintf("rf_tune_%02d.rds", task_id))
)
message(
  "Selected mtry=", selected$mtry,
  ", min.node.size=", selected$min_node_size,
  ", sample.fraction=", selected$sample_fraction
)
