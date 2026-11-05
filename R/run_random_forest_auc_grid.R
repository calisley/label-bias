library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(fs)
library(here)

here::i_am("R/run_random_forest_auc_grid.R")

# Reuse the exact data construction and feature definitions of the regression
# RF grid without executing its command-line entry point.
Sys.setenv(RUBRICS_RF_SOURCE_ONLY = "true")
source(path("R", "run_random_forest_grid.R"))
Sys.unsetenv("RUBRICS_RF_SOURCE_ONLY")

parse_args <- function(args) {
  opts <- list(
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
    out_dir = path("replication_data", "rf_auc_grid"),
    skip_existing = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    if (args[[i]] == "--task-id") {
      i <- i + 1; opts$task_id <- args[[i]]
    } else if (args[[i]] == "--out-dir") {
      i <- i + 1; opts$out_dir <- args[[i]]
    } else if (args[[i]] == "--skip-existing") {
      opts$skip_existing <- TRUE
    } else stop("Unknown argument: ", args[[i]])
    i <- i + 1
  }
  opts$task_id <- as.integer(opts$task_id)
  if (is.na(opts$task_id) || opts$task_id < 1) stop("Missing or invalid task id.")
  opts
}

rf_auc_grid <- function() {
  expand_grid(
    feats_fit = c("feat_ks", "feat_gp", "feat_rb"),
    variant_fit = c("orig", "neut"),
    iter_ran = seq_len(20)
  ) %>% mutate(task_id = row_number(), .before = 1)
}

analyze_rf_auc <- function(
    feats_fit, variant_fit, iter_ran, rf_data, calibrated = TRUE) {
  feats <- rf_data[[feats_fit]]
  df_g <- rf_data[[paste0("df_", variant_fit)]] %>%
    mutate(cf_fold = sample(rep(seq_len(K_outer), length.out = n())))
  df_g_domestic <- df_g %>% filter(!demo_international, !demo_race_w_m)

  p_orig <- p_male <- p_fema <- double(nrow(df_g))
  p_race <- double(nrow(df_g_domestic))
  auc_true_fold <- auc_presented_fold <- auc_race_fold <- double(K_outer)

  for (k in seq_len(K_outer)) {
    train_idx <- which(df_g$cf_fold != k)
    test_idx <- which(df_g$cf_fold == k)
    train_idx_race <- which(df_g_domestic$cf_fold != k)
    test_idx_race <- which(df_g_domestic$cf_fold == k)

    gender_forest <- fit_calibrated_ranger_classif(
      reformulate(feats, "demo_male"), df_g[train_idx, ], nfolds = K_inner,
      calibrate = calibrated
    )
    race_forest <- fit_calibrated_ranger_classif(
      reformulate(feats, "demo_white"), df_g_domestic[train_idx_race, ],
      nfolds = K_inner, calibrate = calibrated
    )
    p_orig[test_idx] <- predict(gender_forest, df_g[test_idx, ])
    p_male[test_idx] <- predict(gender_forest, rf_data$df_male[test_idx, ])
    p_fema[test_idx] <- predict(gender_forest, rf_data$df_fema[test_idx, ])
    p_race[test_idx_race] <- predict(race_forest, df_g_domestic[test_idx_race, ])

    auc_true_fold[k] <- auc_rank(df_g$demo_male[test_idx], p_orig[test_idx])
    auc_presented_fold[k] <- auc_rank(
      c(rep(TRUE, length(test_idx)), rep(FALSE, length(test_idx))),
      c(p_male[test_idx], p_fema[test_idx])
    )
    auc_race_fold[k] <- auc_rank(
      df_g_domestic$demo_white[test_idx_race], p_race[test_idx_race]
    )
  }

  tibble(
    feats_fit,
    variant_fit,
    iter_ran,
    auc_true_gender = mean(auc_true_fold),
    auc_true_gender_std.err = sd(auc_true_fold) / sqrt(K_outer),
    auc_presented_gender = mean(auc_presented_fold),
    auc_presented_gender_std.err = sd(auc_presented_fold) / sqrt(K_outer),
    auc_true_race = mean(auc_race_fold),
    auc_true_race_std.err = sd(auc_race_fold) / sqrt(K_outer),
    gender_pred_diff = mean(p_male - p_fema),
    gender_pred_diff_std.err = sd(p_male - p_fema) / sqrt(length(p_male)),
    gender_pred_diff_orig_male = mean((p_male - p_fema)[df_g$demo_male]),
    gender_pred_diff_orig_fema = mean((p_male - p_fema)[!df_g$demo_male])
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  grid <- rf_auc_grid()
  if (opts$task_id > nrow(grid)) stop("Task id exceeds grid size ", nrow(grid), ".")
  dir_create(opts$out_dir)
  out_file <- path(opts$out_dir, sprintf("rf_auc_%03d.rds", opts$task_id))
  if (opts$skip_existing && file_exists(out_file)) return(invisible(NULL))

  job <- slice(grid, opts$task_id)
  set.seed(5000L + opts$task_id)
  result <- analyze_rf_auc(
    job$feats_fit[[1]], job$variant_fit[[1]], job$iter_ran[[1]], load_rf_data())
  result <- result %>% mutate(task_id = opts$task_id, .before = 1)
  write_rds(result, out_file)
  message("Wrote ", out_file)
}

if (Sys.getenv("RUBRICS_RF_AUC_SOURCE_ONLY", unset = "false") != "true") {
  main()
}
