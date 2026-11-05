library(glue)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(fs)
library(here)

here::i_am("R/run_random_forest_grid.R")

source(path("R", "utils.R"))

K_outer <- 10
K_inner <- 10

parse_args <- function(args) {
  opts <- list(
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
    out_dir = path("replication_data", "rf_grid"),
    hyperparameters = NULL,
    skip_existing = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--task-id") {
      i <- i + 1
      opts$task_id <- args[[i]]
    } else if (arg == "--out-dir") {
      i <- i + 1
      opts$out_dir <- args[[i]]
    } else if (arg == "--hyperparameters") {
      i <- i + 1
      opts$hyperparameters <- args[[i]]
    } else if (arg == "--skip-existing") {
      opts$skip_existing <- TRUE
    } else {
      stop("Unknown argument: ", arg)
    }

    i <- i + 1
  }

  if (is.na(opts$task_id) || opts$task_id == "") {
    stop("Missing task id. Use --task-id or set SLURM_ARRAY_TASK_ID.")
  }

  opts$task_id <- as.integer(opts$task_id)
  if (is.na(opts$task_id) || opts$task_id < 1) {
    stop("Task id must be a positive integer.")
  }

  opts
}

rf_grid <- function() {
  expand_grid(
    feats = c("feat_ks", "feat_rb", "feat_gp"),
    variant = c("orig", "neut"),
    bias = c(0, 0.5, 1, 1.5, 2, 2.5),
    iter = seq(20)
  ) %>%
    mutate(task_id = row_number())
}

load_rf_data <- function(n_embedding_dims = 250) {
  # Requires confidential source records, extracted features, and embeddings.
  # The private workflow returns `df_orig`, `df_neut`, `df_male`, and
  # `df_fema`, plus the feature vectors `feat_ks`, `feat_gp`, `feat_gp_demo`,
  # and `feat_rb`. No record-level input is distributed.
  stop("Random-forest fitting requires confidential source data.")
}
analyze_rf <- function(
    b, feats, variant, iter, rf_data, hyperparameters = NULL,
    calibrated = TRUE) {
  print(
    glue(
      "{Sys.time()}: bias {b}, features {feats}, variant {variant}, iter {iter}"
    )
  )

  feat_spec <- feats

  feats <- switch(
    feats,
    "feat_ks" = rf_data$feat_ks,
    "feat_gp" = rf_data$feat_gp,
    "feat_gp_demo" = rf_data$feat_gp_demo,
    "feat_rb" = rf_data$feat_rb
  )

  df_b <- switch(
    variant,
    "orig" = rf_data$df_orig,
    "neut" = rf_data$df_neut
  )

  df_b <- df_b %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_male - 1)
    )

  rmse_biased <- double(K_outer)
  rmse_unbiased <- double(K_outer)
  oos_preds <- double(nrow(df_b))
  gender_marginalized <- double(nrow(df_b))
  audit_diffs <- double(nrow(df_b))
  pred_male_all <- double(nrow(df_b))
  pred_fema_all <- double(nrow(df_b))
  prop_male <- mean(df_b$demo_male)

  for (k in seq_len(K_outer)) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx <- which(df_b$cf_fold == k)

    m_k <- fit_calibrated_ranger(
      reformulate(feats, "biased_score"),
      data = df_b[train_idx, ],
      nfolds = K_inner,
      mtry = hyperparameters$mtry,
      min.node.size = hyperparameters$min_node_size,
      sample.fraction = hyperparameters$sample_fraction,
      calibrate = calibrated
    )

    biased_label <- df_b[test_idx, ]$biased_score
    unbiased_label <- df_b[test_idx, ]$score

    oos_preds[test_idx] <- pred_orig <- predict(
      m_k,
      newdata = df_b[test_idx, ]
    )
    rmse_biased[k] <- sqrt(mean((pred_orig - biased_label)^2))
    rmse_unbiased[k] <- sqrt(mean((pred_orig - unbiased_label)^2))

    male_counterfactual <- rf_data$df_male[test_idx, ]
    female_counterfactual <- rf_data$df_fema[test_idx, ]
    # When demographic gender is an explicit model feature, it is part of the
    # presentation counterfactual: male-presenting inputs receive TRUE and
    # female-presenting inputs receive FALSE.  Otherwise demo_male remains
    # untouched, as in the baseline document-only audit.
    if (feat_spec == "feat_gp_demo") {
      male_counterfactual$demo_male <- TRUE
      female_counterfactual$demo_male <- FALSE
    }

    pred_male <- predict(m_k, newdata = male_counterfactual)
    pred_fema <- predict(m_k, newdata = female_counterfactual)

    audit_diffs[test_idx] <- as.numeric(pred_male - pred_fema)
    pred_male_all[test_idx] <- as.numeric(pred_male)
    pred_fema_all[test_idx] <- as.numeric(pred_fema)
    gender_marginalized[test_idx] <- prop_male * pred_male +
      (1 - prop_male) * pred_fema
  }

  df_b <- df_b %>%
    mutate(
      oos_preds = as.numeric(oos_preds),
      gender_marginalized = as.numeric(gender_marginalized),
      audit_diffs = as.numeric(audit_diffs),
      pred_male = as.numeric(pred_male_all),
      pred_fema = as.numeric(pred_fema_all)
    )

  summarize_admitted <- function(df, score_col, selection_rule) {
    df %>%
      filter(
        .data[[score_col]] > quantile(.data[[score_col]], 0.8, na.rm = TRUE)
      ) %>%
      summarize(
        prop_male = mean(demo_male),
        mean_score = mean(score),
        mean_gpa = mean(gpa),
        mean_gre_quant = mean(gre, na.rm = TRUE),
        mean_gre_verbal = mean(gre_verbal, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        selection_rule = selection_rule,
        admit_score_col = score_col
      )
  }

  res <- summarize_admitted(
    df = df_b,
    score_col = "oos_preds",
    selection_rule = "original_prediction"
  ) %>%
    mutate(
      b_fit = b,
      feats_fit = feat_spec,
      variant_fit = variant,
      iter_ran = iter,
      audit_bias_cf_calib_est = mean(audit_diffs),
      audit_bias_cf_calib_std.err = sd(audit_diffs) / sqrt(length(audit_diffs)),
      rmse_biased_est = mean(rmse_biased),
      rmse_biased_std.err = sd(rmse_biased) / sqrt(K_outer),
      rmse_unbiased_est = mean(rmse_unbiased),
      rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer)
    )

  if (feat_spec == "feat_gp") {
    res_gm <- summarize_admitted(
      df = df_b,
      score_col = "gender_marginalized",
      selection_rule = "gender_marginalized_prediction"
    ) %>%
      mutate(
        b_fit = b,
        feats_fit = feat_spec,
        variant_fit = variant,
        iter_ran = iter,
        audit_bias_cf_calib_est = 0,
        audit_bias_cf_calib_std.err = 0,
        rmse_biased_est = mean(rmse_biased),
        rmse_biased_std.err = sd(rmse_biased) / sqrt(K_outer),
        rmse_unbiased_est = mean(rmse_unbiased),
        rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer)
      )

    res <- bind_rows(res, res_gm)
  }

  res
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  grid <- rf_grid()

  if (opts$task_id > nrow(grid)) {
    stop("Task id ", opts$task_id, " exceeds grid size ", nrow(grid), ".")
  }

  dir_create(opts$out_dir)

  out_file <- path(opts$out_dir, sprintf("rf_%03d.rds", opts$task_id))
  if (opts$skip_existing && file_exists(out_file)) {
    message("Skipping existing shard: ", out_file)
    return(invisible(NULL))
  }

  job <- slice(grid, opts$task_id)

  hyperparameters <- NULL
  if (!is.null(opts$hyperparameters)) {
    tuned <- read_rds(opts$hyperparameters)
    hyperparameters <- tuned %>%
      filter(feats_fit == job$feats[[1]], variant_fit == job$variant[[1]])
    if ("b_fit" %in% names(hyperparameters)) {
      hyperparameters <- hyperparameters %>%
        filter(b_fit == job$bias[[1]])
    }
    if (nrow(hyperparameters) != 1) {
      stop("Expected one tuned RF specification for ", job$feats[[1]], "/", job$variant[[1]])
    }
    hyperparameters <- as.list(
      hyperparameters[1, c("mtry", "min_node_size", "sample_fraction")]
    )
  }

  # Deterministic per-task seeding so retries reproduce the same shard.
  set.seed(1000L + opts$task_id)

  rf_data <- load_rf_data()
  res <- analyze_rf(
    b = job$bias[[1]],
    feats = job$feats[[1]],
    variant = job$variant[[1]],
    iter = job$iter[[1]],
    rf_data = rf_data,
    hyperparameters = hyperparameters
  ) %>%
    mutate(
      task_id = opts$task_id,
      tuned_mtry = if (is.null(hyperparameters)) NA_integer_ else hyperparameters$mtry,
      tuned_min_node_size = if (is.null(hyperparameters)) NA_integer_ else hyperparameters$min_node_size,
      tuned_sample_fraction = if (is.null(hyperparameters)) NA_real_ else hyperparameters$sample_fraction,
      .before = 1
    )

  write_rds(res, out_file)
  message("Wrote ", out_file)
}

if (!identical(Sys.getenv("RUBRICS_RF_SOURCE_ONLY", unset = "false"), "true")) {
  main()
}
