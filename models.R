library(glue)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(glmnet)
library(glmnetUtils)
library(readr)
library(fs)
library(parallel)

here::i_am("models.R")

source(path("R", "utils.R"))

set.seed(1)

K_outer <- 10
K_inner <- 10

EMBEDDING_MODEL <- "text-embedding-3-large"

################################################################################
# Data assembly -----------------------------------------------
#
# This script documents the estimation workflow but cannot run from the
# replication package: it requires confidential source records, extracted
# features, and embeddings. `clean_analysis.R` reproduces every released
# figure from disclosure-safe aggregate outputs.
#
# In the private workflow, four aligned analysis frames are assembled:
# `df_orig` (original materials), `df_neut` (redacted materials), and
# `df_male` / `df_fema` (paired counterfactual presentations). Each contains
# one application per row, outcomes, demographic indicators, structured
# features, and embeddings. Only aggregate results are released.

# Summary outputs -------------------------------------------------------------
# `df` is the confidential primary analysis frame assembled above.

top20 <- df %>% arrange(desc(score)) %>% slice_head(n = round(nrow(df) * .2))
df_domestic <- df %>% filter(!demo_international, !demo_race_w_m)
top20_domestic <- df_domestic %>%
  arrange(desc(score)) %>%
  slice_head(n = round(nrow(df_domestic) * .2))
df_summary <- tibble(
  true_prop_female = 100 * mean(!top20$demo_male),
  true_prop_nonwhite = 100 * mean(!top20_domestic$demo_white),
  true_avg_score_domestic = mean(top20_domestic$score),
  true_avg_score = mean(top20$score),
  true_avg_gpa = mean(top20$gpa, na.rm = TRUE),
  true_avg_gre_quant = mean(top20$gre, na.rm = TRUE),
  true_avg_gre_verbal = mean(top20$gre_verbal, na.rm = TRUE)
)
write_rds(df_summary, path("replication_data", "df_summary.rds"))

# Release histogram-bin counts only; do not release exact individual scores.
df_score_dist <- df %>%
  filter(!is.na(score)) %>%
  transmute(score = floor(score + .5)) %>%
  count(score, name = "n_applicants")
write_rds(df_score_dist, path("replication_data", "df_score_dist.rds"))

# Feature stability is computed inline from the paired counterfactual frames.

feature_names <- grep("^(resume|essay|letter)_", names(df_male), value = TRUE)
feature_names <- feature_names[!grepl("^(resume|essay|letter)_(n|emb)", feature_names)]
feature_names <- feature_names[!grepl("_m$|_mTRUE", feature_names)]
feature_names <- feature_names[!grepl("^letter_.*(_mean|_n_rated)", feature_names)]
feature_matrix <- model.matrix(reformulate(feature_names),
  data = bind_rows(df_male[, feature_names], df_fema[, feature_names]))[, -1, drop = FALSE]

n_applicants <- nrow(df_male)

delta <- feature_matrix[seq_len(n_applicants), , drop = FALSE] -
  feature_matrix[n_applicants + seq_len(n_applicants), , drop = FALSE]

stab_rb <- tibble(
  feature_id = sprintf("F%03d", seq_len(ncol(delta))), n_pairs = colSums(!is.na(delta)),
  mean_delta = colMeans(delta, na.rm = TRUE), sd_delta = apply(delta, 2, sd, na.rm = TRUE),
  change_rate = colMeans(delta != 0, na.rm = TRUE),
  source = case_when(str_detect(colnames(delta), "^resume_") ~ "Resume",
    str_detect(colnames(delta), "^essay_") ~ "Essay", TRUE ~ "Letters")
) %>% mutate(se_delta = sd_delta / sqrt(n_pairs),
  significant = abs(mean_delta) > 2 * se_delta) %>%
  arrange(factor(source, levels = c("Essay", "Resume", "Letters")),
    desc(abs(mean_delta)))

write_rds(stab_rb, path("replication_data", "feature_stability.rds"))

# To produce the alternative-embedding outputs, set this before confidential
# frame assembly: EMBEDDING_MODEL <- "mixedbread-ai/mxbai-embed-large-v1".
# The identical main `run_grid(0)`, marginalized-selection, and `analyze_gender`
# blocks below then write modeling_mixedbread.rds, marginalization_mixedbread.rds,
# and aucs_mixedbread.rds instead of their primary-embedding filenames.

feat_ks <- df_orig %>%
  select(
    matches(
      "^(gre|gpa|toefl|ielts|ranking|joint|letter|essay|resume|transcript)_"
    )
  ) %>%
  colnames()

feat_gp <- feat_ks %>%
  str_subset(
    "^(gre|gpa|toefl|ielts|ranking|joint)_|^(letter|essay|resume|transcript)_(n|emb)"
  ) %>%
  str_subset("^transcript_n_", negate = TRUE)

feat_rb <- feat_ks %>%
  str_subset("^(letter|essay|resume|transcript)_(n|emb)", negate = TRUE)


# Fit simple linear regression for score diffs ----------------------------

m_baseline <- lm(reformulate(c(feat_rb, "demo_male"), "score"), data = df_orig)

# Release only the gender estimate: the full coefficient vector names every
# rubric feature.
gender_ci <- confint(m_baseline, "demo_maleTRUE")

model_stats <- tibble(
  term = "demo_male",
  estimate = coef(m_baseline)[["demo_maleTRUE"]],
  conf_low = gender_ci[[1]],
  conf_high = gender_ci[[2]]
)

write_rds(model_stats, path("replication_data", "baseline_model_stats.rds"))

# Main models ------------------------------------------------------------

analyze <- function(b, feats, variant, iter, alpha = 0, ...) {
  print(
    glue(
      "{Sys.time()}: Bias {b}, features {feats}, variant {variant} iter {iter}",
      " alpha {alpha}"
    )
  )
  
  feat_spec <- feats
  
  feats <- switch(
    feats,
    "feat_ks" = feat_ks,
    "feat_gp" = feat_gp,
    "feat_rb" = feat_rb
  )
  
  df_b <- switch(
    variant,
    "orig" = df_orig,
    "neut" = df_neut
  )
  
  df_b <- df_b %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_male - 1)
    )
  
  prop_male <- mean(df_b$demo_male)
  
  rmse_biased <- double(K_outer)
  rmse_unbiased <- double(K_outer)
  r_sq_biased <- double(K_outer)
  r_sq_unbiased <- double(K_outer)
  prop_bias_k <- double(K_outer)
  
  oos_preds <- double(nrow(df_b))
  gender_marginalized <- double(nrow(df_b))
  audit_diffs <- double(nrow(df_b))
  
  for (k in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx <- which(df_b$cf_fold == k)
    
    m_k <- fit_calibrated_cv_glmnet(
      reformulate(feats, "biased_score"),
      data = df_b[train_idx, ],
      family = "gaussian",
      alpha = alpha,
      nfolds = K_inner
    )
    
    biased_label <- df_b[test_idx, ]$biased_score
    
    unbiased_label <- df_b[test_idx, ]$score
    
    oos_preds[test_idx] <- pred_orig <- predict(
      m_k,
      newdata = df_b[test_idx, ]
    )
    
    rmse_biased[k] <- sqrt(mean((pred_orig - biased_label)^2))
    rmse_unbiased[k] <- sqrt(mean((pred_orig - unbiased_label)^2))
    r_sq_biased[k] <- 1 - var(pred_orig - biased_label) / var(biased_label)
    r_sq_unbiased[k] <- 1 - var(pred_orig - unbiased_label) /
      var(unbiased_label)
    
    pred_male <- predict(
      m_k,
      newdata = df_male[test_idx, ]
    )
    
    pred_fema <- predict(
      m_k,
      newdata = df_fema[test_idx, ]
    )
    
    audit_diffs[test_idx] <- as.numeric(pred_male - pred_fema)
    gender_marginalized[test_idx] <- prop_male * pred_male +
      (1 - prop_male) * pred_fema
    
    R_male_test <- model.matrix(
      reformulate(feats), data = df_male[test_idx, ]
    )[, -1]
    
    R_fema_test <- model.matrix(
      reformulate(feats), data = df_fema[test_idx, ]
    )[, -1]
    
    s_k <- colMeans(R_male_test - R_fema_test)
    
    R_train <- model.matrix(reformulate(feats), data = df_b[train_idx, ])[, -1]
    Y_train  <- df_b[train_idx, ]$score
    
    VarR_k  <- cov(R_train)
    CovRY_k <- cov(R_train, Y_train)
    
    prop_bias_k[k] <- as.numeric(t(s_k) %*% MASS::ginv(VarR_k) %*% CovRY_k)
  }
  
  df_b <- df_b %>%
    mutate(
      oos_preds = as.numeric(oos_preds),
      gender_marginalized = as.numeric(gender_marginalized),
      audit_diffs = as.numeric(audit_diffs)
    )
  
  # Calculate descriptive statistics about the admitted class
  summarize_admitted <- function(df, score_col, selection_rule) {
    df %>%
      filter(
        .data[[score_col]] > quantile(.data[[score_col]], 0.8, na.rm = TRUE)
      ) %>%
      summarize(
        prop_male          = mean(demo_male),
        mean_score         = mean(score),
        mean_gpa           = mean(gpa),
        mean_gre_quant     = mean(gre, na.rm = TRUE),
        mean_gre_verbal    = mean(gre_verbal, na.rm = TRUE),
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
  )
  
  # Return model results
  res <- res %>%
    mutate(
      b_fit = b,
      feats_fit = feat_spec,
      variant_fit = variant,
      alpha_fit = alpha, 
      audit_bias_cf_calib_est = mean(audit_diffs),
      audit_bias_cf_calib_std.err = sd(audit_diffs) / sqrt(length(audit_diffs)),
      rmse_biased_est = mean(rmse_biased),
      rmse_biased_std.err = sd(rmse_biased) / sqrt(K_outer),
      rmse_unbiased_est = mean(rmse_unbiased),
      rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer),
      prop_baseline_bias_est     = mean(prop_bias_k),
      prop_baseline_bias_std.err = sd(prop_bias_k) / sqrt(K_outer)
    )
  
  if (feat_spec == "feat_gp") {
    res_gm <- summarize_admitted(
      df = df_b,
      score_col = "gender_marginalized",
      selection_rule = "gender_marginalized_prediction"
    ) %>%
      #' For the gender marginalized spec, return zero, as the scores are
      #' constant for men and women by construction
      mutate(
        audit_bias_cf_calib_est = 0,
        audit_bias_cf_calib_std.err = 0,
        b_fit = b,
        feats_fit = feat_spec,
        variant_fit = variant
      )
    
    res <- bind_rows(res, res_gm)
  }
  
  res
  
}

run_grid <- function(alpha_val) {
  param_grid <- expand_grid(
    feats   = c("feat_rb", "feat_gp", "feat_ks"),
    variant = c("orig", "neut"),
    bias    = c(0, 0.5, 1, 1.5, 2, 2.5),
    alpha   = alpha_val,
    iter    = seq(20)
  )
  res_list <- mclapply(seq_len(nrow(param_grid)), function(i) {
    set.seed(600000L + i)
    RhpcBLASctl::blas_set_num_threads(1)
    row <- param_grid[i, ]
    analyze(row$bias, row$feats, row$variant, row$iter, row$alpha)
  }, mc.cores = 16)
  bind_rows(res_list)
}

# alpha = 0
res <- run_grid(0)
res %>%
  filter(selection_rule == "original_prediction") %>%
  write_rds(path("replication_data", "modeling.rds"))

res %>%
  filter(selection_rule == "gender_marginalized_prediction") %>%
  write_rds(path("replication_data", "marginalization.rds"))

# alpha = 1
res_lasso <- run_grid(1)
res_lasso %>%
  filter(selection_rule == "original_prediction") %>%
  write_rds(path("replication_data", "modeling_lasso.rds"))
res_lasso %>%
  filter(selection_rule == "gender_marginalized_prediction") %>%
  write_rds(path("replication_data", "marginalization_lasso.rds"))

# Main Models -- AUC ------------------------------------------------------

analyze_gender <- function(feats_name, variant, iter, ...) {
  print(
    glue(
      "{Sys.time()}: Gender, features {feats_name}, variant {variant}, iter {iter}"
    )
  )
  
  feats <- switch(
    feats_name,
    "feat_ks" = feat_ks,
    "feat_gp" = feat_gp,
    "feat_rb" = feat_rb
  )
  
  df_g <- switch(
    variant,
    "orig" = df_orig,
    "neut" = df_neut
  )
  
  df_g <- df_g %>%
    mutate(cf_fold = sample(rep(1:K_outer, length.out = n())))
  
  df_g_domestic <- df_g %>% 
    filter(!demo_international, !demo_race_w_m)
  
  p_orig <- p_male <- p_fema <- double(nrow(df_g))
  p_race <- double(nrow(df_g_domestic))
  auc_true_fold <- auc_presented_fold <- double(K_outer)
  auc_true_race_fold <- double(K_outer)
  
  for (k in 1:K_outer) {
    train_idx <- which(df_g$cf_fold != k)
    test_idx  <- which(df_g$cf_fold == k)
    
    train_idx_race <- which(df_g_domestic$cf_fold != k)
    test_idx_race  <- which(df_g_domestic$cf_fold == k)
    
    m_k <- fit_calibrated_cv_glmnet(
      reformulate(feats, "demo_male"),
      data   = df_g[train_idx, ],
      family = "binomial",
      alpha  = 0,
      nfolds = K_inner
    )
    
    m_k_race <- fit_calibrated_cv_glmnet(
      reformulate(feats, "demo_white"),
      data   = df_g_domestic[train_idx_race, ],
      family = "binomial",
      alpha  = 0,
      nfolds = K_inner
    )
    
    p_orig[test_idx] <- predict(
      m_k, newdata = df_g[test_idx, ],   type = "response"
    )
    p_male[test_idx] <- predict(
      m_k, newdata = df_male[test_idx, ], type = "response"
    )
    p_fema[test_idx] <- predict(
      m_k, newdata = df_fema[test_idx, ], type = "response"
    )
    p_race[test_idx_race] <- predict(
      m_k_race, newdata = df_g_domestic[test_idx_race, ], type = "response"
    )
    
    auc_true_fold[k] <- auc_rank(df_g$demo_male[test_idx], p_orig[test_idx])
    auc_presented_fold[k] <- auc_rank(
      c(rep(TRUE, length(test_idx)), rep(FALSE, length(test_idx))),
      c(p_male[test_idx], p_fema[test_idx])
    )
    auc_true_race_fold[k] <- auc_rank(
      df_g_domestic$demo_white[test_idx_race],p_race[test_idx_race]
    )
  }
  
  tibble(
    feats_fit    = feats_name,
    iter_ran     = iter,
    variant_fit  = variant,
    auc_true_gender           = mean(auc_true_fold),
    auc_true_gender_std.err   = sd(auc_true_fold) / sqrt(K_outer),
    auc_presented_gender      = mean(auc_presented_fold),
    auc_presented_gender_std.err = sd(auc_presented_fold) / sqrt(K_outer),
    auc_true_race = mean(auc_true_race_fold),
    auc_true_race_std.err = sd(auc_true_race_fold) / sqrt(K_outer),
    gender_pred_diff          = mean(p_male - p_fema),
    gender_pred_diff_std.err  = sd(p_male - p_fema) / sqrt(length(p_male)),
    gender_pred_diff_orig_male = mean((p_male - p_fema)[df_g$demo_male]),
    gender_pred_diff_orig_fema = mean((p_male - p_fema)[!df_g$demo_male])
  )
}

param_grid_gender <- expand_grid(
  feats_name = c("feat_ks", "feat_gp", "feat_rb"),
  variant = c("orig", "neut"),
  iter = seq(20)
)

res_gender_list <- mclapply(seq_len(nrow(param_grid_gender)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- param_grid_gender[i, ]
  analyze_gender(row$feats_name, row$variant, row$iter)
}, mc.cores = 12)

res_gender <- bind_rows(res_gender_list)

res_gender %>%
  write_rds(path("replication_data", "aucs.rds"))

# Additional Models ------------------------------------------------------

analyze_ks <- function(b, iter) {
  print(glue("{Sys.time()}: bias {b}, iter: {iter}"))
  
  df_b <- df_orig %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_male - 1)
    )
  
  mse_ks <- double(K_outer)
  mse_rb <- double(K_outer)
  resid_var_ks <- double(K_outer)
  resid_var_rb <- double(K_outer)
  cov_ks <-  double(K_outer)
  
  preds_ks <- double(nrow(df_b))
  preds_rb <- double(nrow(df_b))
  preds_rb_true <- double(nrow(df_b))
  
  true_score <- df_b$score
  biased_score <- df_b$biased_score
  
  for (k in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx <- which(df_b$cf_fold == k)
    
    m_ks <- fit_calibrated_cv_glmnet(
      reformulate(feat_ks, "biased_score"),
      data = df_b[train_idx, ],
      family = "gaussian",
      alpha = 0,
      nfolds = K_inner
    )
    
    m_rb <- fit_calibrated_cv_glmnet(
      reformulate(feat_rb, "biased_score"),
      data = df_b[train_idx, ],
      family = "gaussian",
      alpha = 0,
      nfolds = K_inner
    )
    
    m_rb_true <- fit_calibrated_cv_glmnet(
      reformulate(feat_rb, "score"),
      data = df_b[train_idx, ],
      family = "gaussian",
      alpha = 0,
      nfolds = K_inner
    )
    
    preds_ks[test_idx] <- predict(
      m_ks,
      newdata = df_b[test_idx, ]
    )
    
    preds_rb[test_idx] <- predict(
      m_rb,
      newdata = df_b[test_idx, ]
    )
    
    preds_rb_true[test_idx] <- predict(
      m_rb_true,
      newdata = df_b[test_idx, ]
    )
    
    mse_ks[k] <- mean((preds_ks[test_idx] - true_score[test_idx])^2)
    mse_rb[k] <- mean((preds_rb[test_idx] - true_score[test_idx])^2)
    
    resid_var_ks[k] <- mean((preds_ks[test_idx] - biased_score[test_idx])^2)
    resid_var_rb[k] <- mean((preds_rb[test_idx] - biased_score[test_idx])^2)
  }
  
  for (k in seq(K_outer)) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx <- which(df_b$cf_fold == k)
    
    preds_ks_rb <- df_b %>%
      mutate(preds_ks = preds_ks) %>%
      slice(train_idx) %>%
      fit_calibrated_cv_glmnet(
        reformulate(feat_rb, "preds_ks"),
        data = .,
        family = "gaussian",
        alpha = 0,
        nfolds = K_inner
      ) %>%
      predict(., newdata = df_b[test_idx, ])
    
    cov_ks[k] <- cov(
      preds_ks[test_idx] - preds_ks_rb,
      true_score[test_idx] - preds_rb_true[test_idx]
    )
  }
  
  tibble(
    iter_ran = iter,
    mse_ks_est = mean(mse_ks),
    mse_ks_std.err = sd(mse_ks) / sqrt(K_outer),
    mse_rb_est = mean(mse_rb),
    mse_rb_std.err = sd(mse_rb) / sqrt(K_outer),
    resid_var_ks_est = mean(resid_var_ks),
    resid_var_ks_std.err = sd(resid_var_ks) / sqrt(K_outer),
    resid_var_rb_est = mean(resid_var_rb),
    resid_var_rb_std.err = sd(resid_var_rb) / sqrt(K_outer),
    cov_ks_est = mean(cov_ks),
    cov_ks_std.err = sd(cov_ks) / sqrt(K_outer)
  )
}

res_ks <- expand_grid(
  bias = c(0, 0.5, 1, 1.5, 2, 2.5),
  iter = seq(20)
) %>%
  rowwise() %>%
  mutate(res = list(analyze_ks(bias, iter))) %>%
  ungroup() %>%
  unnest(res)

res_ks %>%
  write_rds(path("replication_data", "ks_theory.rds"))

# Additional Analyses -----------------------------------------------------

## Race instead of Gender
analyze_race <- function(b, feats, variant, iter, ...) {
  print(
    glue(
      "{Sys.time()}: Race, bias {b}, features {feats}, variant {variant}, iter {iter}"
    )
  )
  
  feat_spec <- feats
  
  feats <- switch(
    feats,
    "feat_ks" = feat_ks,
    "feat_gp" = feat_gp,
    "feat_rb" = feat_rb
  )
  
  df_b <- switch(
    variant,
    "orig" = df_orig,
    "neut" = df_neut
  )
  
  df_b <- df_b %>% 
    filter(!demo_international, !demo_race_w_m)
  
  df_b <- df_b %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      # +b for white applicants, -b for non-white applicants
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_white - 1)
    )
  
  rmse_biased   <- double(K_outer)
  rmse_unbiased <- double(K_outer)
  oos_preds     <- double(nrow(df_b))
  
  for (k in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx  <- which(df_b$cf_fold == k)
    
    m_k <- fit_calibrated_cv_glmnet(
      reformulate(feats, "biased_score"),
      data   = df_b[train_idx, ],
      family = "gaussian",
      alpha  = 0,
      nfolds = K_inner
    )
    
    biased_label   <- df_b[test_idx, ]$biased_score
    unbiased_label <- df_b[test_idx, ]$score
    
    oos_preds[test_idx] <- pred_orig <- predict(m_k, newdata = df_b[test_idx, ])
    
    rmse_biased[k]   <- sqrt(mean((pred_orig - biased_label)^2))
    rmse_unbiased[k] <- sqrt(mean((pred_orig - unbiased_label)^2))
  }
  
  df_b <- mutate(df_b, oos_preds = as.numeric(oos_preds))
  
  df_b %>%
    filter(oos_preds > quantile(oos_preds, 0.8, na.rm = TRUE)) %>%
    summarize(
      prop_nonwhite = mean(!demo_white),
      prop_male       = mean(demo_male),
      mean_score      = mean(score),
      mean_gpa        = mean(gpa),
      mean_gre_quant  = mean(gre, na.rm = TRUE),
      mean_gre_verbal = mean(gre_verbal, na.rm = TRUE)
    ) %>%
    mutate(
      b_fit                 = b,
      feats_fit             = feat_spec,
      variant_fit           = variant,
      iter_ran              = iter,
      rmse_biased_est       = mean(rmse_biased),
      rmse_biased_std.err   = sd(rmse_biased) / sqrt(K_outer),
      rmse_unbiased_est     = mean(rmse_unbiased),
      rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer)
    )
}

param_grid_race <- expand_grid(
  feats   = c("feat_ks", "feat_rb", "feat_gp"),
  variant = c("orig"),
  bias    = c(0, 0.5, 1, 1.5, 2, 2.5),
  iter    = seq(20)
)

res_race_list <- mclapply(seq_len(nrow(param_grid_race)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- param_grid_race[i, ]
  analyze_race(row$bias, row$feats, row$variant, row$iter)
}, mc.cores = 16)

res_race <- bind_rows(res_race_list)

res_race %>%
  write_rds(path("replication_data", "modeling_race.rds"))

# Random Forest -----------------------------------------------------------

## Fitting code is released in R/run_random_forest_bias_tuned_grid.R.

# Orthogonalization -------------------------------------------------------

analyze_orth <- function(b, variant, iter, alpha = 0) {
  print(glue("{Sys.time()}: Orth bias {b}, variant {variant}, iter {iter}"))

  feats <- feat_gp
  df_b <- switch(
    variant,
    "orig" = df_orig,
    "neut" = df_neut
  )

  df_b <- df_b %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_male - 1)
    ) %>%
    # Dropped due to fit stability issues once features are orthogonalized
    select(
      -toefl_under_80,
      -toefl_under_90,
      -toefl_under_100,
      -ielts_under_65,
      -gpa_m,
      -joint_hbs_joint_m,
      -joint_hls_joint_m
    )

  feats <- setdiff(feats, c(
    "toefl_under_80", "toefl_under_90", "toefl_under_100",
    "ielts_under_65", "gpa_m",
    "joint_hbs_joint_m", "joint_hls_joint_m"
  ))

  oos_preds_orth   <- double(nrow(df_b))
  audit_diffs_orth <- double(nrow(df_b))

  for (k in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx  <- which(df_b$cf_fold == k)

    train_group <- df_b$demo_male[train_idx]
    group_means <- df_b[train_idx, feats] %>%
      mutate(demo_male = train_group) %>%
      group_by(demo_male) %>%
      summarize(across(everything(), mean), .groups = "drop")

    center_by_group <- function(data) {
      data <- data %>% mutate(across(all_of(feats), as.numeric))
      idx_true  <- which(data$demo_male)
      idx_false <- which(!data$demo_male)
      means_true  <- group_means %>% filter(demo_male) %>% select(-demo_male)
      means_false <- group_means %>% filter(!demo_male) %>% select(-demo_male)
      data[idx_true, feats]  <- sweep(as.matrix(data[idx_true, feats]), 2,
                                      unlist(means_true), "-")
      data[idx_false, feats] <- sweep(as.matrix(data[idx_false, feats]), 2,
                                      unlist(means_false), "-")
      data
    }

    df_b_orth    <- center_by_group(df_b)
    df_male_orth <- center_by_group(df_male)
    df_fema_orth <- center_by_group(df_fema)

    m_k_orth <- fit_calibrated_cv_glmnet(
      reformulate(feats, "biased_score"),
      data = df_b_orth[train_idx, ],
      family = "gaussian",
      alpha = alpha,
      nfolds = K_inner
    )

    p_oos  <- predict(m_k_orth, newdata = df_b_orth[test_idx, ])
    p_male <- predict(m_k_orth, newdata = df_male_orth[test_idx, ])
    p_fema <- predict(m_k_orth, newdata = df_fema_orth[test_idx, ])

    # No longer a hard stop() -- the known destabilizing columns are
    # removed, so this should be rare/absent. If it still fires, warn
    # (don't kill the parallel worker) but let it surface as a real NA
    # in the final audit estimate below rather than being averaged away
    if (any(is.na(p_oos)) || any(is.na(p_male)) || any(is.na(p_fema))) {
      message(glue(
        "NA predictions after stability fix: bias {b}, variant {variant}, ",
        "iter {iter}, fold {k}, n_na={sum(is.na(p_oos))}"
      ))
    }

    oos_preds_orth[test_idx]   <- p_oos
    audit_diffs_orth[test_idx] <- as.numeric(p_male - p_fema)
  }

  df_b <- mutate(df_b, oos_preds_orth = as.numeric(oos_preds_orth))

  df_b %>%
    filter(oos_preds_orth > quantile(oos_preds_orth, 0.8, na.rm = TRUE)) %>%
    summarize(
      prop_male       = mean(demo_male),
      mean_score      = mean(score),
      mean_gpa        = mean(gpa),
      mean_gre_quant  = mean(gre),
      mean_gre_verbal = mean(gre_verbal),
      .groups = "drop"
    ) %>%
    mutate(
      selection_rule = "orthogonalized_prediction",
      b_fit = b,
      feats_fit = "feat_gp",
      variant_fit = variant,
      alpha_fit = alpha,
      iter_ran = iter,
      n_na_audit = sum(is.na(audit_diffs_orth)),
      audit_bias_cf_calib_est = mean(audit_diffs_orth),
      audit_bias_cf_calib_std.err = sd(audit_diffs_orth) / sqrt(length(audit_diffs_orth))
    )
}

param_grid <- expand_grid(
  variant = c("orig", "neut"),
  bias    = c(0, 0.5, 1, 1.5, 2, 2.5),
  iter    = seq(20)
)

res_list <- mclapply(seq_len(nrow(param_grid)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- param_grid[i, ]
  analyze_orth(row$bias, row$variant, row$iter)
}, mc.cores = 16)

res_orth <- bind_rows(res_list)

res_orth %>% write_rds(path("replication_data", "orthogonal.rds"))

# K sweep -----------------------------------------------------------------

analyze_ksweep <- function(k, b, iter) {
  print(glue("{Sys.time()}: k-sweep k={k}, bias={b}, iter={iter}"))
  
  feats_k <- sample(feat_rb, size = k)
  
  df_b <- df_orig %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      biased_score = score + rnorm(n(), b, 1) * (2 * demo_male - 1)
    )
  
  rmse_unbiased <- double(K_outer)
  oos_preds     <- double(nrow(df_b))
  audit_diffs   <- double(nrow(df_b))
  
  for (fold in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != fold)
    test_idx  <- which(df_b$cf_fold == fold)
    
    m_k <- fit_calibrated_cv_glmnet(
      reformulate(feats_k, "biased_score"),
      data   = df_b[train_idx, ],
      family = "gaussian",
      alpha  = 0,
      nfolds = K_inner
    )
    
    unbiased_label <- df_b[test_idx, ]$score
    
    oos_preds[test_idx] <- pred_orig <- predict(m_k, newdata = df_b[test_idx,
    ])
    rmse_unbiased[fold] <- sqrt(mean((pred_orig - unbiased_label)^2))
    
    pred_male <- predict(m_k, newdata = df_male[test_idx, ])
    pred_fema <- predict(m_k, newdata = df_fema[test_idx, ])
    audit_diffs[test_idx] <- as.numeric(pred_male - pred_fema)
  }
  
  df_b <- mutate(df_b, oos_preds = as.numeric(oos_preds))
  
  df_b %>%
    filter(oos_preds > quantile(oos_preds, 0.8, na.rm = TRUE)) %>%
    summarize(
      prop_male       = mean(demo_male),
      mean_score      = mean(score),
      mean_gpa        = mean(gpa),
      mean_gre_quant  = mean(gre, na.rm = TRUE),
      mean_gre_verbal = mean(gre_verbal, na.rm = TRUE)
    ) %>%
    mutate(
      k_fit    = k,
      b_fit    = b,
      iter_ran = iter,
      feats_used = list(feats_k),
      audit_bias_cf_calib_est     = mean(audit_diffs),
      audit_bias_cf_calib_std.err = sd(audit_diffs) /
        sqrt(length(audit_diffs)),
      rmse_unbiased_est     = mean(rmse_unbiased),
      rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer)
    )
}

k_grid <- unique(round(length(feat_rb) * c(0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8,
                                           1)))

param_grid_ksweep <- expand_grid(
  k    = k_grid,
  b    = c(0, .5, 1, 1.5, 2, 2.5),   # unbiased + a clearly-biased point; add more if runtime
  iter = seq(20)      # resample which k items get kept each iter
)

res_ksweep <- mclapply(seq_len(nrow(param_grid_ksweep)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- param_grid_ksweep[i, ]
  analyze_ksweep(row$k, row$b, row$iter)
}, mc.cores = 16) %>%
  bind_rows()

res_ksweep %>%
  write_rds(path("replication_data", "modeling_feature_sweep.rds"))


# Bias on gender-correlated feature.  -------------------------------------

baseline <- mean(df$demo_male)

logical_cols <- names(df)[sapply(df, is.logical)]
logical_cols <- setdiff(logical_cols, grep("demo_", logical_cols, value = TRUE, ignore.case = TRUE))

res <- map_dfr(logical_cols, function(col) {
  x <- df[[col]]
  n_true <- sum(x, na.rm = TRUE)
  if (n_true < 30 || n_true > nrow(df) - 30) return(NULL)
  prop_male_true <- mean(df$demo_male[x], na.rm = TRUE)
  tibble(col = col, n_true = n_true, prop_male_true = prop_male_true, skew = prop_male_true - baseline)
})

res %>%
  mutate(score = skew * sqrt(n_true)) %>%
  arrange(desc(score)) %>%
  print(n = 40)

cat("baseline prop male:", baseline, "\n")

#' Business majors are disproportionately male (.615 of business majors are men)
#' Adding bias favoring business majors

analyze_correlated <- function(b, feats, variant, bias_col, iter, alpha = 0, ...) {
  print(
    glue(
      "{Sys.time()}: Correlated bias {b}, features {feats}, variant {variant},",
      " bias_col {bias_col}, iter {iter}"
    )
  )
  
  feat_spec <- feats
  
  feats <- switch(
    feats,
    "feat_ks" = feat_ks,
    "feat_gp" = feat_gp,
    "feat_rb" = feat_rb
  )
  
  df_b <- switch(
    variant,
    "orig" = df_orig,
    "neut" = df_neut
  )
  
  df_b <- df_b %>%
    mutate(
      cf_fold = sample(rep(1:K_outer, length.out = n())),
      bias_group_composite = transcript_major_business |
        transcript_major_engineering | transcript_major_economics,
      biased_score = score + rnorm(n(), b, 1) * (2 * bias_group_composite - 1),
    )
  
  bias_gpa_cor <- cor(as.numeric(df_b[[bias_col]]), as.numeric(df_b$demo_male))
  
  exempt_vars <- c(
    "transcript_major_business",
    "transcript_major_engineering",
    "transcript_major_economics"
  )
  
  penalty_vec <- as.numeric(!(feats %in% exempt_vars))  # 0 for exempt, 1 otherwise
  
  rmse_unbiased <- double(K_outer)
  oos_preds     <- double(nrow(df_b))
  audit_diffs   <- double(nrow(df_b))
  
  for (k in 1:K_outer) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx  <- which(df_b$cf_fold == k)
    
    m_k <- fit_calibrated_cv_glmnet(
      reformulate(feats, "biased_score"),
      data   = df_b[train_idx, ],
      family = "gaussian",
      alpha  = alpha,
      nfolds = K_inner,
      penalty.factor = penalty_vec
    )
    
    unbiased_label <- df_b[test_idx, ]$score
    
    oos_preds[test_idx] <- pred_orig <- predict(m_k, newdata = df_b[test_idx, ])
    rmse_unbiased[k] <- sqrt(mean((pred_orig - unbiased_label)^2))
    
    pred_male <- predict(m_k, newdata = df_male[test_idx, ])
    pred_fema <- predict(m_k, newdata = df_fema[test_idx, ])
    audit_diffs[test_idx] <- as.numeric(pred_male - pred_fema)
  }
  
  df_b <- mutate(df_b, oos_preds = as.numeric(oos_preds))
  
  mean_score_pred_male_est <- mean(df_b$oos_preds[df_b$demo_male])
  mean_score_pred_fema_est <- mean(df_b$oos_preds[!df_b$demo_male])
  
  df_b %>%
    filter(oos_preds > quantile(oos_preds, 0.8, na.rm = TRUE)) %>%
    summarize(
      prop_male       = mean(demo_male),
      mean_score      = mean(score),
      mean_gpa        = mean(gpa),
      mean_gre_quant  = mean(gre, na.rm = TRUE),
      mean_gre_verbal = mean(gre_verbal, na.rm = TRUE)
    ) %>%
    mutate(
      feats           = feat_spec,
      variant         = variant,
      bias            = b,
      bias_col        = bias_col,
      iter_ran        = iter,
      bias_gpa_cor    = bias_gpa_cor,
      audit_bias_cf_calib_est     = mean(audit_diffs),
      audit_bias_cf_calib_std.err = sd(audit_diffs) / sqrt(length(audit_diffs)),
      rmse_unbiased_est     = mean(rmse_unbiased),
      rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer),
      mean_score_pred_male_est = mean_score_pred_male_est,
      mean_score_pred_fema_est = mean_score_pred_fema_est,
      mean_score_pred_gap_est  = mean_score_pred_male_est - mean_score_pred_fema_est
    )
}

param_grid_correlated <- expand_grid(
  feats    = c("feat_rb", "feat_gp", "feat_ks"),
  variant  = c("orig", "neut"),
  bias_col = "bias_group_composite",
  bias     = c(0, 0.5, 1, 1.5, 2, 2.5),
  iter     = seq(20)
)

res_correlated_list <- mclapply(seq_len(nrow(param_grid_correlated)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- param_grid_correlated[i, ]
  analyze_correlated(row$bias, row$feats, row$variant, row$bias_col, row$iter)
}, mc.cores = 16)

res_correlated <- bind_rows(res_correlated_list)

res_correlated %>%
  write_rds(path("replication_data", "modeling_attribute_bias.rds"))


# First-year GPA outcome ------------------------------------------------------
# `df_first_year` is assembled privately with the same predictors as the main
# analysis and an observed first-year GPA column, `grade_outcome`. Unlike the
# audit analyses, it has no paired presentation frames and therefore computes
# no counterfactual audit differences. These applications are from a different 
# year.

top20_first_year <- df_first_year %>%
  arrange(desc(grade_outcome)) %>%
  slice_head(n = round(nrow(df_first_year) * .2))

df_summary_first_year <- tibble(
  true_prop_female = mean(!top20_first_year$demo_male),
  true_avg_score = mean(top20_first_year$grade_outcome),
  true_avg_gpa = mean(top20_first_year$gpa, na.rm = TRUE),
  true_avg_gre_quant = mean(top20_first_year$gre, na.rm = TRUE),
  true_avg_gre_verbal = mean(top20_first_year$gre_verbal, na.rm = TRUE)
)
write_rds(df_summary_first_year, path("replication_data", "df_summary_first_year.rds"))

analyze_first_year <- function(b, feats, iter, alpha = 0) {
  feat_spec <- feats
  features <- switch(feats, feat_ks = feat_ks, feat_gp = feat_gp, feat_rb = feat_rb)
  
  df_b <- df_first_year %>% mutate(
    cf_fold = sample(rep(seq_len(K_outer), length.out = n())),
    biased_gpa = grade_outcome + rnorm(n(), b, .5) * (2 * demo_male - 1)
  )
  
  rmse_biased <- rmse_unbiased <- double(K_outer)
  oos_preds <- double(nrow(df_b))
  
  for (k in seq_len(K_outer)) {
    train_idx <- which(df_b$cf_fold != k)
    test_idx <- which(df_b$cf_fold == k)
    model <- fit_calibrated_cv_glmnet(reformulate(features, "biased_gpa"),
                                      data = df_b[train_idx, ], family = "gaussian", alpha = alpha, nfolds = K_inner)
    oos_preds[test_idx] <- predict(model, newdata = df_b[test_idx, ]) %>% as.numeric()
    rmse_biased[k] <- sqrt(mean((oos_preds[test_idx] - df_b$biased_gpa[test_idx])^2))
    rmse_unbiased[k] <- sqrt(mean((oos_preds[test_idx] - df_b$grade_outcome[test_idx])^2))
  }
  
  df_b <- mutate(df_b, oos_preds = oos_preds)
  admitted <- filter(df_b, oos_preds > quantile(oos_preds, .8, na.rm = TRUE))
  
  tibble(
    prop_male = mean(admitted$demo_male), mean_score = mean(admitted$grade_outcome),
    mean_gpa = mean(admitted$gpa, na.rm = TRUE),
    mean_gre_quant = mean(admitted$gre, na.rm = TRUE),
    mean_gre_verbal = mean(admitted$gre_verbal, na.rm = TRUE),
    selection_rule = "original_prediction", admit_score_col = "oos_preds",
    b_fit = b, feats_fit = feat_spec, variant_fit = "orig", alpha_fit = alpha,
    iter = iter, rmse_biased_est = mean(rmse_biased),
    rmse_biased_std.err = sd(rmse_biased) / sqrt(K_outer),
    rmse_unbiased_est = mean(rmse_unbiased),
    rmse_unbiased_std.err = sd(rmse_unbiased) / sqrt(K_outer)
  )
}

first_year_grid <- expand_grid(
  feats = c("feat_rb", "feat_gp", "feat_ks"), b = c(0, .2, .4, .6, .8, 1), iter = seq_len(20)
)

res_first_year <- mclapply(seq_len(nrow(first_year_grid)), function(i) {
  set.seed(600000L + i)
  RhpcBLASctl::blas_set_num_threads(1)
  row <- first_year_grid[i, ]
  analyze_first_year(row$b, row$feats, row$iter)
}, mc.cores = 16) %>% bind_rows()

write_rds(res_first_year, path("replication_data", "modeling_first_year.rds"))
