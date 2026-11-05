# R/utils.R
#
# Shared utility functions sourced by models.R and any other
# downstream scripts that need them.
#
# Data-loading functions below require institutional data access and are
# provided as stubs for transparency. See the paper for a full description
# of the data pipeline.
#
# Fully implemented functions:
#   fit_calibrated_cv_glmnet(formula, data, family, alpha, nfolds, ...)
#   predict.calibrated_cv_glmnet(object, newdata, type)
#   auc_rank(labels, scores)
#   inject_bias(df, demographic_col, demographic_val, p, seed)
#   clip_prob(p, eps)

library(here)
library(fs)
library(tidyverse)
library(vroom)

# Data-loading stubs ----------------------------------------------------------
# The following functions require access to the institutional data pipeline
# and cannot be run from the replication package. They are included to
# document the data loading and processing steps. See the paper for a
# full description of the data pipeline.

emb_path <- function(source, model = "text-embedding-3-large") {
  # Returns the file path to the raw embeddings CSV for a given source
  # (essays, resumes, letters, transcripts_structured)
  # and embedding model.
  stop("Requires institutional data access.")
}

audit_emb_path <- function(source, variant, model = "text-embedding-3-large") {
  # Returns the file path to the raw embeddings CSV for a given source and
  # audit variant (original, male, female, neutral).
  stop("Requires institutional data access.")
}

audit_manifest_path <- function(source) {
  # Returns the file path to the audit manifest CSV for a given source,
  # which tracks which application files have been processed.
  stop("Requires institutional data access.")
}

list_file_stems <- function(dir, glob) {
  # Lists the file stems (filenames without extensions) in a directory
  # matching a given glob pattern.
  stop("Requires institutional data access.")
}

list_file_names <- function(dir, glob) {
  # Lists the filenames in a directory matching a given glob pattern.
  stop("Requires institutional data access.")
}

read_embedding_filenames <- function(source, variant,
                                     model = "text-embedding-3-large") {
  # Reads the source_filename column from an embeddings CSV for a given
  # source and audit variant, returning a character vector of filenames.
  stop("Requires institutional data access.")
}

load_audit_complete_ids <- function(sources = c("essays", "resumes", "letters"),
                                    variants = c("original", "male", "female",
                                                 "neutral"),
                                    model = "text-embedding-3-large") {
  # Returns a vector of applicant IDs for which all audit study materials
  # (templates, filled variants, and embeddings) are complete across all
  # sources and variants.
  stop("Requires institutional data access.")
}

load_audit_structured <- function(source, variant) {
  # Loads pre-computed structured rubric features for a given source
  # (essays, resumes, letters) and audit variant (original, male, female,
  # neutral) from the institutional feature cache.
  # Returns: tibble with applicant identifiers and rubric feature columns.
  stop("Requires institutional data access.")
}

join_essays <- function(emb) {
  # Joins raw essay embeddings to applicant identifiers using an institutional
  # file map linking essay filenames to applicants.
  # Returns: emb tibble with applicant identifier column added.
  stop("Requires institutional data access.")
}

join_resumes <- function(emb) {
  # Joins raw resume embeddings to applicant identifiers using an institutional
  # file map linking resume filenames to applicants.
  # Returns: emb tibble with applicant identifier column added.
  stop("Requires institutional data access.")
}

join_letters <- function(emb) {
  # Joins raw letter of recommendation embeddings to applicant identifiers
  # using an institutional file map linking letter filenames to applicants.
  # Returns: emb tibble with applicant identifier column added.
  stop("Requires institutional data access.")
}

join_transcripts <- function(emb) {
  # Joins raw transcript embeddings to applicant identifiers using an
  # institutional file map linking transcript filenames to applicants.
  # Returns: emb tibble with applicant identifier column added.
  stop("Requires institutional data access.")
}

JOIN_FNS <- list(
  # Dispatch table mapping each source to its join function.
  # Requires institutional data access.
  essays                 = join_essays,
  resumes                = join_resumes,
  letters                = join_letters,
  transcripts_structured = join_transcripts
)

load_embeddings <- function(source, n_dims = NULL,
                            model = "text-embedding-3-large",
                            variant = NULL) {
  # Loads embeddings for a given source from the institutional data pipeline.
  # Reads raw embeddings CSVs, joins to applicant identifiers via
  # source-specific file maps, mean-pools across multiple documents per
  # applicant, and slices to n_dims dimensions.
  # Returns: tibble with columns applicant_id, emb_1..emb_n, n_docs.
  stop("Requires institutional data access.")
}

load_swapped_embeddings <- function(source, n_dims = NULL,
                                    model = "text-embedding-3-large",
                                    df) {
  # Loads gender-swapped embeddings for use in the correspondence experiment.
  # For each applicant, returns the embedding of the opposite-gender variant
  # of their materials (male applicants get female variant embeddings and
  # vice versa).
  # Returns: tibble with columns applicant_id, emb_1..emb_n, n_docs.
  stop("Requires institutional data access.")
}


# Fully implemented functions ---------------------------------------------

## inject_bias -----------------------------------------------------------------
# Randomly flips accept -> reject with probability p for applicants where
# demographic_col == demographic_val.
# Adds a new column admit_biased, leaving meta_admitted untouched.

inject_bias <- function(df, demographic_col, demographic_val, p, seed = 42) {
  set.seed(seed)
  n <- nrow(df)
  df %>%
    mutate(score_biased = score +
             rnorm(n, mean = 0, sd = 1) +
             (.data[[demographic_col]] == demographic_val) * (p / 2) -
             (.data[[demographic_col]] != demographic_val) * (p / 2))
}

# Calibrated glmnet models ----------------------------------------------------

clip_prob <- function(p, eps = 1e-6) {
  p <- as.numeric(p)
  pmin(pmax(p, eps), 1 - eps)
}

fit_calibrated_cv_glmnet <- function(
    formula,
    data,
    family = "gaussian",
    alpha = 0,
    nfolds = 5,
    ...) {
  # `...` is passed to cv.glmnet, e.g. penalty.factor in analyze_correlated().
  base_model <- cv.glmnet(
    formula,
    data = data,
    family = family,
    alpha = alpha,
    nfolds = nfolds,
    ...
  )
  
  if (family == "binomial") {
    p_train <- predict(
      base_model,
      newdata = data,
      type = "response",
      s = "lambda.min",
      na.action = na.pass
    ) %>%
      as.numeric() %>%
      clip_prob()
    
    y_name <- all.vars(formula)[1]
    
    calib_df <- tibble(
      y = data[[y_name]],
      logit_p = qlogis(p_train)
    )
    
    calib_model <- glm(
      y ~ logit_p,
      data = calib_df,
      family = binomial()
    )
    
    model_type <- "binomial"
  } else {
    pred_train <- predict(
      base_model,
      newdata = data,
      s = "lambda.min",
      na.action = na.pass
    ) %>%
      as.numeric()
    
    y_name <- all.vars(formula)[1]
    
    calib_df <- tibble(
      y = data[[y_name]],
      pred = pred_train
    )
    
    calib_model <- lm(
      y ~ pred,
      data = calib_df
    )
    
    model_type <- "gaussian"
  }
  
  structure(
    list(
      base_model = base_model,
      calib_model = calib_model,
      family = model_type
    ),
    class = "calibrated_cv_glmnet"
  )
}

predict.calibrated_cv_glmnet <- function(
    object,
    newdata,
    type = "response",
    ...) {
  if (object$family == "binomial") {
    p_raw <- predict(
      object$base_model,
      newdata = newdata,
      type = "response",
      s = "lambda.min",
      na.action = na.pass
    ) %>%
      as.numeric() %>%
      clip_prob()
    
    logit_raw <- qlogis(p_raw)
    
    eta_calib <- as.numeric(cbind(1, logit_raw) %*% coef(object$calib_model))
    
    if (type == "link") {
      return(eta_calib)
    }
    
    if (type == "response") {
      return(plogis(eta_calib))
    }
    
    stop("For binomial models, type must be 'response' or 'link'.")
  }
  
  pred_raw <- predict(
    object$base_model,
    newdata = newdata,
    s = "lambda.min",
    na.action = na.pass
  ) %>%
    as.numeric()
  
  pred_calib <- as.numeric(cbind(1, pred_raw) %*% coef(object$calib_model))
  
  pred_calib
}

auc_rank <- function(labels, scores) {
  as.numeric(pROC::auc(pROC::roc(labels, scores, direction = "<")))
}

fit_calibrated_ranger <- function(
    formula,
    data,
    nfolds = 10,
    num.trees = 500,
    mtry = NULL,
    min.node.size = NULL,
    sample.fraction = NULL,
    calibrate = TRUE) {
  ranger_args <- purrr::compact(list(
    formula = formula,
    data = data,
    num.trees = num.trees,
    num.threads = 4,
    mtry = mtry,
    min.node.size = min.node.size,
    sample.fraction = sample.fraction
  ))
  
  y_name <- all.vars(formula)[1]

  if (!calibrate) {
    return(structure(
      list(base_model = do.call(ranger::ranger, ranger_args)),
      class = "uncalibrated_ranger"
    ))
  }
  
  # Out-of-fold predictions, used only to fit the calibration line -- an
  # uncalibrated forest shrinks toward the mean, which would make it look
  # less biased than glmnet for reasons that have nothing to do with the
  # forest itself.
  folds <- sample(rep(seq_len(nfolds), length.out = nrow(data)))
  oof_pred <- double(nrow(data))
  
  for (k in seq_len(nfolds)) {
    fold_args <- ranger_args
    fold_args$data <- data[folds != k, ]
    fit_k <- do.call(
      ranger::ranger,
      fold_args
    )
    oof_pred[folds == k] <- predict(
      fit_k,
      data = data[folds == k, ]
    )$predictions
  }
  
  calib_df <- tibble(
    y = data[[y_name]],
    pred = oof_pred
  )
  
  calib_model <- lm(y ~ pred, data = calib_df)
  
  base_model <- do.call(ranger::ranger, ranger_args)
  
  structure(
    list(
      base_model = base_model,
      calib_model = calib_model
    ),
    class = "calibrated_ranger"
  )
}

predict.calibrated_ranger <- function(object, newdata, ...) {
  pred_raw <- predict(object$base_model, data = newdata)$predictions
  pred_calib <- as.numeric(cbind(1, pred_raw) %*% coef(object$calib_model))
  pred_calib
}

predict.uncalibrated_ranger <- function(object, newdata, ...) {
  predict(object$base_model, data = newdata)$predictions
}

# Calibrated probability random forest for binary outcomes.  The calibration
# model is fitted on inner-CV predictions, mirroring fit_calibrated_ranger().
fit_calibrated_ranger_classif <- function(formula, data, nfolds = 10,
                                          num.trees = 500,
                                          calibrate = TRUE) {
  y_name <- all.vars(formula)[1]
  train_data <- data
  train_data[[y_name]] <- factor(train_data[[y_name]])
  if (nlevels(train_data[[y_name]]) != 2) {
    stop("Random-forest classification requires a binary outcome: ", y_name)
  }

  positive_class <- levels(train_data[[y_name]])[[2]]

  if (!calibrate) {
    base_model <- ranger::ranger(
      formula, data = train_data, probability = TRUE,
      num.trees = num.trees, num.threads = 4
    )
    return(structure(
      list(base_model = base_model, positive_class = positive_class),
      class = "uncalibrated_ranger_classif"
    ))
  }
  folds <- sample(rep(seq_len(nfolds), length.out = nrow(train_data)))
  oof_prob <- double(nrow(train_data))

  for (k in seq_len(nfolds)) {
    fit_k <- ranger::ranger(
      formula, data = train_data[folds != k, ], probability = TRUE,
      num.trees = num.trees, num.threads = 4
    )
    oof_prob[folds == k] <- predict(
      fit_k, data = train_data[folds == k, ]
    )$predictions[, positive_class]
  }

  calibration_data <- tibble(
    y = as.integer(train_data[[y_name]] == positive_class),
    logit_prob = qlogis(clip_prob(oof_prob))
  )
  calibration_model <- glm(y ~ logit_prob, data = calibration_data,
                           family = binomial())
  base_model <- ranger::ranger(
    formula, data = train_data, probability = TRUE,
    num.trees = num.trees, num.threads = 4
  )

  structure(
    list(
      base_model = base_model,
      calibration_model = calibration_model,
      positive_class = positive_class
    ),
    class = "calibrated_ranger_classif"
  )
}

predict.calibrated_ranger_classif <- function(object, newdata, ...) {
  raw_prob <- predict(object$base_model, data = newdata)$predictions[
    , object$positive_class
  ]
  logit_prob <- qlogis(clip_prob(raw_prob))
  plogis(as.numeric(cbind(1, logit_prob) %*% coef(object$calibration_model)))
}

predict.uncalibrated_ranger_classif <- function(object, newdata, ...) {
  predict(object$base_model, data = newdata)$predictions[, object$positive_class]
}

