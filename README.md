# Replication Materials for "Mitigating Label Bias with Interpretable Rubric Embeddings"

Run the plotting script from this directory:

```sh
Rscript --vanilla clean_analysis.R
```

It reads the aggregate results in `replication_data/` and
writes the paper figures to `plots/`, including the appendix figures.

## Software

The plotting script requires R 4.4+ and the R packages `fs`,
`here`, `glue`, `ggplot2`, `patchwork`, `cowplot`, `scales`, and `tidyverse`.

The included model and random-forest scripts additionally use `glmnet`,
`glmnetUtils`, `pROC`, `ranger`, `RhpcBLASctl`, and `vroom`; `MASS` ships with standard R installations. 

## Contents and data availability

`clean_analysis.R` is fully runnable. 
`models.R` documents the updated regularized-regression workflow. The random-forest scripts document the bias-specific tuning workflow (`R/tune_random_forest_by_bias.R`), its 36-specification tuning aggregation, the 720-task tuned fitting grid (`R/run_random_forest_bias_tuned_grid.R`), AUC analysis, and aggregation.Those estimation scripts require confidential source records, extracted features, and embeddings, so they do not run as provided.

The released `.rds` files contain only pre-computed aggregate model outputs, summary statistics, binned score frequencies, and feature-level paired-change statistics. 
