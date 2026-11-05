library(glue)
library(fs)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)
library(tidyverse)

here::i_am("clean_analysis.R")

# Data --------------------------------------------------------------------

df_summary    <- read_rds(path("replication_data", "df_summary.rds"))
df_score_dist <- read_rds(path("replication_data", "df_score_dist.rds"))

res <- read_rds(path("replication_data", "modeling.rds")) %>%
  rename(feats = feats_fit, bias = b_fit, variant = variant_fit) %>% 
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    mean_gpa_est = mean(mean_gpa),
    mean_gpa_std.err = sd(mean_gpa) / sqrt(n()),
    mean_gre_quant_est = mean(mean_gre_quant),
    mean_gre_quant_std.err = sd(mean_gre_quant) / sqrt(n()),
    mean_gre_verbal_est = mean(mean_gre_verbal),
    mean_gre_verbal_std.err = sd(mean_gre_verbal) / sqrt(n()),
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) + (1 + 1 / n()) *
        var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
    rmse_biased_std.err = sqrt(
      mean(rmse_biased_std.err^2) + (1 + 1 / n()) *
        var(rmse_biased_est)
    ),
    rmse_biased_est = mean(rmse_biased_est),
    rmse_unbiased_std.err = sqrt(
      mean(rmse_unbiased_std.err^2) + (1 + 1 / n()) *
        var(rmse_unbiased_est)
    ),
    rmse_unbiased_est = mean(rmse_unbiased_est)
  ) %>%
  ungroup()

res_margin <- read_rds(path("replication_data", "marginalization.rds")) %>%
  rename(feats = feats_fit, bias = b_fit, variant = variant_fit) %>% 
  mutate(spec = "marginalized_gender") %>%
  filter(variant == "orig") %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    mean_gpa_est = mean(mean_gpa),
    mean_gpa_std.err = sd(mean_gpa) / sqrt(n()),
    mean_gre_quant_est = mean(mean_gre_quant),
    mean_gre_quant_std.err = sd(mean_gre_quant) / sqrt(n()),
    mean_gre_verbal_est = mean(mean_gre_verbal),
    mean_gre_verbal_std.err = sd(mean_gre_verbal) / sqrt(n()),
    rmse_biased_std.err = sqrt(
      mean(rmse_biased_std.err^2) + (1 + 1 / n()) *
        var(rmse_biased_est)
    ),
    rmse_biased_est = mean(rmse_biased_est),
    rmse_unbiased_std.err = sqrt(
      mean(rmse_unbiased_std.err^2) + (1 + 1 / n()) *
        var(rmse_unbiased_est)
    ),
    rmse_unbiased_est = mean(rmse_unbiased_est),
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) + (1 + 1 / n()) *
        var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
  ) %>%
  ungroup()

res_orth <- read_rds(path("replication_data", "orthogonal.rds")) %>%
  rename(feats = feats_fit, variant = variant_fit, bias = b_fit) %>% 
  mutate(spec = "orthogonalized") %>%
  filter(variant == "orig", feats == "feat_gp") %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    mean_gpa_est = mean(mean_gpa),
    mean_gpa_std.err = sd(mean_gpa) / sqrt(n()),
    mean_gre_quant_est = mean(mean_gre_quant),
    mean_gre_quant_std.err = sd(mean_gre_quant) / sqrt(n()),
    mean_gre_verbal_est = mean(mean_gre_verbal),
    mean_gre_verbal_std.err = sd(mean_gre_verbal) / sqrt(n()),
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) + (1 + 1 / n()) *
        var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
  ) %>%
  ungroup()

res_gender <- read_rds(path("replication_data", "aucs.rds")) %>%
  rename(feats = feats_fit, variant = variant_fit) %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec) %>%
  summarize(
    auc_true_gender_est = mean(auc_true_gender),
    auc_true_gender_std.err = sqrt(
      mean(auc_true_gender_std.err^2) + (1 + 1 / n()) * var(auc_true_gender)
    ),
    auc_presented_gender_est = mean(auc_presented_gender),
    auc_presented_gender_std.err = sqrt(
      mean(auc_presented_gender_std.err^2) +
        (1 + 1 / n()) * var(auc_presented_gender)
    ),
    auc_true_race_est = mean(auc_true_race),
    auc_true_race_std.err = sqrt(
      mean(auc_true_race_std.err^2) + (1 + 1 / n()) * var(auc_true_race)
    ),
    gender_pred_diff_est = mean(gender_pred_diff),
    gender_pred_diff_std.err = sqrt(
      mean(gender_pred_diff_std.err^2) +
        (1 + 1 / n()) * var(gender_pred_diff)
    )
  )

res_gender_ridge <- res_gender

res_ks <- read_rds(path("replication_data", "ks_theory.rds")) %>%
  mutate(
    empirical_est = mse_ks_est - mse_rb_est,
    empirical_std.err = sqrt(mse_ks_std.err^2 + mse_rb_std.err^2),
    theory_est = resid_var_rb_est - resid_var_ks_est - 2 * cov_ks_est,
    theory_std.err = sqrt(
      resid_var_rb_std.err^2 + resid_var_ks_std.err^2 + 4 * cov_ks_std.err^2
      )
  ) 

# Plotting Config ---------------------------------------------------------

theme_set(theme_bw())
theme_update(
  axis.title.x     = element_text(size = 8),
  axis.title.y     = element_text(size = 8),
  axis.text        = element_text(size = 8),
  legend.text      = element_text(size = 8),
  panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.line        = element_blank(),
  strip.background = element_rect(fill = "white", color = NA),
  strip.text       = element_text(size = 8),
  legend.position  = "bottom",
  legend.margin    = margin(t = -8),
  plot.margin      = margin(1, 1, 1, 1)
)

panel_theme <- theme(
  legend.position = "none",
  plot.margin     = margin(2, 4, 0, 4),
  axis.title.x    = element_text(size = 8, margin = margin(t = 1))
)

update_geom_defaults("point", list(size = 1))

method_colors <- c(
  "feat_gp orig"        = "#F8766D",
  "feat_gp neut"        = "#B79F00",
  "feat_ks orig"        = "#00BA38",
  "feat_rb orig"        = "#00BFC4",
  "marginalized_gender" = "#619CFF",
  "orthogonalized"     = "#F564E3"
)

# Figure 1 ----------------------------------------------------------------

## Figure 1(a): Causal Penalty (Dense Only) -----------------------------------

p_audit <- res %>%
  filter(spec == "feat_gp orig") %>%
  left_join(res_gender, by = "spec") %>%
  ggplot(
    aes(x = bias, y = audit_bias_cf_calib_est, color = spec, fill = spec)
  ) +
  geom_abline(
    data = ~ filter(res, spec == "feat_gp orig", bias == 0) %>%
      left_join(res_gender, by = "spec") %>%
      distinct(spec, gender_pred_diff_est, audit_bias_cf_calib_est),
    aes(
      slope = 2 * gender_pred_diff_est,
      intercept = audit_bias_cf_calib_est,
      color = spec
    ),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est - 1 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est + 1 * audit_bias_cf_calib_std.err,
    ),
    color = "transparent",
    alpha = 0.33,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est - 2 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est + 2 * audit_bias_cf_calib_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = expression("" ~ widehat(bias)[T](h) ~ "")
  ) +
  theme(legend.position = "none") +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin,",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Figure 1(b): Pct female in admitted class  ---------------------------------

true_prop_female <- df_summary$true_prop_female

p_pct_women <- res %>%
  filter(spec == "feat_gp orig") %>%
  select(bias, mean_prop_male_est, mean_prop_male_std.err, spec) %>%
  mutate(
    pct_female = 100 * (1 - mean_prop_male_est)
  ) %>%
  ggplot(aes(x = bias, y = pct_female, color = spec, fill = spec)) +
  geom_hline(yintercept = true_prop_female, linetype = "solid") +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = pct_female - 100 * mean_prop_male_std.err,
      ymax = pct_female + 100 * mean_prop_male_std.err,
    ),
    color = "transparent",
    alpha = 0.33,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - 200 * mean_prop_male_std.err,
      ymax = pct_female + 200 * mean_prop_male_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Pct. female (admitted class)",
    color = NULL
  ) +
  theme(
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Figure 1(c): Mean score of admitted class -------------------------------

true_avg_score <- df_summary$true_avg_score

p_avg_score <- res %>%
  filter(spec == "feat_gp orig") %>%
  select(bias, mean_score_est, mean_score_std.err, spec) %>%
  ggplot(aes(x = bias, y = mean_score_est, color = spec, fill = spec)) +
  geom_line() +
  geom_hline(yintercept = true_avg_score, linetype = "solid") +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 1 * mean_score_std.err,
      ymax = mean_score_est + 1 * mean_score_std.err,
    ),
    color = "transparent",
    alpha = 0.33,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean score (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(limits = c(14, 17)) +
  theme(legend.position = "none") +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Marg.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Synthesis ---------------------------------------------------------------

fig1 <- (p_audit + p_pct_women +  p_avg_score) &
  panel_theme &
  theme(axis.title.x = element_text(size = 8, margin = margin(t = 10)))

ggsave(
  path("plots", "audit_dense.pdf"),
  fig1,
  width = 5.5,
  height = 2.05
)


# Fig 3 Panels ------------------------------------------------------------

## Figure 3(a): Causal Penalty (All) -----------------------------------

p_audit_full <- res %>%
  filter(grepl("orig", spec)) %>%
  left_join(filter(res_gender), by = "spec") %>%
  bind_rows(res_orth) %>% 
  mutate(
    spec = factor(
      spec,
      #' Drop marginalized and redacted but keep them in factor level to preserve
      #' line color
      levels = c(
        "feat_gp orig", 
        "feat_gp neut", 
        "feat_ks orig", 
        "feat_rb orig",
        "marginalized_gender", 
        "orthogonalized"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = bias,
      y = audit_bias_cf_calib_est,
      group = spec,
      colour = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = 0) +
  geom_abline(
    data = ~ filter(., bias == 0) %>%
      distinct(spec, gender_pred_diff_est, audit_bias_cf_calib_est),
    aes(
      slope = 2 * gender_pred_diff_est, 
      intercept = audit_bias_cf_calib_est, 
      color = spec
    ),
    linetype = "dashed"
  ) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est - 1 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est + 1 * audit_bias_cf_calib_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est - 2 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est + 2 * audit_bias_cf_calib_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = expression("" ~ widehat(bias)[T](h) ~ "")
  ) +
  theme(legend.position = "none") +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Figure 3(b): Pct female in admitted class  ----------------------------------

p_pct_women_full <- res %>%
  filter(grepl("orig", spec) | spec == "feat_gp neut") %>%
  bind_rows(res_margin) %>%
  bind_rows(res_orth) %>%
  mutate(
    spec = factor(
      spec,
      levels = c(
        "feat_gp orig", "feat_gp neut", "feat_ks orig", "feat_rb orig",
        "marginalized_gender", "orthogonalized"
      )
    )
  ) %>%
  select(bias, mean_prop_male_est, mean_prop_male_std.err, spec) %>%
  mutate(
    pct_female = 100 * (1 - mean_prop_male_est)
  ) %>%
  ggplot(aes(x = bias, y = pct_female, color = spec, fill = spec)) +
  geom_hline(yintercept = true_prop_female) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = pct_female - 100 * mean_prop_male_std.err,
      ymax = pct_female + 100 * mean_prop_male_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - 200 * mean_prop_male_std.err,
      ymax = pct_female + 200 * mean_prop_male_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Pct. female (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Figure 3(c): Mean score of admitted class -------------------------------

p_avg_score_full <- res %>%
  filter(grepl("orig", spec) | spec == "feat_gp neut") %>%
  bind_rows(res_margin) %>%
  bind_rows(res_orth) %>%
  mutate(
    spec = factor(
      spec,
      levels = c(
        "feat_gp orig", "feat_gp neut", "feat_ks orig", "feat_rb orig",
        "marginalized_gender", "orthogonalized"
      )
    )
  ) %>%
  select(bias, mean_score_est, mean_score_std.err, spec) %>%
  ggplot(aes(x = bias, y = mean_score_est, color = spec, fill = spec)) +
  geom_hline(yintercept = true_avg_score) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 1 * mean_score_std.err,
      ymax = mean_score_est + 1 * mean_score_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean score (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(limits = c(14, 17)) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Synthesis ---------------------------------------------------------------

legend_grob <- cowplot::get_plot_component(
  p_pct_women_full +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box",
  return_all = TRUE
)[[3]]

fig3_panels <- (p_audit_full + p_pct_women_full + p_avg_score_full) &
  panel_theme &
  theme(axis.title.x = element_text(size = 8, margin = margin(t = 10)))

fig3 <- cowplot::plot_grid(
  legend_grob,
  patchworkGrob(fig3_panels),
  ncol        = 1,
  rel_heights = c(.05, 0.95)
)

ggsave(
  path("plots", "audit_all.pdf"),
  fig3,
  width  = 5.5,
  height = 2 / .93 + .05,
  device = cairo_pdf
)

# Figure 4 ----------------------------------------------------------------

## Figure 4(a): RMSE Biased -----------------------------------------------

p_rmse_biased <- res %>%
  filter(grepl("orig", spec)) %>%
  select(spec, bias, rmse_biased_est, rmse_biased_std.err) %>%
  ggplot(
    aes(
      x = bias,
      y = rmse_biased_est,
      group = spec,
      colour = spec,
      fill = spec
    )
  ) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = rmse_biased_est - 1 * rmse_biased_std.err,
      ymax = rmse_biased_est + 1 * rmse_biased_std.err
    ),
    color = "transparent", alpha = 0.3, show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_biased_est - 2 * rmse_biased_std.err,
      ymax = rmse_biased_est + 2 * rmse_biased_std.err
    ),
    color = "transparent", alpha = 0.1, show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = expression("RMSE, biased label (" ~ Y * "'" ~ ")"), color = NULL
  ) +
  guides(color = guide_legend(position = "inside")) +
  theme(
    legend.position.inside = c(0.4, 0.7),
    legend.justification   = c(0.5, 0),
    legend.text            = element_text(size = 7),
    legend.key.size        = unit(0.4, "cm"),
    legend.margin          = margin(0, 0, 0, 0),
    legend.background      = element_blank(),
    legend.key             = element_blank()
  ) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Figure 4(b): RMSE Unbiased -----------------------------------------------

p_rmse_unbiased <- res %>%
  filter(grepl("orig", spec)) %>%
  select(spec, bias, rmse_unbiased_est, rmse_unbiased_std.err) %>%
  ggplot(
    aes(
      x = bias,
      y = rmse_unbiased_est,
      group = spec,
      colour = spec,
      fill = spec
    )
  ) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = rmse_unbiased_est - 1 * rmse_unbiased_std.err,
      ymax = rmse_unbiased_est + 1 * rmse_unbiased_std.err
    ),
    color = "transparent", alpha = 0.3, show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_unbiased_est - 2 * rmse_unbiased_std.err,
      ymax = rmse_unbiased_est + 2 * rmse_unbiased_std.err
    ),
    color = "transparent", alpha = 0.1, show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = expression("RMSE, unbiased label (" ~ Y ~ ")"), color = NULL
  ) +
  guides(color = guide_legend(position = "inside")) +
  theme(
    legend.position = "none"
  ) + scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

## Synthesis ---------------------------------------------------------------

fig4 <- cowplot::plot_grid(
  patchworkGrob(
    (p_rmse_biased & theme(legend.position = "inside")) +
      (p_rmse_unbiased & panel_theme)
  ),
  textGrob(
    expression("Male advantage (" ~ b ~ ")"),
    gp = gpar(fontsize = 8),
    hjust = 0.5
  ),
  ncol = 1,
  rel_heights = c(.93, 0.07)
)

ggsave(
  path("plots", "rmse.pdf"),
  fig4,
  width  = 5.5 * 2 / 3,
  height = 2,
  device = cairo_pdf
)


# Appendix ----------------------------------------------------------------

## Score Distribution ------------------------------------------------------

p_score_dist <- df_score_dist %>%
  ggplot(aes(x = score, weight = n_applicants)) +
  geom_histogram(color = "black", fill = "gray", binwidth = 1, boundary = 0.5) +
  scale_x_continuous(limits = c(3.5, 20.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Applicant Score",
    y = "Number of Applicants"
  ) +
  theme(
    plot.margin = margin(1, 6, 1, 1)
  )

ggsave(
  path("plots", "score_dist.pdf"),
  plot = p_score_dist,
  height = 4,
  width = 5.5 / 1.5,
  device = cairo_pdf
)

## Kitchen Sink Theory ------------------------------------------------------

p_theory <- res_ks %>%
  pivot_longer(
    cols = c(ends_with("_est"), ends_with("_std.err")),
    names_pattern = "(.+)_(est|std.err)",
    names_to = c("stat", ".value")
  ) %>%
  filter(stat %in% c("empirical", "theory"), bias < 4) %>%
  group_by(bias, stat) %>%
  summarize(
    std.err = sqrt(mean(std.err^2) + (1 + 1 / n()) * var(est)),
    est     = mean(est)
  ) %>%
  ggplot(aes(x = bias, y = est, color = stat, fill = stat)) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = est - 1 * std.err,
      ymax = est + 1 * std.err,
    ),
    color = "transparent",
    alpha = 0.33,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = est - 2 * std.err,
      ymax = est + 2 * std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  scale_color_discrete(labels = c(
    "empirical" = "Empirical",
    "theory"    = "Theoretical"
  )) +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Difference in MSE",
    color = NULL
  ) +
  guides(color = guide_legend(position = "inside")) +
  theme(
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  scale_fill_discrete(guide = "none")

ggsave(
  path("plots", "diff_in_mse.pdf"),
  p_theory,
  width  = 5.5 / 1.5,
  height = 4,
  device = cairo_pdf
)

## Class composition detailed ----------------------------------------------

true_avg_gpa <- df_summary$true_avg_gpa
true_avg_gre_quant <- df_summary$true_avg_gre_quant
true_avg_gre_verbal <- df_summary$true_avg_gre_verbal

### Figure: Mean GPA of admitted class --------------------------------------

p_avg_gpa_full <- res %>%
  filter(grepl("orig", spec) | spec == "feat_gp neut") %>%
  bind_rows(res_margin) %>%
  bind_rows(res_orth) %>%
  mutate(
    spec = factor(
      spec,
      levels = c(
        "feat_gp orig", "feat_gp neut", "feat_ks orig", "feat_rb orig",
        "marginalized_gender", "orthogonalized"
      )
    )
  ) %>%
  select(bias, mean_gpa_est, mean_gpa_std.err, spec) %>%
  ggplot(aes(x = bias, y = mean_gpa_est, color = spec, fill = spec)) +
  geom_hline(yintercept = true_avg_gpa) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = mean_gpa_est - 1 * mean_gpa_std.err,
      ymax = mean_gpa_est + 1 * mean_gpa_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_gpa_est - 2 * mean_gpa_std.err,
      ymax = mean_gpa_est + 2 * mean_gpa_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean GPA (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(limits = c(3.5, 3.75)) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

### Figure: Mean GRE Quant of admitted class --------------------------------

p_avg_gre_quant_full <- res %>%
  filter(grepl("orig", spec) | spec == "feat_gp neut") %>%
  bind_rows(res_margin) %>%
  bind_rows(res_orth) %>% 
  mutate(
    spec = factor(
      spec,
      levels = c(
        "feat_gp orig", "feat_gp neut", "feat_ks orig", "feat_rb orig",
        "marginalized_gender", "orthogonalized"
      )
    )
  ) %>%
  select(bias, mean_gre_quant_est, mean_gre_quant_std.err, spec) %>%
  ggplot(aes(x = bias, y = mean_gre_quant_est, color = spec, fill = spec)) +
  geom_hline(yintercept = true_avg_gre_quant) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = mean_gre_quant_est - 1 * mean_gre_quant_std.err,
      ymax = mean_gre_quant_est + 1 * mean_gre_quant_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_gre_quant_est - 2 * mean_gre_quant_std.err,
      ymax = mean_gre_quant_est + 2 * mean_gre_quant_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Mean GRE Quant. (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(limits = c(160, 164)) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

### Figure: Mean GRE Verbal of admitted class -------------------------------

p_avg_gre_verbal_full <- res %>%
  filter(grepl("orig", spec) | spec == "feat_gp neut") %>%
  bind_rows(res_margin) %>%
  bind_rows(res_orth) %>%
  mutate(
    spec = factor(
      spec,
      levels = c(
        "feat_gp orig", "feat_gp neut", "feat_ks orig", "feat_rb orig",
        "marginalized_gender", "orthogonalized"
      )
    )
  ) %>%
  select(bias, mean_gre_verbal_est, mean_gre_verbal_std.err, spec) %>%
  ggplot(aes(x = bias, y = mean_gre_verbal_est, color = spec, fill = spec)) +
  geom_hline(yintercept = true_avg_gre_verbal) +
  geom_line() +
  geom_point() +
  geom_ribbon(
    aes(
      ymin = mean_gre_verbal_est - 1 * mean_gre_verbal_std.err,
      ymax = mean_gre_verbal_est + 1 * mean_gre_verbal_std.err,
    ),
    color = "transparent",
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_gre_verbal_est - 2 * mean_gre_verbal_std.err,
      ymax = mean_gre_verbal_est + 2 * mean_gre_verbal_std.err,
    ),
    color = "transparent",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean GRE Verbal (admitted class)",
    color = NULL
  ) +
  scale_y_continuous(limits = c(160, 164)) +
  scale_color_manual(
    values = method_colors,
    labels = c(
      "feat_gp orig"        = "Black-box Embs.",
      "feat_gp neut"        = "Redacted Embs.",
      "feat_ks orig"        = "Kitchen Sink",
      "feat_rb orig"        = "Rubric Embs.",
      "marginalized_gender" = "Margin.",
      "orthogonalized"     = "Orth."
    )
  ) +
  scale_fill_manual(values = method_colors, guide = "none")

### Synthesis ---------------------------------------------------------------

legend_grob_cov <- cowplot::get_plot_component(
  p_avg_gre_verbal_full +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box",
  return_all = TRUE
)[[3]]

fig_cov_panels <- (p_avg_gpa_full + p_avg_gre_quant_full +
                     p_avg_gre_verbal_full) &
  panel_theme &
  theme(axis.title.x = element_text(size = 8, margin = margin(t = 10)))

fig_cov <- cowplot::plot_grid(
  legend_grob_cov,
  patchworkGrob(fig_cov_panels),
  ncol        = 1,
  rel_heights = c(.05, 0.95)
)

ggsave(
  path("plots", "audit_covariates.pdf"),
  fig_cov,
  width  = 5.5,
  height = 2 / .93 + .05,
  device = cairo_pdf
)

# Tables ------------------------------------------------------------------

res_gender %>%
  select(
    spec,
    auc_true_gender_est,
    auc_true_gender_std.err,
    gender_pred_diff_est,
    gender_pred_diff_std.err
  ) %>%
  print()

res_gender %>%
  select(
    spec,
    auc_true_race_est,
    auc_true_race_std.err
  ) %>%
  print()



# Robustness checks --------------------------------------------------------

method_labels <- c(
  "feat_gp orig"        = "Black-box Embs.",
  "feat_gp neut"        = "Redacted Embs.",
  "feat_ks orig"        = "Kitchen Sink",
  "feat_rb orig"        = "Rubric Embs.",
  "marginalized_gender" = "Margin.",
  "orthogonalized"     = "Orth."
)

method_levels <- names(method_colors)

dir_create(path("plots", "appendix"))

## Alternative model families ---------------------------------------------

res_models <- tibble()
res_models_margin <- tibble()
res_models_gender <- tibble()

for (model in c("lasso", "rf")) {
  
  res_model <- read_rds(
    path("replication_data", glue("modeling_{model}.rds"))
  ) %>%
    rename(
      feats   = feats_fit,
      variant = variant_fit,
      bias    = b_fit
    ) %>%
    filter(selection_rule == "original_prediction") %>%
    mutate(spec = paste(feats, variant)) %>%
    group_by(spec, bias) %>%
    summarize(
      mean_prop_male_est = mean(prop_male),
      mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
      
      mean_score_est = mean(mean_score),
      mean_score_std.err = sd(mean_score) / sqrt(n()),
      
      audit_bias_cf_calib_std.err = sqrt(
        mean(audit_bias_cf_calib_std.err^2) +
          (1 + 1 / n()) * var(audit_bias_cf_calib_est)
      ),
      audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
      
      rmse_biased_std.err = sqrt(
        mean(rmse_biased_std.err^2) +
          (1 + 1 / n()) * var(rmse_biased_est)
      ),
      rmse_biased_est = mean(rmse_biased_est),
      
      rmse_unbiased_std.err = sqrt(
        mean(rmse_unbiased_std.err^2) +
          (1 + 1 / n()) * var(rmse_unbiased_est)
      ),
      rmse_unbiased_est = mean(rmse_unbiased_est),
      
      .groups = "drop"
    ) %>%
    mutate(model = model)
  
  if (model == "rf") {
    
    margin_file <- "modeling_rf.rds"
    
    res_gender_model <- read_rds(
      path("replication_data", "aucs_rf.rds")
    ) %>%
      rename(
        feats   = feats_fit,
        variant = variant_fit
      ) %>%
      mutate(spec = paste(feats, variant)) %>%
      group_by(spec) %>%
      summarize(
        gender_pred_diff_est = mean(gender_pred_diff),
        .groups = "drop"
      )
    
  } else {
    
    margin_file <- "marginalization_lasso.rds"
    
    res_gender_model <- res_gender_ridge %>%
      select(spec, gender_pred_diff_est)
    
  }
  
  res_margin_model <- read_rds(
    path("replication_data", margin_file)
  ) %>%
    rename(
      feats   = feats_fit,
      variant = variant_fit,
      bias    = b_fit
    ) %>%
    filter(
      selection_rule == "gender_marginalized_prediction",
      variant == "orig"
    ) %>%
    mutate(spec = "marginalized_gender") %>%
    group_by(spec, bias) %>%
    summarize(
      mean_prop_male_est = mean(prop_male),
      mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
      
      mean_score_est = mean(mean_score),
      mean_score_std.err = sd(mean_score) / sqrt(n()),
      
      .groups = "drop"
    ) %>%
    mutate(model = model)
  
  res_models <- bind_rows(
    res_models,
    res_model
  )
  
  res_models_margin <- bind_rows(
    res_models_margin,
    res_margin_model
  )
  
  res_models_gender <- bind_rows(
    res_models_gender,
    res_gender_model %>% mutate(model = model)
  )
  
}


## Alternative embedding model --------------------------------------------

res_emb <- read_rds(
  path("replication_data", "modeling_mixedbread.rds")
) %>%
  rename(
    feats   = feats_fit,
    variant = variant_fit,
    bias    = b_fit
  ) %>%
  filter(selection_rule == "original_prediction") %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) +
        (1 + 1 / n()) * var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
    
    rmse_biased_std.err = sqrt(
      mean(rmse_biased_std.err^2) +
        (1 + 1 / n()) * var(rmse_biased_est)
    ),
    rmse_biased_est = mean(rmse_biased_est),
    
    rmse_unbiased_std.err = sqrt(
      mean(rmse_unbiased_std.err^2) +
        (1 + 1 / n()) * var(rmse_unbiased_est)
    ),
    rmse_unbiased_est = mean(rmse_unbiased_est),
    
    .groups = "drop"
  ) %>%
  mutate(model = "mixedbread")

res_margin_emb <- read_rds(
  path("replication_data", "marginalization_mixedbread.rds")
) %>%
  rename(
    feats   = feats_fit,
    variant = variant_fit,
    bias    = b_fit
  ) %>%
  filter(
    selection_rule == "gender_marginalized_prediction",
    variant == "orig"
  ) %>%
  mutate(spec = "marginalized_gender") %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    .groups = "drop"
  ) %>%
  mutate(model = "mixedbread")

res_gender_emb <- read_rds(
  path("replication_data", "aucs_mixedbread.rds")
) %>%
  rename(feats = feats_name) %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec) %>%
  summarize(
    auc_true_gender_est = mean(auc_true_gender),
    auc_true_gender_std.err = sqrt(
      mean(auc_true_gender_std.err^2) +
        (1 + 1 / n()) * var(auc_true_gender)
    ),
    
    auc_presented_gender_est = mean(auc_presented_gender),
    auc_presented_gender_std.err = sqrt(
      mean(auc_presented_gender_std.err^2) +
        (1 + 1 / n()) * var(auc_presented_gender)
    ),
    
    auc_true_race_est = mean(auc_true_race),
    auc_true_race_std.err = sqrt(
      mean(auc_true_race_std.err^2) +
        (1 + 1 / n()) * var(auc_true_race)
    ),
    
    gender_pred_diff_est = mean(gender_pred_diff),
    gender_pred_diff_std.err = sqrt(
      mean(gender_pred_diff_std.err^2) +
        (1 + 1 / n()) * var(gender_pred_diff)
    ),
    
    .groups = "drop"
  )

## Combined model and embedding results ------------------------------------

res_models <- bind_rows(
  res_models,
  res_emb
) %>%
  mutate(
    model = factor(
      model,
      levels = c("lasso", "rf", "mixedbread"),
      labels = c("LASSO", "Random forest", "Mixedbread")
    ),
    spec = factor(spec, levels = method_levels)
  )

res_models_margin <- bind_rows(
  res_models_margin,
  res_margin_emb
) %>%
  mutate(
    model = factor(
      model,
      levels = c("lasso", "rf", "mixedbread"),
      labels = c("LASSO", "Random forest", "Mixedbread")
    ),
    spec = factor(spec, levels = method_levels)
  )

res_models_gender <- bind_rows(
  res_models_gender,
  res_gender_emb %>%
    select(spec, gender_pred_diff_est) %>%
    mutate(model = "mixedbread")
) %>%
  mutate(
    model = factor(
      model,
      levels = c("lasso", "rf", "mixedbread"),
      labels = c("LASSO", "Random forest", "Mixedbread")
    ),
    spec = factor(spec, levels = method_levels)
  )

res_models_audit <- res_models %>%
  filter(spec %in% c(
    "feat_gp orig",
    "feat_ks orig",
    "feat_rb orig"
  )) %>%
  left_join(
    res_models_gender,
    by = c("model", "spec")
  )

res_models_selection <- res_models %>%
  filter(spec %in% c(
    "feat_gp orig",
    "feat_gp neut",
    "feat_ks orig",
    "feat_rb orig"
  )) %>%
  bind_rows(res_models_margin)


### Audit bias -------------------------------------------------------------

p_audit_models <- res_models_audit %>%
  ggplot(
    aes(
      x = bias,
      y = audit_bias_cf_calib_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = 0) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est -
        2 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est +
        2 * audit_bias_cf_calib_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est -
        audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est +
        audit_bias_cf_calib_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_abline(
    data = res_models_audit %>% filter(bias == 0),
    aes(
      slope = 2 * gender_pred_diff_est,
      intercept = audit_bias_cf_calib_est,
      color = spec
    ),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  facet_grid(model ~ ., switch = "y") +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = NULL,
    y = expression(widehat(bias)[T](h)),
    color = NULL
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  )


### Percentage female ------------------------------------------------------

p_pct_women_models <- res_models_selection %>%
  mutate(
    pct_female = 100 * (1 - mean_prop_male_est),
    pct_female_std.err = 100 * mean_prop_male_std.err
  ) %>%
  ggplot(
    aes(
      x = bias,
      y = pct_female,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_prop_female) +
  geom_ribbon(
    aes(
      ymin = pct_female - 2 * pct_female_std.err,
      ymax = pct_female + 2 * pct_female_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - pct_female_std.err,
      ymax = pct_female + pct_female_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  facet_grid(model ~ .) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Pct. female (admitted class)",
    color = NULL
  ) +
  theme(strip.text.y = element_blank())


### Mean selected score ----------------------------------------------------

p_avg_score_models <- res_models_selection %>%
  ggplot(
    aes(
      x = bias,
      y = mean_score_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_avg_score) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - mean_score_std.err,
      ymax = mean_score_est + mean_score_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  facet_grid(model ~ .) +
  coord_cartesian(ylim = c(14, 17)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Mean score (admitted class)",
    color = NULL
  ) +
  theme(strip.text.y = element_blank())


### Synthesis --------------------------------------------------------------

legend_grob_models <- cowplot::get_plot_component(
  p_pct_women_models +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

fig_models_panels <- (
  p_audit_models +
    p_pct_women_models +
    p_avg_score_models
) &
  panel_theme

fig_models <- cowplot::plot_grid(
  legend_grob_models,
  patchworkGrob(fig_models_panels),
  ncol        = 1,
  rel_heights = c(.07, .93)
)

ggsave(
  path("plots", "appendix", "robustness_models.pdf"),
  fig_models,
  width  = 5.5,
  height = 6.3,
  device = cairo_pdf
)

## Combined RMSE results ---------------------------------------------------

res_models_rmse <- res_models %>%
  filter(
    model %in% c("LASSO", "Random forest"),
    spec %in% c(
      "feat_gp orig",
      "feat_ks orig",
      "feat_rb orig"
    )
  ) %>%
  mutate(model = droplevels(model))

p_rmse_biased_models <- res_models_rmse %>%
  ggplot(
    aes(
      x = bias,
      y = rmse_biased_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_biased_est - 2 * rmse_biased_std.err,
      ymax = rmse_biased_est + 2 * rmse_biased_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_biased_est - rmse_biased_std.err,
      ymax = rmse_biased_est + rmse_biased_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  facet_grid(model ~ ., switch = "y") +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = expression("RMSE, biased label (" ~ Y * "'" ~ ")"),
    color = NULL
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  )

p_rmse_unbiased_models <- res_models_rmse %>%
  ggplot(
    aes(
      x = bias,
      y = rmse_unbiased_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_unbiased_est - 2 * rmse_unbiased_std.err,
      ymax = rmse_unbiased_est + 2 * rmse_unbiased_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = rmse_unbiased_est - rmse_unbiased_std.err,
      ymax = rmse_unbiased_est + rmse_unbiased_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  facet_grid(model ~ .) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = expression("RMSE, unbiased label (" ~ Y ~ ")"),
    color = NULL
  ) +
  theme(strip.text.y = element_blank())

legend_grob_rmse <- cowplot::get_plot_component(
  p_rmse_biased_models +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

fig_rmse_panels <- (
  p_rmse_biased_models +
    p_rmse_unbiased_models
) &
  panel_theme

fig_rmse_models <- cowplot::plot_grid(
  legend_grob_rmse,
  patchworkGrob(fig_rmse_panels),
  ncol        = 1,
  rel_heights = c(.07, .93)
)

ggsave(
  path("plots", "appendix", "robustness_rmse.pdf"),
  fig_rmse_models,
  width  = 5.5,
  height = 4.3,
  device = cairo_pdf
)


# Alternative attributes and outcomes --------------------------------------

## Race -------------------------------------------------------------------

res_race <- read_rds(
  path("replication_data", "modeling_race.rds")
) %>%
  rename(
    feats   = feats_fit,
    variant = variant_fit,
    bias    = b_fit
  ) %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_nonwhite_est = mean(prop_nonwhite),
    mean_prop_nonwhite_std.err = sd(prop_nonwhite) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    .groups = "drop"
  ) %>%
  filter(spec %in% c(
    "feat_gp orig",
    "feat_ks orig",
    "feat_rb orig"
  )) %>%
  mutate(spec = factor(spec, levels = method_levels))

true_prop_nonwhite <- df_summary$true_prop_nonwhite
true_avg_score_domestic <- df_summary$true_avg_score_domestic

p_pct_nonwhite <- res_race %>%
  mutate(
    pct_nonwhite = 100 * mean_prop_nonwhite_est,
    pct_nonwhite_std.err = 100 * mean_prop_nonwhite_std.err
  ) %>%
  ggplot(
    aes(
      x = bias,
      y = pct_nonwhite,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_prop_nonwhite) +
  geom_ribbon(
    aes(
      ymin = pct_nonwhite - 2 * pct_nonwhite_std.err,
      ymax = pct_nonwhite + 2 * pct_nonwhite_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_nonwhite - pct_nonwhite_std.err,
      ymax = pct_nonwhite + pct_nonwhite_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    title = "Race",
    x = expression("White advantage (" ~ b ~ ")"),
    y = "Pct. nonwhite (selected class)",
    color = NULL
  )

p_avg_score_race <- res_race %>%
  ggplot(
    aes(
      x = bias,
      y = mean_score_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_avg_score_domestic) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - mean_score_std.err,
      ymax = mean_score_est + mean_score_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(14, 17)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    title = " ",
    x = expression("White advantage (" ~ b ~ ")"),
    y = "Mean score (selected class)",
    color = NULL
  )


## Observed first-year grades ---------------------------------------------

res_2024 <- read_rds(
  path("replication_data", "modeling_first_year.rds")
) %>%
  rename(
    feats   = feats_fit,
    bias    = b_fit,
    variant = variant_fit
  ) %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    .groups = "drop"
  ) %>%
  filter(spec %in% c(
    "feat_gp orig",
    "feat_ks orig",
    "feat_rb orig"
  )) %>%
  mutate(spec = factor(spec, levels = method_levels))

df_summary_2024 <- read_rds(
  path("replication_data", "df_summary_first_year.rds")
)

true_prop_female_2024 <- df_summary_2024$true_prop_female
true_avg_score_2024 <- df_summary_2024$true_avg_score

p_pct_women_2024 <- res_2024 %>%
  mutate(
    pct_female = 100 * (1 - mean_prop_male_est),
    pct_female_std.err = 100 * mean_prop_male_std.err
  ) %>%
  ggplot(
    aes(
      x = bias,
      y = pct_female,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = 100 * true_prop_female_2024) +
  geom_ribbon(
    aes(
      ymin = pct_female - 2 * pct_female_std.err,
      ymax = pct_female + 2 * pct_female_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - pct_female_std.err,
      ymax = pct_female + pct_female_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    title = "First-year quantitative GPA",
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Pct. female (selected class)",
    color = NULL
  )

p_avg_score_2024 <- res_2024 %>%
  ggplot(
    aes(
      x = bias,
      y = mean_score_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_avg_score_2024) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - mean_score_std.err,
      ymax = mean_score_est + mean_score_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(3.5, 4)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    title = " ",
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Mean quantitative GPA (selected class)",
    color = NULL
  )


## Synthesis --------------------------------------------------------------

legend_grob_outcomes <- cowplot::get_plot_component(
  p_pct_nonwhite +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

p_pct_nonwhite <- p_pct_nonwhite + labs(title = NULL)
p_avg_score_race <- p_avg_score_race + labs(title = NULL)
p_pct_women_2024 <- p_pct_women_2024 + labs(title = NULL)
p_avg_score_2024 <- p_avg_score_2024 + labs(title = NULL)

fig_outcomes_panels <- (
  wrap_elements(
    full = textGrob(
      "Race",
      gp = gpar(fontsize = 8)
    )
  ) +
    p_pct_nonwhite +
    p_avg_score_race +
    wrap_elements(
      full = textGrob(
        "First-year\nquantitative GPA",
        gp = gpar(fontsize = 8)
      )
    ) +
    p_pct_women_2024 +
    p_avg_score_2024 +
    plot_layout(
      ncol = 3,
      widths = c(.4, 1, 1)
    )
) &
  panel_theme

fig_outcomes <- cowplot::plot_grid(
  legend_grob_outcomes,
  patchworkGrob(fig_outcomes_panels),
  ncol        = 1,
  rel_heights = c(.07, .93)
)

ggsave(
  path("plots", "appendix", "robustness_outcomes.pdf"),
  fig_outcomes,
  width  = 5.5,
  height = 4.5,
  device = cairo_pdf
)

# Major-based advantage ----------------------------------------------------

res_correlated <- read_rds(
  path("replication_data", "modeling_attribute_bias.rds")
) %>%
  mutate(spec = paste(feats, variant)) %>%
  group_by(spec, bias) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) +
        (1 + 1 / n()) * var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
    
    .groups = "drop"
  ) %>%
  filter(spec %in% c(
    "feat_gp orig",
    "feat_ks orig",
    "feat_rb orig"
  )) %>%
  mutate(spec = factor(spec, levels = method_levels))

p_audit_correlated <- res_correlated %>%
  ggplot(
    aes(
      x = bias,
      y = audit_bias_cf_calib_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = 0) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est -
        2 * audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est +
        2 * audit_bias_cf_calib_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = audit_bias_cf_calib_est -
        audit_bias_cf_calib_std.err,
      ymax = audit_bias_cf_calib_est +
        audit_bias_cf_calib_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = NULL,
    y = expression(widehat(bias)[T](h)),
    color = NULL
  )

p_pct_women_correlated <- res_correlated %>%
  mutate(
    pct_female = 100 * (1 - mean_prop_male_est),
    pct_female_std.err = 100 * mean_prop_male_std.err
  ) %>%
  ggplot(
    aes(
      x = bias,
      y = pct_female,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_prop_female) +
  geom_ribbon(
    aes(
      ymin = pct_female - 2 * pct_female_std.err,
      ymax = pct_female + 2 * pct_female_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - pct_female_std.err,
      ymax = pct_female + pct_female_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = expression("Major-based advantage (" ~ b ~ ")"),
    y = "Pct. female (admitted class)",
    color = NULL
  )

p_avg_score_correlated <- res_correlated %>%
  ggplot(
    aes(
      x = bias,
      y = mean_score_est,
      color = spec,
      fill = spec
    )
  ) +
  geom_hline(yintercept = true_avg_score) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err
    ),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - mean_score_std.err,
      ymax = mean_score_est + mean_score_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(14, 17)) +
  scale_color_manual(
    values = method_colors,
    labels = method_labels
  ) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Mean score (admitted class)",
    color = NULL
  )

legend_grob_correlated <- cowplot::get_plot_component(
  p_pct_women_correlated +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 7),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

fig_correlated_panels <- (
  p_audit_correlated +
    p_pct_women_correlated +
    p_avg_score_correlated
) &
  panel_theme

fig_correlated <- cowplot::plot_grid(
  legend_grob_correlated,
  patchworkGrob(fig_correlated_panels),
  ncol        = 1,
  rel_heights = c(.08, .92)
)

ggsave(
  path("plots", "appendix", "audit_all_correlated.pdf"),
  fig_correlated,
  width  = 5.5,
  height = 2.4,
  device = cairo_pdf
)

# Rubric feature-count sweep -----------------------------------------------

res_ksweep <- read_rds(
  path("replication_data", "modeling_feature_sweep.rds")
) %>%
  group_by(k_fit, b_fit) %>%
  summarize(
    mean_prop_male_est = mean(prop_male),
    mean_prop_male_std.err = sd(prop_male) / sqrt(n()),
    
    mean_score_est = mean(mean_score),
    mean_score_std.err = sd(mean_score) / sqrt(n()),
    
    audit_bias_cf_calib_std.err = sqrt(
      mean(audit_bias_cf_calib_std.err^2) +
        (1 + 1 / n()) * var(audit_bias_cf_calib_est)
    ),
    audit_bias_cf_calib_est = mean(audit_bias_cf_calib_est),
    
    mean_gpa_est = mean(mean_gpa),
    mean_gpa_std.err = sd(mean_gpa) / sqrt(n()),
    
    mean_gre_quant_est = mean(mean_gre_quant),
    mean_gre_quant_std.err = sd(mean_gre_quant) / sqrt(n()),
    
    mean_gre_verbal_est = mean(mean_gre_verbal),
    mean_gre_verbal_std.err = sd(mean_gre_verbal) / sqrt(n()),
    
    .groups = "drop"
  ) %>%
  rename(bias = b_fit)

res_ksweep_plot <- res_ksweep %>%
  mutate(
    k_fit = factor(k_fit, levels = sort(unique(k_fit))),
    pct_female = 100 * (1 - mean_prop_male_est),
    pct_female_std.err = 100 * mean_prop_male_std.err
  )


## Audit bias -------------------------------------------------------------

p_audit_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = audit_bias_cf_calib_est,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = 0) +
  # geom_ribbon(
  #   aes(
  #     ymin = audit_bias_cf_calib_est -
  #       2 * audit_bias_cf_calib_std.err,
  #     ymax = audit_bias_cf_calib_est +
  #       2 * audit_bias_cf_calib_std.err
  #   ),
  #   alpha = 0.05,
  #   color = NA,
  #   show.legend = FALSE
  # ) +
  geom_line() +
  geom_point() +
  scale_fill_discrete(guide = "none") +
  labs(
    x = NULL,
    y = expression(widehat(bias)[T](h)),
    color = "k"
  )


## Gender composition -----------------------------------------------------

p_pct_female_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = pct_female,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = true_prop_female) +
  geom_ribbon(
    aes(
      ymin = pct_female - 2 * pct_female_std.err,
      ymax = pct_female + 2 * pct_female_std.err
    ),
    alpha = 0.05,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = pct_female - 1 * pct_female_std.err,
      ymax = pct_female + 1 * pct_female_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1)
  ) +
  coord_cartesian(ylim = c(0, 65)) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Pct. female (admitted class)",
    color = "k"
  )


## Mean selected score ----------------------------------------------------

p_quality_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = mean_score_est,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = true_avg_score) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 2 * mean_score_std.err,
      ymax = mean_score_est + 2 * mean_score_std.err
    ),
    alpha = 0.05,
    color = NA,
    show.legend = FALSE
  ) +
  geom_ribbon(
    aes(
      ymin = mean_score_est - 1 * mean_score_std.err,
      ymax = mean_score_est + 1 * mean_score_std.err
    ),
    alpha = 0.3,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(14, 17)) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = NULL,
    y = "Mean score (admitted class)",
    color = "k"
  )


## Synthesis --------------------------------------------------------------

legend_grob_ksweep <- cowplot::get_plot_component(
  p_pct_female_ksweep +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

fig_ksweep_panels <- (
  p_audit_ksweep +
    p_pct_female_ksweep +
    p_quality_ksweep
) &
  panel_theme

fig_ksweep <- cowplot::plot_grid(
  legend_grob_ksweep,
  patchworkGrob(fig_ksweep_panels),
  ncol        = 1,
  rel_heights = c(.08, .92)
)



## Other measures of cohort quality ---------------------------------------

p_avg_gpa_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = mean_gpa_est,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = true_avg_gpa) +
  geom_ribbon(
    aes(
      ymin = mean_gpa_est - 2 * mean_gpa_std.err,
      ymax = mean_gpa_est + 2 * mean_gpa_std.err
    ),
    alpha = 0.05,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(3.5, 3.75)) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = NULL,
    y = "Mean GPA (admitted class)",
    color = "k"
  )

p_avg_gre_quant_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = mean_gre_quant_est,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = true_avg_gre_quant) +
  geom_ribbon(
    aes(
      ymin = mean_gre_quant_est - 2 * mean_gre_quant_std.err,
      ymax = mean_gre_quant_est + 2 * mean_gre_quant_std.err
    ),
    alpha = 0.05,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(161, 164)) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = expression("Male advantage (" ~ b ~ ")"),
    y = "Mean GRE Quant. (admitted class)",
    color = "k"
  )

p_avg_gre_verbal_ksweep <- res_ksweep_plot %>%
  ggplot(
    aes(
      x = bias,
      y = mean_gre_verbal_est,
      color = k_fit,
      fill = k_fit
    )
  ) +
  geom_hline(yintercept = true_avg_gre_verbal) +
  geom_ribbon(
    aes(
      ymin = mean_gre_verbal_est - 2 * mean_gre_verbal_std.err,
      ymax = mean_gre_verbal_est + 2 * mean_gre_verbal_std.err
    ),
    alpha = 0.05,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(160, 164)) +
  scale_fill_discrete(guide = "none") +
  labs(
    x = NULL,
    y = "Mean GRE Verbal (admitted class)",
    color = "k"
  )

## Synthesis --------------------------------------------------------------

legend_grob_ksweep <- cowplot::get_plot_component(
  p_pct_female_ksweep +
    guides(color = guide_legend(nrow = 1)) +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.margin   = margin(0, 0, 0, 0),
      legend.key      = element_blank()
    ),
  pattern = "guide-box-bottom"
)

fig_ksweep_panels <- (
  p_audit_ksweep +
    p_pct_female_ksweep +
    p_quality_ksweep +
    p_avg_gpa_ksweep +
    p_avg_gre_quant_ksweep +
    p_avg_gre_verbal_ksweep +
    plot_layout(ncol = 3)
) &
  panel_theme

fig_ksweep <- cowplot::plot_grid(
  legend_grob_ksweep,
  patchworkGrob(fig_ksweep_panels),
  ncol        = 1,
  rel_heights = c(.05, .95)
)

ggsave(
  path("plots", "appendix", "ksweep.pdf"),
  fig_ksweep,
  width  = 5.5,
  height = 4.5,
  device = cairo_pdf
)

# Feature stability --------------------------------------------------------

# This release contains only feature-level paired-change summaries.  It does
# not include applicant records, counterfactual records, or fitted models.
stab_rb <- read_rds(path("replication_data", "feature_stability.rds")) %>%
  mutate(
    source = factor(source, levels = c("Essay", "Resume", "Letters")),
    feature_id = factor(feature_id, levels = rev(feature_id))
  )

p_shift <- stab_rb %>%
  ggplot(aes(x = mean_delta, y = feature_id)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_linerange(
    aes(
      xmin = mean_delta - 2 * se_delta,
      xmax = mean_delta + 2 * se_delta
    ),
    orientation = "y",
    linewidth = 0.3
  ) +
  geom_linerange(
    aes(
      xmin = mean_delta - se_delta,
      xmax = mean_delta + se_delta
    ),
    orientation = "y",
    linewidth = 0.9
  ) +
  geom_point(size = 1) +
  facet_grid(
    source ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    x = "Male minus female difference",
    y = "Rubric feature"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(size = 8),
    axis.title.y     = element_text(size = 8),
  )

ggsave(
  path("plots", "appendix", "feature_shift.pdf"),
  p_shift,
  height = max(4, 0.11 * nrow(stab_rb) + 1),
  width  = 5.5,
  device = cairo_pdf
)

# Inline feature-stability statistics --------------------------------------

stab_rb %>%
  summarize(
    n_features = n(),
    n_changed = sum(change_rate > 0, na.rm = TRUE),
    n_significant = sum(significant, na.rm = TRUE)
  ) %>%
  print()

# Appendix: Tables ------------------------------------------------------------------

## Gender AUC --------------------------------------------------------------

res_gender_emb %>%
  select(
    spec,
    auc_true_gender_est,
    auc_true_gender_std.err,
    gender_pred_diff_est,
    gender_pred_diff_std.err
  ) %>%
  print()

# Adjusted gender score gap ------------------------------------------------

read_rds(path("replication_data", "baseline_model_stats.rds")) %>%
  print()
