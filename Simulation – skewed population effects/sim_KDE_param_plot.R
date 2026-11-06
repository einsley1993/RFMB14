################################################################################
# Simulation results: two separate plots
# Plot 1: Shared metrics (Raw difference, MAE, RMSE, Mean P_orig, Overestimation)
# Plot 2: Method-specific metrics (Coverage, PI Width), KDE vs Parametric, side by side
# Columns: 4 parameters (eta, tau2_c, beta_mod, k)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

scenario    <- "skewed"   # "skewed" (main design) or "normal" (supplement)
file_suffix <- if (scenario == "skewed") "_skewed" else ""

kde_raw <- read.csv(paste0("sim_results_KDE", file_suffix, ".csv"),
                    stringsAsFactors = FALSE)
par_raw <- read.csv(paste0("sim_results_param", file_suffix, ".csv"),
                    stringsAsFactors = FALSE)


if ("dist" %in% names(kde_raw)) kde_raw <- filter(kde_raw, dist == "normal")
if ("dist" %in% names(par_raw)) par_raw <- filter(par_raw, dist == "normal")

model_levels <- c("dmaMod", "mrMod2", "metaCorr", "mrMod3")
model_labels <- c("MA-U", "MR-U", "MA-C", "MR-C")
model_colors <- c("MA-U" = "#D55E00", "MR-U" = "#E69F00",
                  "MA-C" = "#0072B2", "MR-C" = "#009E73")
model_shapes <- c("MA-U" = 16, "MR-U" = 17,
                  "MA-C" = 15, "MR-C" = 18)

all_metrics <- c("Raw difference", "MAE", "RMSE",
                 "Mean P_orig (%)", "Overestimation (%)",
                 "Coverage (%)", "PI Width")

param_levels <- c(
  "eta ~ '(publication bias)'",
  "tau[c]^2 ~ '(heterogeneity)'",
  "beta[mod] ~ '(moderator strength)'",
  "k ~ '(no. MA studies)'"
)

# Marginal summaries 

make_marginal <- function(df, x_var, param_label) {
  df |>
    group_by(across(all_of(c(x_var, "model")))) |>
    summarise(
      `Raw difference`     = mean(raw_diff_mean,    na.rm = TRUE),
      MAE                  = mean(mae_mean,         na.rm = TRUE),
      RMSE                 = mean(rmse_mean,        na.rm = TRUE),
      `Mean P_orig (%)`    = mean(mean_p_orig_mean, na.rm = TRUE) * 100,
      `Overestimation (%)` = mean(pct_overest_mean, na.rm = TRUE) * 100,
      `Coverage (%)`       = mean(coverage_mean,    na.rm = TRUE) * 100,
      `PI Width`           = mean(pi_width_mean,    na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      Model     = factor(model, levels = model_levels,
                         labels = model_labels),
      x_value   = .data[[x_var]],
      parameter = param_label
    ) |>
    select(-model, -all_of(x_var))
}

build_marginals <- function(df) {
  bind_rows(
    make_marginal(df, "eta",      param_levels[1]),
    make_marginal(df, "tau2_c",   param_levels[2]),
    make_marginal(df, "beta_mod", param_levels[3]),
    make_marginal(df, "k",        param_levels[4])
  )
}

kde_marg <- build_marginals(kde_raw)
par_marg <- build_marginals(par_raw)

# Theme 

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position   = "bottom",
    legend.key.width  = unit(1.5, "cm"),
    strip.text.x      = element_text(size = 9,  face = "bold"),
    strip.text.y      = element_text(size = 8,  face = "bold"),
    panel.grid.minor  = element_blank(),
    plot.caption      = element_text(size = 7, color = "grey40",
                                     hjust = 0),
    axis.text.x       = element_text(size = 7),
    axis.text.y       = element_text(size = 7)
  )

# Plot 1: Shared metrics 

shared_levels <- c("Raw difference", "MAE", "RMSE",
                   "Mean P_orig (%)", "Overestimation (%)")

shared_data <- kde_marg |>
  pivot_longer(cols = all_of(all_metrics),
               names_to = "Metric", values_to = "Value") |>
  filter(Metric %in% shared_levels) |>
  mutate(
    Metric    = factor(Metric, levels = shared_levels),
    parameter = factor(parameter, levels = param_levels)
  )

ref_shared <- tribble(
  ~Metric,              ~ref,
  "Raw difference",        0,
  "Overestimation (%)",   50
) |>
  mutate(Metric = factor(Metric, levels = shared_levels)) |>
  crossing(parameter = param_levels) |>
  mutate(parameter = factor(parameter, levels = param_levels))

p1 <- ggplot(shared_data,
             aes(x = x_value, y = Value,
                 color = Model, shape = Model, group = Model)) +
  geom_hline(data = ref_shared, aes(yintercept = ref),
             linetype = "solid", color = "grey60", linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  facet_grid(
    rows     = vars(Metric),
    cols     = vars(parameter),
    scales   = "free",
    labeller = labeller(Metric = label_value,
                        parameter = label_parsed)
  ) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  scale_x_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL, color = "Model", shape = "Model",
       caption = paste(
         "Each panel: marginal effect of one parameter,",
         "averaged over all remaining parameters.",
         "Grey lines: nominal targets",
         "(Raw difference = 0; Overestimation = 50%)."
       )) +
  base_theme

ggsave(paste0("sim_metrics1", file_suffix, ".jpeg"), p1, width = 8, height = 10)


# Plot 2: Method-specific metrics 

method_levels <- c(
  "Coverage (%) [KDE]", "Coverage (%) [Parametric]",
  "PI Width [KDE]", "PI Width [Parametric]"
)

make_method_long <- function(df, suffix) {
  df |>
    pivot_longer(cols = all_of(all_metrics),
                 names_to = "Metric", values_to = "Value") |>
    filter(Metric %in% c("Coverage (%)", "PI Width")) |>
    mutate(Metric = paste0(Metric, " [", suffix, "]"))
}

method_data <- bind_rows(
  make_method_long(kde_marg, "KDE"),
  make_method_long(par_marg, "Parametric")
) |>
  mutate(
    Metric    = factor(Metric, levels = method_levels),
    parameter = factor(parameter, levels = param_levels)
  )

ref_method <- tribble(
  ~Metric,                     ~ref,
  "Coverage (%) [KDE]",         95,
  "Coverage (%) [Parametric]",  95
) |>
  mutate(Metric = factor(Metric, levels = method_levels)) |>
  crossing(parameter = param_levels) |>
  mutate(parameter = factor(parameter, levels = param_levels))

p2 <- ggplot(method_data,
             aes(x = x_value, y = Value,
                 color = Model, shape = Model, group = Model)) +
  geom_hline(data = ref_method, aes(yintercept = ref),
             linetype = "solid", color = "grey60", linewidth = 0.4,
             inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  facet_grid(
    rows     = vars(Metric),
    cols     = vars(parameter),
    scales   = "free",
    labeller = labeller(Metric = label_value,
                        parameter = label_parsed)
  ) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  scale_x_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL, color = "Model", shape = "Model",
       caption = paste(
         "Each panel: marginal effect of one parameter,",
         "averaged over all remaining parameters.",
         "Grey line: nominal target (Coverage = 95%).",
         "PI Width: MA-U and MA-C overlap."
       )) +
  base_theme

ggsave(paste0("sim_metrics2", file_suffix, ".jpeg"), p2, width = 8, height = 8)
