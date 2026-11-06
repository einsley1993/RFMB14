################################################################################
# Simulation Study: 4-model comparison (skewed population effects)
#
#
# Models:
#   dmaMod   – intercept-only naive MA (no moderator, no pub-bias correction)
#   mrMod2   – meta-regression with z1, z2, z3, uncorrected
#   metaCorr – intercept-only IPPW-corrected MA (no moderator)
#   mrMod3   – meta-regression with z1, z2, z3, IPPW-corrected (primary model)
#
#
# Evaluation metrics:
#   1. Raw difference (REP pooled mean - MA predicted mean)
#   2. Mean P_orig
#   3. MAE
#   4. RMSE
#   5. % Overestimation (MA prediction > REP estimate)
#   6. Coverage of 95% KDE-based prediction intervals
#   7. Mean PI width
################################################################################

devtools::load_all("/Users/yili/Desktop/Meta Analysis and Replication Forecasting/robumeta-master")
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

set.seed(3104)

# Fixed parameters 
beta0  <- 0.20   # true intercept
sigma2 <- 0.04   # within-study variance (SE = 0.20 on Cohen's d scale)
beta_scales <- c(1, 0.5, 0.75, 0.5, 0.25)  # beta1:beta2:beta3:beta4:beta5 ratios relative to beta_mod (beta4, beta5 are unmeasured)
omega2 <- 0.01   # variance of random slopes (SD ~ 0.10)
delta  <- 0.05   # systematic REP population shift (scaled by beta_scales)
n_rep  <- 30     # REP studies per replication
B      <- 1000    # replications per condition 

# Simulation grid 
sim_grid <- expand.grid(
  k        = c(20, 50, 100),
  tau2_c   = c(0.04, 0.09, 0.25),
  eta      = c(1, 5, 10),
  beta_mod = c(0, 0.20, 0.40),
  stringsAsFactors = FALSE
)

# Data generation 
# Skewed residual population effects u_i with mean 0 and variance tau2_c:
# standardized chi-square (df = 4), skewness sqrt(8 / 4).
r_resid_skewed <- function(n, tau2_c) {
  sqrt(tau2_c) * (rchisq(n, df = 4) - 4) / sqrt(8)
}

# z1, z2, z3 are measured (used by all fitted models); 
# z4, z5 are unmeasured (generated here to drive theta, but never passed to any fitted model).
generate_ma <- function(k, tau2_c, eta, beta0, beta1, beta2, beta3, beta4, beta5,
                        sigma2, omega2) {
  n_ma   <- max(k * 20L, 500L)
  z1     <- rbinom(n_ma, 1L, 0.5)
  z2     <- rnorm(n_ma)
  z3     <- rnorm(n_ma)
  z4     <- rbinom(n_ma, 1L, 0.5)
  z5     <- rnorm(n_ma)

  b1_i   <- rnorm(n_ma, beta1, sqrt(omega2))
  b2_i   <- rnorm(n_ma, beta2, sqrt(omega2))
  b3_i   <- rnorm(n_ma, beta3, sqrt(omega2))
  b4_i   <- rnorm(n_ma, beta4, sqrt(omega2))
  b5_i   <- rnorm(n_ma, beta5, sqrt(omega2))

  theta  <- beta0 + z1 * b1_i + z2 * b2_i + z3 * b3_i +
              z4 * b4_i + z5 * b5_i +
              r_resid_skewed(n_ma, tau2_c)
  yi     <- rnorm(n_ma, theta, sqrt(sigma2))
  pval   <- 2 * (1 - pnorm(abs(yi) / sqrt(sigma2)))
  affirm <- (yi > 0) & (pval < 0.05)

  keep <- runif(n_ma) <= ifelse(affirm, 1, 1 / eta)
  idx  <- which(keep)
  if (length(idx) < k) return(NULL)

  idx <- idx[seq_len(k)]
  tibble(
    yi       = yi[idx],
    vi       = sigma2,
    z1       = z1[idx],
    z2       = z2[idx],
    z3       = z3[idx],
    z4       = z4[idx],
    z5       = z5[idx],
    affirm   = affirm[idx],
    study_id = seq_len(k)
  )
}

generate_rep <- function(n_rep, tau2_c, beta0, beta1, beta2, beta3, beta4, beta5,
                         sigma2, omega2, delta, beta_scales) {
  z1    <- rbinom(n_rep, 1L, 0.5)
  z2    <- rnorm(n_rep)
  z3    <- rnorm(n_rep)
  z4    <- rbinom(n_rep, 1L, 0.5)
  z5    <- rnorm(n_rep)

  b1_i  <- rnorm(n_rep, beta1 + delta * beta_scales[1], sqrt(omega2))
  b2_i  <- rnorm(n_rep, beta2 + delta * beta_scales[2], sqrt(omega2))
  b3_i  <- rnorm(n_rep, beta3 + delta * beta_scales[3], sqrt(omega2))
  b4_i  <- rnorm(n_rep, beta4 + delta * beta_scales[4], sqrt(omega2))
  b5_i  <- rnorm(n_rep, beta5 + delta * beta_scales[5], sqrt(omega2))

  theta <- beta0 + z1 * b1_i + z2 * b2_i + z3 * b3_i +
             z4 * b4_i + z5 * b5_i +
             r_resid_skewed(n_rep, tau2_c)
  yi    <- rnorm(n_rep, theta, sqrt(sigma2))
  tibble(yi = yi, vi = sigma2, z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5)
}

# Single replication 
run_one_rep <- function(k, tau2_c, eta, beta1, n_rep,
                        beta0, beta2, beta3, beta4, beta5, sigma2,
                        omega2, delta, beta_scales) {

  dma <- generate_ma(k, tau2_c, eta, beta0, beta1, beta2, beta3, beta4, beta5,
                     sigma2, omega2)
  if (is.null(dma)) return(NULL)
  dr <- generate_rep(n_rep, tau2_c, beta0, beta1, beta2, beta3, beta4, beta5,
                     sigma2, omega2, delta, beta_scales)

  X_rep <- cbind(1, dr$z1, dr$z2, dr$z3)

  # Step 1: Fit dmaMod (intercept-only, naive) 
  dmaMod <- tryCatch(
    robu(yi ~ 1,
         data         = dma,
         studynum     = as.factor(study_id),
         var.eff.size = vi,
         modelweights = "CORR",  # one estimate per study: standard RE weights
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(dmaMod)) return(NULL)

  tau2_ma <- dmaMod$mod_info$tau.sq
  mu_ma   <- as.numeric(dmaMod$b.r)
  se_ma   <- as.numeric(dmaMod$reg_table$SE)
  df_ma   <- as.numeric(dmaMod$reg_table$dfs)
  if (is.na(tau2_ma) || is.na(mu_ma)) return(NULL)

  # Step 2: Fit mrMod2 (MR with z1 + z2 + z3, naive) 
  mrMod2 <- tryCatch(
    robu(yi ~ z1 + z2 + z3,
         data         = dma,
         studynum     = as.factor(study_id),
         var.eff.size = vi,
         modelweights = "CORR",  # one estimate per study: standard RE weights
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(mrMod2)) return(NULL)

  tau2_mr      <- mrMod2$mod_info$tau.sq
  beta_hat_mr2 <- as.numeric(mrMod2$b.r)
  if (is.na(tau2_mr) || any(is.na(beta_hat_mr2))) return(NULL)

  # Step 3: Fit metaCorr (intercept-only, IPPW-corrected) 
  w_ippw_ma <- ifelse(dma$affirm, 1, eta) / (dma$vi + tau2_ma)

  metaCorr <- tryCatch(
    robu(yi ~ 1,
         data         = dma,
         studynum     = as.factor(study_id),
         userweights  = w_ippw_ma,
         var.eff.size = vi,
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(metaCorr)) return(NULL)

  mu_mc <- as.numeric(metaCorr$b.r)
  se_mc <- as.numeric(metaCorr$reg_table$SE)
  df_mc <- as.numeric(metaCorr$reg_table$dfs)
  if (is.na(mu_mc)) return(NULL)

  # Step 4: Fit mrMod3 (MR with z1 + z2 + z3, IPPW-corrected) 
  w_ippw_mr <- ifelse(dma$affirm, 1, eta) / (dma$vi + tau2_mr)

  mrMod3 <- tryCatch(
    robu(yi ~ z1 + z2 + z3,
         data         = dma,
         studynum     = as.factor(study_id),
         userweights  = w_ippw_mr,
         var.eff.size = vi,
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(mrMod3)) return(NULL)

  beta_hat_mr3 <- as.numeric(mrMod3$b.r)
  if (any(is.na(beta_hat_mr3))) return(NULL)

  # Step 5: Fit REP models for raw difference (Metric 1) 
  dr$study_id <- seq_len(n_rep)

  repMod_int <- tryCatch(
    robu(yi ~ 1,
         data         = dr,
         studynum     = as.factor(study_id),
         var.eff.size = vi,
         modelweights = "CORR",
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(repMod_int)) return(NULL)

  repMod_mr <- tryCatch(
    robu(yi ~ z1 + z2 + z3,
         data         = dr,
         studynum     = as.factor(study_id),
         var.eff.size = vi,
         modelweights = "CORR",
         small        = TRUE),
    error = function(e) NULL
  )
  if (is.null(repMod_mr)) return(NULL)

  mu_rep_int   <- as.numeric(repMod_int$b.r)
  beta_hat_rep <- as.numeric(repMod_mr$b.r)

  pred_rep_at_mean <- beta_hat_rep[1] +
    mean(dr$z1) * beta_hat_rep[2] +
    mean(dr$z2) * beta_hat_rep[3] +
    mean(dr$z3) * beta_hat_rep[4]

  pred_mr2_at_mean <- beta_hat_mr2[1] +
    mean(dma$z1) * beta_hat_mr2[2] +
    mean(dma$z2) * beta_hat_mr2[3] +
    mean(dma$z3) * beta_hat_mr2[4]

  pred_mr3_at_mean <- beta_hat_mr3[1] +
    mean(dma$z1) * beta_hat_mr3[2] +
    mean(dma$z2) * beta_hat_mr3[3] +
    mean(dma$z3) * beta_hat_mr3[4]

  raw_diff_dma <- mu_rep_int - mu_ma
  raw_diff_mr2 <- pred_rep_at_mean - pred_mr2_at_mean
  raw_diff_mc  <- mu_rep_int - mu_mc
  raw_diff_mr3 <- pred_rep_at_mean - pred_mr3_at_mean

  # Step 6: Standard parametric prediction intervals 
  # PI = point_est +/- t_{df} * sqrt(SE^2 + tau2 + vi)

  # MA-U (dmaMod): intercept-only, shared PI for all REP studies
  crit_dma     <- qt(0.975, df = df_ma)
  totalVar_dma <- se_ma^2 + tau2_ma + dr$vi
  pred_dma     <- rep(mu_ma, n_rep)
  pi_lo_dma    <- mu_ma - crit_dma * sqrt(totalVar_dma)
  pi_hi_dma    <- mu_ma + crit_dma * sqrt(totalVar_dma)

  # MA-C (metaCorr): intercept-only, tau2 from naive model
  crit_mc     <- qt(0.975, df = df_mc)
  totalVar_mc <- se_mc^2 + tau2_ma + dr$vi
  pred_mc     <- rep(mu_mc, n_rep)
  pi_lo_mc    <- mu_mc - crit_mc * sqrt(totalVar_mc)
  pi_hi_mc    <- mu_mc + crit_mc * sqrt(totalVar_mc)

  # MR-U (mrMod2): prediction and SE 
  pred_mr2   <- as.numeric(X_rep %*% beta_hat_mr2)
  se_mr2     <- sqrt(pmax(rowSums((X_rep %*% mrMod2$VR.r) * X_rep), 0))
  df_mr2     <- min(mrMod2$reg_table$dfs)
  crit_mr2   <- qt(0.975, df = df_mr2)
  totalVar_mr2 <- se_mr2^2 + tau2_mr + dr$vi
  pi_lo_mr2  <- pred_mr2 - crit_mr2 * sqrt(totalVar_mr2)
  pi_hi_mr2  <- pred_mr2 + crit_mr2 * sqrt(totalVar_mr2)

  # MR-C (mrMod3): prediction and SE, tau2 from naive MR
  pred_mr3   <- as.numeric(X_rep %*% beta_hat_mr3)
  se_mr3     <- sqrt(pmax(rowSums((X_rep %*% mrMod3$VR.r) * X_rep), 0))
  df_mr3     <- min(mrMod3$reg_table$dfs)
  crit_mr3   <- qt(0.975, df = df_mr3)
  totalVar_mr3 <- se_mr3^2 + tau2_mr + dr$vi
  pi_lo_mr3  <- pred_mr3 - crit_mr3 * sqrt(totalVar_mr3)
  pi_hi_mr3  <- pred_mr3 + crit_mr3 * sqrt(totalVar_mr3)

  # Step 7: P_orig (Mathur & VanderWeele, 2020)
  # P_orig_j = 2 * [1 - Phi(|y_j - pred_j| / sqrt(tau2 + SE(pred_j)^2 + v_j))]
  # tau2 for the corrected models is borrowed from the uncorrected models.
  p_orig <- function(pred, tau2, se_pred) {
    2 * (1 - pnorm(abs(dr$yi - pred) / sqrt(tau2 + se_pred^2 + dr$vi)))
  }

  p_orig_dma <- p_orig(pred_dma, tau2_ma, se_ma)
  p_orig_mc  <- p_orig(pred_mc,  tau2_ma, se_mc)
  p_orig_mr2 <- p_orig(pred_mr2, tau2_mr, se_mr2)
  p_orig_mr3 <- p_orig(pred_mr3, tau2_mr, se_mr3)

  # Aggregate per model 
  make_row <- function(model, raw_diff, pred, p_orig, lo, hi) {
    tibble(
      model       = model,
      raw_diff    = raw_diff,
      mean_p_orig = mean(p_orig),
      mae         = mean(abs(dr$yi - pred)),
      rmse        = sqrt(mean((dr$yi - pred)^2)),
      pct_overest = mean(pred > dr$yi),
      coverage    = mean(dr$yi >= lo & dr$yi <= hi),
      pi_width    = mean(hi - lo)
    )
  }

  bind_rows(
    make_row("dmaMod",   raw_diff_dma, pred_dma, p_orig_dma, pi_lo_dma, pi_hi_dma),
    make_row("mrMod2",   raw_diff_mr2, pred_mr2, p_orig_mr2, pi_lo_mr2, pi_hi_mr2),
    make_row("metaCorr", raw_diff_mc,  pred_mc,  p_orig_mc,  pi_lo_mc,  pi_hi_mc),
    make_row("mrMod3",   raw_diff_mr3, pred_mr3, p_orig_mr3, pi_lo_mr3, pi_hi_mr3)
  )
}

# Simulation loop 
model_levels <- c("dmaMod", "mrMod2", "metaCorr", "mrMod3")
results_list <- vector("list", nrow(sim_grid))

for (s in seq_len(nrow(sim_grid))) {

  k        <- sim_grid$k[s]
  tau2_c   <- sim_grid$tau2_c[s]
  eta      <- sim_grid$eta[s]
  beta_mod <- sim_grid$beta_mod[s]
  beta1    <- beta_mod * beta_scales[1]
  beta2    <- beta_mod * beta_scales[2]
  beta3    <- beta_mod * beta_scales[3]
  beta4    <- beta_mod * beta_scales[4]
  beta5    <- beta_mod * beta_scales[5]

  cat(sprintf("[%02d/%02d]  k = %3d  tau2_c = %.2f  eta = %d  beta_mod = %.2f\n",
              s, nrow(sim_grid), k, tau2_c, eta, beta_mod))

  reps <- map(
    seq_len(B),
    ~ run_one_rep(k, tau2_c, eta, beta1, n_rep,
                  beta0, beta2, beta3, beta4, beta5, sigma2,
                  omega2, delta, beta_scales)
  )

  reps_ok  <- compact(reps)
  n_valid  <- length(reps_ok)
  n_failed <- B - n_valid

  if (n_failed > 0)
    message(sprintf("  %d/%d replications failed (model error or insufficient retained studies)",
                    n_failed, B))

  if (n_valid == 0) {
    warning(sprintf("Condition %d: all replications failed — skipping", s))
    next
  }

  results_list[[s]] <- bind_rows(reps_ok) |>
    mutate(model = factor(model, levels = model_levels)) |>
    group_by(model) |>
    summarise(
      across(everything(),
             list(mean = \(x) mean(x, na.rm = TRUE),
                  sd   = \(x) sd(x,   na.rm = TRUE))),
      .groups = "drop"
    ) |>
    mutate(k = k, tau2_c = tau2_c, eta = eta,
           beta_mod = beta_mod, n_valid = n_valid)
}

sim_results <- bind_rows(results_list)

# Save results 
write.csv(sim_results, "sim_results_param_skewed.csv", row.names = FALSE)


# Summary tables 

# Summary 1: average by eta × tau2_c (over k and beta_mod)
sim_summary1 <- sim_results |>
  group_by(eta, tau2_c, model) |>
  summarise(
    raw_diff    = mean(raw_diff_mean,    na.rm = TRUE),
    mae         = mean(mae_mean,         na.rm = TRUE),
    rmse        = mean(rmse_mean,        na.rm = TRUE),
    mean_p_orig = mean(mean_p_orig_mean, na.rm = TRUE) * 100,
    pct_overest = mean(pct_overest_mean, na.rm = TRUE) * 100,
    coverage    = mean(coverage_mean,    na.rm = TRUE) * 100,
    pi_width    = mean(pi_width_mean,    na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = model,
    values_from = c(raw_diff, mae, rmse, mean_p_orig, pct_overest, coverage, pi_width),
    names_glue  = "{model}_{.value}"
  ) |>
  select(eta, tau2_c,
         dmaMod_raw_diff, mrMod2_raw_diff, metaCorr_raw_diff, mrMod3_raw_diff,
         dmaMod_mae,      mrMod2_mae,      metaCorr_mae,      mrMod3_mae,
         dmaMod_rmse,     mrMod2_rmse,     metaCorr_rmse,     mrMod3_rmse,
         dmaMod_mean_p_orig, mrMod2_mean_p_orig, metaCorr_mean_p_orig, mrMod3_mean_p_orig,
         dmaMod_pct_overest, mrMod2_pct_overest, metaCorr_pct_overest, mrMod3_pct_overest,
         dmaMod_coverage,    mrMod2_coverage,    metaCorr_coverage,    mrMod3_coverage,
         dmaMod_pi_width,    mrMod2_pi_width,    metaCorr_pi_width,    mrMod3_pi_width) |>
  arrange(eta, tau2_c)

write.csv(sim_summary1, "sim_summary_mean_param_skewed.csv", row.names = FALSE)

# Summary 2. Win-rate: per condition, find best model per metric
sim_wins <- sim_results |>
  group_by(k, tau2_c, eta, beta_mod) |>
  summarise(
    win_raw_diff    = model[which.min(abs(raw_diff_mean))],
    win_mae         = model[which.min(mae_mean)],
    win_rmse        = model[which.min(rmse_mean)],
    win_mean_p_orig = model[which.max(mean_p_orig_mean)],
    win_pct_overest = model[which.min(abs(pct_overest_mean - 0.5))],
    win_coverage    = model[which.min(abs(coverage_mean - 0.95))],
    win_pi_width    = model[which.min(pi_width_mean)],
    .groups = "drop"
  )

sim_winrate <- bind_rows(
  sim_wins |> count(Model = win_raw_diff)    |> mutate(Metric = "RawDiff"),
  sim_wins |> count(Model = win_mae)         |> mutate(Metric = "MAE"),
  sim_wins |> count(Model = win_rmse)        |> mutate(Metric = "RMSE"),
  sim_wins |> count(Model = win_mean_p_orig) |> mutate(Metric = "MeanP_orig"),
  sim_wins |> count(Model = win_pct_overest) |> mutate(Metric = "PctOverest"),
  sim_wins |> count(Model = win_coverage)    |> mutate(Metric = "Coverage"),
  sim_wins |> count(Model = win_pi_width)    |> mutate(Metric = "PIWidth")
) |>
  pivot_wider(names_from = Metric, values_from = n, values_fill = 0) |>
  mutate(Model = factor(Model, levels = model_levels)) |>
  arrange(Model) |>
  select(Model, RawDiff, MAE, RMSE, MeanP_orig, PctOverest, Coverage, PIWidth)

write.csv(sim_winrate, "sim_summary_winrate_param_skewed.csv", row.names = FALSE)

# Summary 3. Win-rate by moderator strength group
beta_sum <- sim_results |>
  distinct(beta_mod) |>
  mutate(beta_sum = beta_mod * sum(beta_scales))

nonzero_median <- median(beta_sum$beta_sum[beta_sum$beta_sum > 0])

sim_wins_beta <- sim_results |>
  mutate(
    beta_sum   = beta_mod * sum(beta_scales),
    beta_group = case_when(
      beta_sum == 0             ~ "None",
      beta_sum <= nonzero_median ~ "Weak",
      TRUE                       ~ "Strong"
    ),
    beta_group = factor(beta_group, levels = c("None", "Weak", "Strong"))
  ) |>
  group_by(beta_group, k, tau2_c, eta, beta_mod) |>
  summarise(
    win_raw_diff     = model[which.min(abs(raw_diff_mean))],
    win_mae          = model[which.min(mae_mean)],
    win_rmse         = model[which.min(rmse_mean)],
    win_mean_p_orig  = model[which.max(mean_p_orig_mean)],
    win_pct_overest = model[which.min(abs(pct_overest_mean - 0.5))],
    win_coverage    = model[which.min(abs(coverage_mean - 0.95))],
    win_pi_width    = model[which.min(pi_width_mean)],
    .groups = "drop"
  )

sim_winrate_beta <- bind_rows(
  sim_wins_beta |> count(beta_group, Model = win_raw_diff)    |> mutate(Metric = "RawDiff"),
  sim_wins_beta |> count(beta_group, Model = win_mae)         |> mutate(Metric = "MAE"),
  sim_wins_beta |> count(beta_group, Model = win_rmse)        |> mutate(Metric = "RMSE"),
  sim_wins_beta |> count(beta_group, Model = win_mean_p_orig) |> mutate(Metric = "MeanP_orig"),
  sim_wins_beta |> count(beta_group, Model = win_pct_overest) |> mutate(Metric = "PctOverest"),
  sim_wins_beta |> count(beta_group, Model = win_coverage)    |> mutate(Metric = "Coverage"),
  sim_wins_beta |> count(beta_group, Model = win_pi_width)    |> mutate(Metric = "PIWidth")
) |>
  pivot_wider(names_from = Metric, values_from = n, values_fill = 0) |>
  mutate(Model = factor(Model, levels = model_levels)) |>
  arrange(beta_group, Model) |>
  select(beta_group, Model, RawDiff, MAE, RMSE, MeanP_orig,
         PctOverest, Coverage, PIWidth)

# write.csv(sim_winrate_beta, "sim_summary_winrate_beta_param_skewed.csv", row.names = FALSE)
