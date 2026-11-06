############################################################
# Replication Forecasting from Meta-Analysis 
############################################################


############################################################
# 1. Setup
############################################################

library(tidyverse)
library(knitr)
library(here)
library(tableone)
library(corrr)
library(MetaUtility)
library(fastDummies)
library(weightr)
library(PublicationBias)
library(xtable)
library(boot)
library(testthat)
library(data.table)
library(ggplot2)
library(Replicate)
library(clubSandwich)
library(devtools)
library(metafor)

devtools::load_all("/Users/yili/Desktop/Meta Analysis and Replication Forecasting/robumeta-master")

setwd("/Users/yili/Desktop/Meta Analysis and Replication Forecasting/MB4")

dma = fread("mb4_dma.csv")
dr0 <- fread("REP_full.csv")
dr0 <- dr0[!is.na(dr0$n_excluded_1), ]
dr <- dr0 %>% mutate(yi = d_calc,
                     vi = d_var_calc,
                     study_id = lab_id)

dma_reduc <- dma %>% select(study_id, yi, vi, pval) %>% mutate(sei = sqrt(vi))
dr_reduc <- dr %>% select(study_id, yi, vi, pval) %>% mutate(sei = sqrt(vi))

dma %>% distinct(study_id)
dr %>% distinct(study_id)

significance_funnel(yi = dma$yi, vi = dma$vi, favor_positive = TRUE)
ggsave("outputs/MB4_dma_funnel.png", width = 6, height = 6, dpi = 500)

summary(dma$yi)

# Selection model
# First, fit a standard random-effects meta-analysis
res <- rma(yi, vi, data = dma, method = "REML")
# unadjusted estimate:
c(res$b, res$ci.lb, res$ci.ub)   
# 0.3296498 0.2644624 0.3948372

# Fit a step-function selection model
# type = "stepfun" with a step at p = 0.025, 0.05, 1
selmodel(res, type = "stepfun", 
         steps = c(.025, .05, .5, 1), 
         ptable = TRUE)

sel_res<-selmodel(res, type = "stepfun", steps = c(.025))

# Summary of the selection model
# selection-adjusted pooled estimate:
c(sel_res$beta, sel_res$ci.lb, sel_res$ci.ub)
# 0.3598673 0.2142890 0.5054456

# selection ratio: 
# delta is the relative likelihood of publication for studies with one-sided 
# p > .025 compared with studies with p <= .025
sel_res$delta  # c(1, delta2)
# 1.288246
c(sel_res$ci.lb.delta[2], sel_res$ci.ub.delta[2]) # 95% CI for delta2
# 0.01813363 2.55835924

# Likelihood test of selection:
# test statistic
sel_res$LRT     
# 0.2471636
# P-value
sel_res$LRTp
# 0.6190794

res$tau2
# 0.01721693
sel_res$tau2
# 0.02728773



############################################################
# 2. Descriptive Statistics
############################################################

cat_vars <- c("response_mode", "exposure_phase", "method", "dependent_measure",
              "participant_design", "infant_type", "scenarios",
              "stimuli", "stimuli_type", "choice_object", "Hamlin_Lab",
              "intent_val", "outcome_val")

##### DMA #####
continuous_stats <- bind_rows(
  dma %>%
    summarise(
      Variable = "Mean Age",
      Level    = NA_character_,
      Value    = sprintf("%.2f (%.2f)", mean(mean_age_1, na.rm = TRUE),
                         sd(mean_age_1, na.rm = TRUE))
    ),
  dma %>%
    summarise(
      Variable = "N Excluded",
      Level    = NA_character_,
      Value    = sprintf("%.2f (%.2f)", mean(n_excluded_1, na.rm = TRUE),
                         sd(n_excluded_1, na.rm = TRUE))
    )
)

categorical_stats <- cat_vars %>%
  map_dfr(~ {
    dma %>%
      count(.data[[.x]]) %>%
      mutate(
        Variable = .x,
        Level    = as.character(.data[[.x]]),
        Value    = sprintf("%d (%.1f%%)", n, n / sum(n) * 100)
      ) %>%
      select(Variable, Level, Value)
  })

desc_table_dma <- bind_rows(continuous_stats, categorical_stats)
write.csv(desc_table_dma, "outputs/MB4_descriptive_dma.csv", row.names = FALSE)


##### DR #####
continuous_stats <- bind_rows(
  dr %>%
    summarise(
      Variable = "Mean Age",
      Level    = NA_character_,
      Value    = sprintf("%.2f (%.2f)", mean(mean_age_1, na.rm = TRUE),
                         sd(mean_age_1, na.rm = TRUE))
    ),
  dr %>%
    summarise(
      Variable = "N Excluded",
      Level    = NA_character_,
      Value    = sprintf("%.2f (%.2f)", mean(n_excluded_1, na.rm = TRUE),
                         sd(n_excluded_1, na.rm = TRUE))
    )
)

categorical_stats <- cat_vars %>%
  map_dfr(~ {
    dr %>%
      count(.data[[.x]]) %>%
      mutate(
        Variable = .x,
        Level    = as.character(.data[[.x]]),
        Value    = sprintf("%d (%.1f%%)", n, n / sum(n) * 100)
      ) %>%
      select(Variable, Level, Value)
  })

desc_table_dr <- bind_rows(continuous_stats, categorical_stats)
write.csv(desc_table_dr, "outputs/MB4_descriptive_dr.csv", row.names = FALSE)


############################################################
# 3. Model 1: MA-U (MAU) — Intercept-only, uncorrected
############################################################

MAU = robu(yi ~ 1,
              data = dma,
              studynum = as.factor(study_id),
              var.eff.size = vi,
              modelweights = "HIER",
              small = TRUE)
MAU
MAU$b.r

k <- nrow(dma)

### DerSimonian-Laird tau2
meta.re = rma.uni(yi = dma$yi, vi = dma$vi, method = "DL")
t2hat.naive = meta.re$tau2
t2hat.naive_dma <- t2hat.naive

### Wang & Lee calibrated estimates
dma_cop1 <- dma %>% select(yi, vi) %>%
  mutate(Qi = (yi - sum(yi / vi) / sum(1 / vi))^2 / vi)

sum_inv_vi <- sum(1 / dma_cop1$vi)
sum_inv_vi2 <- sum(1 / dma_cop1$vi^2)
Q <- sum(dma_cop1$Qi)

tau2 <- (Q - (k - 1)) / (sum_inv_vi - sum_inv_vi2 / sum_inv_vi)

dma_cop2 <- dma_cop1 %>%
  mutate(mui_num = yi / (vi + tau2),
         mui_denom = 1 / (vi + tau2))
mu <- sum(dma_cop2$mui_num) / sum(dma_cop2$mui_denom)
var_mu <- 1 / sum(dma_cop2$mui_denom)

dma_yi_cali_wang <- dma %>%
  select(yi, vi) %>%
  mutate(yi_cali_wang = mu + sqrt(tau2 / (tau2 + vi)) * (yi - mu))

### Wang & Lee prediction interval
dma_order <- dma_yi_cali_wang %>%
  select(yi, vi, yi_cali_wang) %>%
  arrange(yi_cali_wang) %>%
  mutate(i = 1:k,
         p = i / (k + 1),
         yi_cali_wang_lead = case_when(i < k ~ lead(yi_cali_wang),
                                       TRUE ~ 999))
dma_order <- dma_order %>%
  mutate(Q_p = case_when(
    p < 1 / (k + 1) ~ mu + (dma_order$yi_cali_wang[1] - mu) / qnorm(1 / (k + 1)) * qnorm(p),
    p < (i + 1) / (k + 1) ~ yi_cali_wang + (yi_cali_wang_lead - yi_cali_wang) /
      (qnorm((i + 1) / (k + 1)) - qnorm(i / (k + 1))) * (qnorm(p) - qnorm(i / (k + 1))),
    p > k / (k + 1) ~ mu + (yi_cali_wang - mu) / qnorm(k / (k + 1)) * qnorm(p))
  )

quantile(dma_order$Q_p, c(0.025, 0.975))

### KDE prediction distribution. Bandwidth: Sheather–Jones method
h <- bw.SJ(dma_order$yi_cali_wang)

dma_pred_dist <- tibble(x = seq(-2, 2, 0.005), fh = NA)
for (i in 1:nrow(dma_pred_dist)) {
  dma_pred_dist$fh[i] <- sum(dnorm((dma_pred_dist$x[i] - dma_order$yi_cali_wang) / h)) / (k * h)
}

ggplot(dma_pred_dist, aes(x = x, y = fh)) +
  geom_point(size = 0.5, color = "blue") +
  geom_point(data = dma_yi_cali_wang,
             aes(x = yi_cali_wang, y = 0),
             size = 3, shape = 18, color = "orange", alpha = 1) +
  stat_function(fun = dnorm, n = 601, color = "darkgray", size = 1,
                args = list(mean = mu, sd = sqrt(tau2 + var_mu))) +
  scale_x_continuous(breaks = seq(-2, 4, 1)) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1)) +
  ggtitle("Predictive distributions (Wang & Lee 2018 method). Dataset = dma.")


############################################################
# 4. Model 2: MR-U (MRU) — Meta-regression, uncorrected
############################################################

MRU = robu(yi ~ mean_age_1 + n_excluded_1 +
                method + scenarios + stimuli + stimuli_type +
                choice_object + Hamlin_Lab + intent_val + outcome_val,
              data = dma,
              studynum = as.factor(study_id),
              var.eff.size = vi,
              modelweights = "HIER",
              small = TRUE)
MRU

tau2_MRU <- MRU$mod_info$tau.sq

# Design matrix for MA studies (each MA study's own moderator values)
Zi_mat = as.matrix(data.frame(
  X.Intercept. = 1,
  mean_age_1 = dma$mean_age_1,
  n_excluded_1 = dma$n_excluded_1,
  methodreach = dma$method == "reach",
  scenariosgive_take = dma$scenarios == "give_take",
  scenarioshelp_hinder = dma$scenarios == "help_hinder",
  scenariosother = dma$scenarios == "other",
  stimulimovies = dma$stimuli == "movies",
  stimuli_typereal = dma$stimuli_type == "real",
  choice_objectpeople = dma$choice_object == "people",
  choice_objectpuppets = dma$choice_object == "puppets",
  choice_objectshapes = dma$choice_object == "shapes",
  Hamlin_LabTRUE = dma$Hamlin_Lab == "TRUE",
  intent_valopp = dma$intent_val == "opp",
  intent_valpos = dma$intent_val == "pos",
  outcome_valopp = dma$outcome_val == "opp",
  outcome_valrev = dma$outcome_val == "rev",
  outcome_valsame = dma$outcome_val == "same"
))



############################################################
# 5. Model 3: MA-C (MAC) — MA, corrected for pub bias
############################################################
# Fitted as in the simulation (sim_KDE.R / sim_param.R): intercept-only robu
# with IPPW weights; tau2 in the weights is borrowed from MA-U.

eta = 4.70   # selection ratio estimated for Metalab (Mathur & VanderWeele, 2021)

# Affirmative = positive point estimate AND two-sided p < 0.05.
# The pval column in mb4_dma.csv comes from the original studies' own tests
# (it does not match a z-test on yi and vi), so recompute the two-sided
# Wald p-value from yi and vi, as in the simulation, MB1, and the funnel plot.
dma$pval = 2 * (1 - pnorm(abs(dma$yi) / sqrt(dma$vi)))
dma$affirm = dma$yi > 0 & dma$pval < 0.05
table(dma$affirm)

tau2_MAU <- as.vector(MAU$mod_info$tau.sq)
w_MAC    <- ifelse(dma$affirm, 1, eta) / (dma$vi + tau2_MAU)

MAC = robu(yi ~ 1,
              data = dma,
              studynum = as.factor(study_id),
              userweights = w_MAC,
              var.eff.size = vi,
              small = TRUE)
MAC

MAC_est <- as.numeric(MAC$b.r)
MAC_se  <- as.numeric(MAC$reg_table$SE)
MAC_df  <- as.numeric(MAC$reg_table$dfs)


############################################################
# 6. Model 4: MR-C (MRC) — Meta-regression, corrected
############################################################
# IPPW weights; tau2 in the weights is borrowed from MR-U

mods = c("mean_age_1", "n_excluded_1",
         "method", "scenarios", "stimuli", "stimuli_type",
         "choice_object", "Hamlin_Lab", "intent_val", "outcome_val")

gotError = TRUE
while (gotError == TRUE) {
  linPred = paste("yi ~", paste(mods, collapse = " + "))

  tryCatch({
    MRU = robu(eval(parse(text = linPred)),
                  data = dma,
                  studynum = as.factor(study_id),
                  var.eff.size = vi,
                  modelweights = "HIER",
                  small = TRUE)
    gotError = FALSE
  }, error = function(err) {
    gotError <<- TRUE
    cat(paste("\n Removing",
              mods[length(mods)],
              "from moderators and trying again"))
    mods <<- mods[1:(length(mods) - 1)]
  })
}

mods
tau2_MRU <- MRU$mod_info$tau.sq

w_MRC <- ifelse(dma$affirm, 1, eta) / (dma$vi + as.vector(tau2_MRU))

MRC = robu(eval(parse(text = linPred)),
              studynum = as.factor(study_id),
              data = dma,
              userweights = w_MRC,
              var.eff.size = vi,
              small = TRUE)
MRC

# Design matrix for REP studies (each REP study's own moderator values)
Xmat = as.matrix(data.frame(
  X.Intercept. = 1,
  mean_age_1 = dr$mean_age_1,
  n_excluded_1 = dr$n_excluded_1,
  methodreach = dr$method == "reach",
  scenariosgive_take = dr$scenarios == "give_take",
  scenarioshelp_hinder = dr$scenarios == "help_hinder",
  scenariosother = dr$scenarios == "other",
  stimulimovies = dr$stimuli == "movies",
  stimuli_typereal = dr$stimuli_type == "real",
  choice_objectpeople = dr$choice_object == "people",
  choice_objectpuppets = dr$choice_object == "puppets",
  choice_objectshapes = dr$choice_object == "shapes",
  Hamlin_LabTRUE = dr$Hamlin_Lab == "TRUE",
  intent_valopp = dr$intent_val == "opp",
  intent_valpos = dr$intent_val == "pos",
  outcome_valopp = dr$outcome_val == "opp",
  outcome_valrev = dr$outcome_val == "rev",
  outcome_valsame = dr$outcome_val == "same"
))

# MA studies with missing moderators are dropped by robu from the MR models;
# the MR calibrated estimates use the same complete cases
sum(!complete.cases(Zi_mat))
stopifnot(!anyNA(Xmat))

expect_equal(attr(Xmat, "dimnames")[[2]], MRC$labels)
expect_equal(attr(Zi_mat, "dimnames")[[2]], MRC$labels)


############################################################
# 7. KDE Prediction Intervals (all 4 models)
############################################################
# Same construction as the simulation: convolving the weighted
# KDE(calibrated estimates, h) with each REP study's sampling distribution
# N(0, v_j) gives the mixture
#   f(t) = sum_i w_i * phi(t; calib_i, h^2 + v_j) / sum_i w_i,
# with w_i = 1 for the uncorrected models and IPPW weights
# w_i = eta^I(A_i = 0) for the corrected models, so that the corrected KDE
# represents all conducted studies rather than only the published ones.
# The PI is the 2.5th and 97.5th percentiles, evaluated on a grid.
# Bandwidth: Sheather-Jones on each model's own (unweighted) calibrated
# estimates (Silverman's rule of thumb if SJ fails).

grid_lo   <- -5
grid_hi   <-  5
grid_step <-  0.01

compute_pi_KDE <- function(calib, h, sigma_new, w = NULL) {
  if (is.null(w)) w <- rep(1, length(calib))   # equal weights: uncorrected
  t_vec   <- seq(grid_lo, grid_hi, by = grid_step)
  sd_pred <- sqrt(h^2 + sigma_new^2)

  pdf_y <- as.numeric(
    outer(t_vec, calib, function(t, c) dnorm(t, mean = c, sd = sd_pred)) %*%
      (w / sum(w))
  )
  cdf_y <- cumsum(pdf_y * grid_step)

  c(lo = t_vec[which.min(abs(cdf_y - 0.025))],
    hi = t_vec[which.min(abs(cdf_y - 0.975))])
}

bw_SJ_safe <- function(x) tryCatch(bw.SJ(x), error = function(e) bw.nrd0(x))

##### 7a. Intercept-only KDE PI function (MAU, MAC) #####
KDE_interval_func_intercept <- function(mu_hat, tau2, w = NULL) {
  calib <- mu_hat + sqrt(tau2 / (tau2 + dma$vi)) * (dma$yi - mu_hat)
  h_int <- bw_SJ_safe(calib)

  pis <- t(sapply(sqrt(dr$vi),
                  function(s) compute_pi_KDE(calib, h_int, s, w = w)))

  dr_pred <- dr %>%
    select(lab_id, yi, vi) %>%
    mutate(pred_yi_est = mu_hat,
           pred_yi_lo  = pis[, "lo"],
           pred_yi_hi  = pis[, "hi"],
           cover_yi    = as.numeric(yi >= pred_yi_lo & yi <= pred_yi_hi))

  return(list(mean(dr_pred$cover_yi), dr_pred))
}

dr_pred_MAU_out <- KDE_interval_func_intercept(
  mu_hat = as.numeric(MAU$b.r),
  tau2   = tau2_MAU
)
# IPPW weights for the corrected models' KDE (Appendix, Proposition 2)
w_kde <- ifelse(dma$affirm, 1, eta)

dr_pred_MAC_out <- KDE_interval_func_intercept(
  mu_hat = MAC_est,
  tau2   = tau2_MAU,
  w      = w_kde
)

##### 7b. Meta-regression KDE PI function (MRU, MRC) #####
# Calibrated estimates are evaluated at each REP study's own moderator values:
# each MA study's residual is shrunk toward x_j' beta_hat. tau2 from MR-U.
KDE_interval_func <- function(fitted.model, w = NULL) {
  b_hat   <- as.numeric(fitted.model$b.r)
  tau2    <- as.vector(tau2_MRU)
  ok      <- complete.cases(Zi_mat)   # MA studies used by the meta-regression
  shrink  <- sqrt(tau2 / (tau2 + dma$vi[ok]))
  resid   <- dma$yi[ok] - as.numeric(Zi_mat[ok, , drop = FALSE] %*% b_hat)
  pred_dr <- as.numeric(Xmat %*% b_hat)
  w_ok    <- if (is.null(w)) NULL else w[ok]

  # The bandwidth is shift-invariant, so one bandwidth serves all REP studies
  h_mr <- bw_SJ_safe(shrink * resid)

  pis <- t(sapply(seq_len(nrow(dr)), function(j)
    compute_pi_KDE(pred_dr[j] + shrink * resid, h_mr, sqrt(dr$vi[j]),
                   w = w_ok)))

  dr_pred <- dr %>%
    select(lab_id, yi, vi) %>%
    mutate(pred_yi_peak = pred_dr,
           pred_yi_lo   = pis[, "lo"],
           pred_yi_hi   = pis[, "hi"],
           cover_yi     = as.numeric(yi >= pred_yi_lo & yi <= pred_yi_hi))

  return(list(mean(dr_pred$cover_yi), dr_pred))
}

dr_pred_MRU_out <- KDE_interval_func(fitted.model = MRU)
dr_pred_MRC_out <- KDE_interval_func(fitted.model = MRC, w = w_kde)


############################################################
# 8. Assemble Predictions under 4 Models
############################################################

##### 8.1 MA-U (MAU) — KDE PI #####
dr$pred.MAU        <- as.numeric(MAU$b.r)
dr$MAU_pred_yi_lo  <- dr_pred_MAU_out[[2]]$pred_yi_lo
dr$MAU_pred_yi_hi  <- dr_pred_MAU_out[[2]]$pred_yi_hi

##### 8.2 MA-C (MAC) — KDE PI #####
dr$pred.MAC        <- MAC_est
dr$MAC_pred_yi_lo  <- dr_pred_MAC_out[[2]]$pred_yi_lo
dr$MAC_pred_yi_hi  <- dr_pred_MAC_out[[2]]$pred_yi_hi

##### 8.3 MR-U (MRU) — KDE PI #####
dr$pred.MRU   <- dr_pred_MRU_out[[2]]$pred_yi_peak
dr$predLo.MRU <- dr_pred_MRU_out[[2]]$pred_yi_lo
dr$predHi.MRU <- dr_pred_MRU_out[[2]]$pred_yi_hi

summary(dr$pred.MRU)
summary(dr$predLo.MRU)
summary(dr$predHi.MRU)

##### 8.4 MR-C (MRC) — KDE PI #####
dr$pred.MRC   <- dr_pred_MRC_out[[2]]$pred_yi_peak
dr$predLo.MRC <- dr_pred_MRC_out[[2]]$pred_yi_lo
dr$predHi.MRC <- dr_pred_MRC_out[[2]]$pred_yi_hi

summary(dr$pred.MRC)

##### 8.5 Point estimates dataset #####
dr_pred_est <- dr %>%
  select(lab_id, yi, vi,
         pred.MAU,
         pred.MAC,
         pred.MRU,
         pred.MRC)

summary(dr_pred_est$pred.MRU)



# Calibrated estimates comparison
dr_calib_est <- calib_ests(yi = dr$yi,
                           sei = sqrt(dr$d_var_calc),
                           method = "DL")
dr_calib_est <- as.data.frame(dr_calib_est)

dr_pred_est_calib <- cbind(dr$lab_id, dr$d_calc, dr$d_var_calc,
                           dr$pred.MRU, dr_calib_est)
summary(dr$pred.MRU)
summary(dr_calib_est)


##### 8.6 Prediction intervals dataset #####
dr_pred_interval <- dr %>%
  select(lab_id, yi, vi,
         MAU_pred_yi_lo, MAU_pred_yi_hi,
         MAC_pred_yi_lo, MAC_pred_yi_hi,
         predLo.MRU, predHi.MRU,
         predLo.MRC, predHi.MRC) %>%
  mutate(MRU_pred_yi_lo = predLo.MRU,
         MRU_pred_yi_hi = predHi.MRU,
         MRC_pred_yi_lo = predLo.MRC,
         MRC_pred_yi_hi = predHi.MRC)


############################################################
# 8c. Standard Parametric Prediction Intervals (all models)
############################################################
# As in the simulation:
#   PI_j = pred_j +/- t_df * sqrt( SE(pred_j)^2 + tau2 + v_j )
#   - SE(pred_j) from the robust variance of the coefficients
#   - df: Satterthwaite df from robu (MR models: smallest df across
#     coefficients)
#   - tau2 for the corrected models is borrowed from the uncorrected models

##### 8c.1 MA-U (MAU) #####
crit_t_MAU <- qt(.975, df = as.numeric(MAU$reg_table$dfs))
param_totalVar_MAU <- tau2_MAU + as.numeric(MAU$reg_table$SE^2) + dr$vi

dr$param_pred.MAU <- as.numeric(MAU$b.r)
dr$param_MAU_pred_lo <- as.numeric(MAU$b.r) -
  crit_t_MAU * sqrt(param_totalVar_MAU)
dr$param_MAU_pred_hi <- as.numeric(MAU$b.r) +
  crit_t_MAU * sqrt(param_totalVar_MAU)

##### 8c.2 MA-C (MAC) #####
crit_t_MAC <- qt(.975, df = MAC_df)
param_totalVar_MAC <- tau2_MAU + MAC_se^2 + dr$vi

dr$param_pred.MAC <- MAC_est
dr$param_MAC_pred_lo <- MAC_est - crit_t_MAC * sqrt(param_totalVar_MAC)
dr$param_MAC_pred_hi <- MAC_est + crit_t_MAC * sqrt(param_totalVar_MAC)

##### 8c.3–8c.4 MR-U (MRU) and MR-C (MRC) #####
param_PI_MR <- function(fitted.model) {
  pred <- as.numeric(Xmat %*% fitted.model$b.r)
  se   <- sqrt(pmax(rowSums((Xmat %*% fitted.model$VR.r) * Xmat), 0))
  crit <- qt(.975, df = min(fitted.model$reg_table$dfs))
  half <- crit * sqrt(se^2 + as.vector(tau2_MRU) + dr$vi)
  list(pred = pred, lo = pred - half, hi = pred + half)
}

MRU_param <- param_PI_MR(MRU)
dr$param_pred.MRU    <- MRU_param$pred
dr$param_MRU_pred_lo <- MRU_param$lo
dr$param_MRU_pred_hi <- MRU_param$hi

MRC_param <- param_PI_MR(MRC)
dr$param_pred.MRC    <- MRC_param$pred
dr$param_MRC_pred_lo <- MRC_param$lo
dr$param_MRC_pred_hi <- MRC_param$hi

##### 8c.5 Standard PI datasets #####
dr_param_pred_est <- dr %>%
  select(lab_id, yi, vi,
         param_pred.MAU,
         param_pred.MAC,
         param_pred.MRU,
         param_pred.MRC)

dr_param_pred_interval <- dr %>%
  select(lab_id, yi, vi,
         param_MAU_pred_lo, param_MAU_pred_hi,
         param_MAC_pred_lo, param_MAC_pred_hi,
         param_MRU_pred_lo, param_MRU_pred_hi,
         param_MRC_pred_lo, param_MRC_pred_hi)


############################################################
# 9. Outcome Metrics
############################################################

##### 9.1 Raw difference between REP and MA pooled estimates #####
# As in the simulation: REP pooled estimate minus MA estimate.
#   - Intercept-only models: REP pooled mean minus (corrected) MA mean
#   - MR models: REP meta-regression evaluated at the REP covariate means
#     minus MA meta-regression evaluated at the MA covariate means

rep_formula <- yi ~ mean_age_1 + n_excluded_1 + Hamlin_Lab

REP_int = robu(yi ~ 1, data = dr,
                  studynum = as.factor(lab_id),
                  var.eff.size = vi,
                  modelweights = "CORR", small = TRUE)

REP_mr = robu(rep_formula, data = dr,
                 studynum = as.factor(lab_id),
                 var.eff.size = vi,
                 modelweights = "CORR", small = TRUE)

X_rep_full <- model.matrix(rep_formula, data = dr)
stopifnot(ncol(X_rep_full) == length(REP_mr$b.r))

rep_mean    <- as.numeric(REP_int$b.r)
rep_mr_mean <- sum(colMeans(X_rep_full) * as.numeric(REP_mr$b.r))

raw_diff_MAU <- rep_mean - as.numeric(MAU$b.r)
raw_diff_MAC <- rep_mean - MAC_est
ma_mr_means  <- colMeans(Zi_mat[complete.cases(Zi_mat), , drop = FALSE])
raw_diff_MRU <- rep_mr_mean - sum(ma_mr_means * as.numeric(MRU$b.r))
raw_diff_MRC <- rep_mr_mean - sum(ma_mr_means * as.numeric(MRC$b.r))

c(MAU = raw_diff_MAU, MRU = raw_diff_MRU,
  MAC = raw_diff_MAC, MRC = raw_diff_MRC)


##### 9.2 P_orig (statistical consistency) #####
# P_orig (Mathur & VanderWeele, 2020; Equation 8), as in the simulation:
#   P_orig_j = 2 * [1 - Phi(|y_j - pred_j| / sqrt(tau2 + SE(pred_j)^2 + v_j))]
# tau2 for the corrected models is borrowed from the uncorrected models.

p_orig_func <- function(pred, tau2, se_pred) {
  2 * (1 - pnorm(abs(dr$yi - pred) / sqrt(tau2 + se_pred^2 + dr$vi)))
}

se_MR <- function(fitted.model) {
  sqrt(pmax(rowSums((Xmat %*% fitted.model$VR.r) * Xmat), 0))
}

p_orig_MAU <- p_orig_func(as.numeric(MAU$b.r), tau2_MAU,
                          as.numeric(MAU$reg_table$SE))
p_orig_MRU <- p_orig_func(as.numeric(Xmat %*% MRU$b.r), as.vector(tau2_MRU),
                          se_MR(MRU))
p_orig_MAC <- p_orig_func(MAC_est, tau2_MAU, MAC_se)
p_orig_MRC <- p_orig_func(as.numeric(Xmat %*% MRC$b.r), as.vector(tau2_MRU),
                          se_MR(MRC))

summary(p_orig_MAU)
summary(p_orig_MRU)
summary(p_orig_MAC)
summary(p_orig_MRC)

dr_pred_est_porig <- dr_pred_est %>%
  mutate(p_orig_MAU = p_orig_MAU,
         p_orig_MRU = p_orig_MRU,
         p_orig_MAC = p_orig_MAC,
         p_orig_MRC = p_orig_MRC)

apply(dr_pred_est_porig[,8:11], 2, mean)


##### 9.3 MAE #####
dr_pred_est %>% summarize(
  MAE_pred.MAU = mean(abs(yi - pred.MAU)),
  MAE_pred.MRU = mean(abs(yi - pred.MRU)),
  MAE_pred.MAC = mean(abs(yi - pred.MAC)),
  MAE_pred.MRC = mean(abs(yi - pred.MRC)))

##### 9.4 RMSE #####
dr_pred_est %>% summarize(
  RMSE_pred.MAU = sqrt(mean((yi - pred.MAU)^2)),
  RMSE_pred.MRU = sqrt(mean((yi - pred.MRU)^2)),
  RMSE_pred.MAC = sqrt(mean((yi - pred.MAC)^2)),
  RMSE_pred.MRC = sqrt(mean((yi - pred.MRC)^2)))

##### 9.5 % Overestimate #####
dr_pred_est %>% summarize(
  percOverest_pred.MAU = mean(yi < pred.MAU),
  percOverest_pred.MRU = mean(yi < pred.MRU),
  percOverest_pred.MAC = mean(yi < pred.MAC),
  percOverest_pred.MRC = mean(yi < pred.MRC))

##### 9.6 Coverage of prediction intervals #####
dr_pred_interval %>% summarize(
  percCover_yi_pred.MAU = mean(ifelse((yi - MAU_pred_yi_hi)*(yi - MAU_pred_yi_lo) <= 0, 1, 0)),
  percCover_yi_pred.MRU = mean(ifelse((yi - MRU_pred_yi_hi)*(yi - MRU_pred_yi_lo) <= 0, 1, 0)),
  percCover_yi_pred.MAC = mean(ifelse((yi - MAC_pred_yi_hi)*(yi - MAC_pred_yi_lo) <= 0, 1, 0)),
  percCover_yi_pred.MRC = mean(ifelse((yi - MRC_pred_yi_hi)*(yi - MRC_pred_yi_lo) <= 0, 1, 0)))

# Width of prediction intervals
dr_pred_interval %>% summarize(
  width_yi_pred.MAU = mean(MAU_pred_yi_hi - MAU_pred_yi_lo),
  width_yi_pred.MRU = mean(MRU_pred_yi_hi - MRU_pred_yi_lo),
  width_yi_pred.MAC = mean(MAC_pred_yi_hi - MAC_pred_yi_lo),
  width_yi_pred.MRC = mean(MRC_pred_yi_hi - MRC_pred_yi_lo))


##### 9.7 Unified KDE Metrics Table #####
kde_metrics_table <- tibble(
  Metric = c("Raw difference", "Mean P_orig", "MAE", "RMSE",
             "% Overestimate", "PI Coverage", "PI Width"),

  `MA-U` = c(
    raw_diff_MAU,
    mean(p_orig_MAU),
    mean(abs(dr_pred_est$yi - dr_pred_est$pred.MAU)),
    sqrt(mean((dr_pred_est$yi - dr_pred_est$pred.MAU)^2)),
    mean(dr_pred_est$yi < dr_pred_est$pred.MAU),
    mean(dr_pred_interval$yi >= dr_pred_interval$MAU_pred_yi_lo &
           dr_pred_interval$yi <= dr_pred_interval$MAU_pred_yi_hi),
    mean(dr_pred_interval$MAU_pred_yi_hi - dr_pred_interval$MAU_pred_yi_lo)
  ),

  `MR-U` = c(
    raw_diff_MRU,
    mean(p_orig_MRU),
    mean(abs(dr_pred_est$yi - dr_pred_est$pred.MRU)),
    sqrt(mean((dr_pred_est$yi - dr_pred_est$pred.MRU)^2)),
    mean(dr_pred_est$yi < dr_pred_est$pred.MRU),
    mean(dr_pred_interval$yi >= dr_pred_interval$MRU_pred_yi_lo &
           dr_pred_interval$yi <= dr_pred_interval$MRU_pred_yi_hi),
    mean(dr_pred_interval$MRU_pred_yi_hi - dr_pred_interval$MRU_pred_yi_lo)
  ),

  `MA-C` = c(
    raw_diff_MAC,
    mean(p_orig_MAC),
    mean(abs(dr_pred_est$yi - dr_pred_est$pred.MAC)),
    sqrt(mean((dr_pred_est$yi - dr_pred_est$pred.MAC)^2)),
    mean(dr_pred_est$yi < dr_pred_est$pred.MAC),
    mean(dr_pred_interval$yi >= dr_pred_interval$MAC_pred_yi_lo &
           dr_pred_interval$yi <= dr_pred_interval$MAC_pred_yi_hi),
    mean(dr_pred_interval$MAC_pred_yi_hi - dr_pred_interval$MAC_pred_yi_lo)
  ),

  `MR-C` = c(
    raw_diff_MRC,
    mean(p_orig_MRC),
    mean(abs(dr_pred_est$yi - dr_pred_est$pred.MRC)),
    sqrt(mean((dr_pred_est$yi - dr_pred_est$pred.MRC)^2)),
    mean(dr_pred_est$yi < dr_pred_est$pred.MRC),
    mean(dr_pred_interval$yi >= dr_pred_interval$MRC_pred_yi_lo &
           dr_pred_interval$yi <= dr_pred_interval$MRC_pred_yi_hi),
    mean(dr_pred_interval$MRC_pred_yi_hi - dr_pred_interval$MRC_pred_yi_lo)
  )
)

print(kde_metrics_table)
write.csv(kde_metrics_table, "outputs/MB4_KDE_PI_metrics.csv", row.names = FALSE)


############################################################
# 10. Outcome Metrics Table — Standard Parametric PIs
############################################################
# Point estimates are the same as in the KDE table, so Raw difference,
# Mean P_orig, MAE, RMSE, and % Overestimate are unchanged; only PI
# coverage and width differ.

param_metrics_table <- tibble(
  Metric = c("Raw difference", "Mean P_orig", "MAE", "RMSE",
             "% Overestimate", "PI Coverage", "PI Width"),

  `MA-U` = c(
    raw_diff_MAU,
    mean(p_orig_MAU),
    mean(abs(dr_param_pred_est$yi - dr_param_pred_est$param_pred.MAU)),
    sqrt(mean((dr_param_pred_est$yi - dr_param_pred_est$param_pred.MAU)^2)),
    mean(dr_param_pred_est$yi < dr_param_pred_est$param_pred.MAU),
    mean(dr_param_pred_interval$yi >= dr_param_pred_interval$param_MAU_pred_lo &
           dr_param_pred_interval$yi <= dr_param_pred_interval$param_MAU_pred_hi),
    mean(dr_param_pred_interval$param_MAU_pred_hi -
           dr_param_pred_interval$param_MAU_pred_lo)
  ),

  `MR-U` = c(
    raw_diff_MRU,
    mean(p_orig_MRU),
    mean(abs(dr_param_pred_est$yi - dr_param_pred_est$param_pred.MRU)),
    sqrt(mean((dr_param_pred_est$yi - dr_param_pred_est$param_pred.MRU)^2)),
    mean(dr_param_pred_est$yi < dr_param_pred_est$param_pred.MRU),
    mean(dr_param_pred_interval$yi >= dr_param_pred_interval$param_MRU_pred_lo &
           dr_param_pred_interval$yi <= dr_param_pred_interval$param_MRU_pred_hi),
    mean(dr_param_pred_interval$param_MRU_pred_hi -
           dr_param_pred_interval$param_MRU_pred_lo)
  ),

  `MA-C` = c(
    raw_diff_MAC,
    mean(p_orig_MAC),
    mean(abs(dr_param_pred_est$yi - dr_param_pred_est$param_pred.MAC)),
    sqrt(mean((dr_param_pred_est$yi - dr_param_pred_est$param_pred.MAC)^2)),
    mean(dr_param_pred_est$yi < dr_param_pred_est$param_pred.MAC),
    mean(dr_param_pred_interval$yi >= dr_param_pred_interval$param_MAC_pred_lo &
           dr_param_pred_interval$yi <= dr_param_pred_interval$param_MAC_pred_hi),
    mean(dr_param_pred_interval$param_MAC_pred_hi -
           dr_param_pred_interval$param_MAC_pred_lo)
  ),

  `MR-C` = c(
    raw_diff_MRC,
    mean(p_orig_MRC),
    mean(abs(dr_param_pred_est$yi - dr_param_pred_est$param_pred.MRC)),
    sqrt(mean((dr_param_pred_est$yi - dr_param_pred_est$param_pred.MRC)^2)),
    mean(dr_param_pred_est$yi < dr_param_pred_est$param_pred.MRC),
    mean(dr_param_pred_interval$yi >= dr_param_pred_interval$param_MRC_pred_lo &
           dr_param_pred_interval$yi <= dr_param_pred_interval$param_MRC_pred_hi),
    mean(dr_param_pred_interval$param_MRC_pred_hi -
           dr_param_pred_interval$param_MRC_pred_lo)
  )
)

print(param_metrics_table)
write.csv(param_metrics_table, "outputs/MB4_param_PI_metrics.csv", row.names = FALSE)


############################################################
# 11. Prediction Interval Plots (KDE PIs)
############################################################

########### Prepare data ###########
dr_pred_est <- dr_pred_est %>% arrange(desc(yi))
dr_pred_interval <- dr_pred_interval %>% arrange(desc(yi))
n_studies <- nrow(dr_pred_est)

study_labels <- paste(seq_len(n_studies), dr_pred_est$lab_id, sep = ": ")

shared_theme <- theme_bw() +
  theme(
    axis.title        = element_blank(),
    legend.position   = "right",
    legend.background = element_rect(fill = "white", color = "black",
                                     linewidth = 0.5),
    legend.key        = element_rect(fill = "white"),
    legend.text       = element_text(size = 9),
    legend.title      = element_blank()
  )

legend_guide <- guide_legend(
  override.aes = list(
    shape    = 15,
    linetype = "blank",
    linewidth = 0,
    size     = 4
  )
)


########### Plot 1: MA-U and MA-C ###########
p_ma <- ggplot() +
  geom_point(data = dr_pred_est,
             aes(x = yi, y = seq_len(n_studies),
                 color = "REP effect size"),
             size = 1) +
  
  geom_vline(aes(xintercept = as.numeric(MAU$b.r),
                 color = "MA-U")) +
  
  geom_vline(aes(xintercept = MAC_est,   # MA-C is now a robu fit
                 color = "MA-C")) +
  
  geom_errorbarh(data = dr_pred_interval,
                 aes(xmax = MAU_pred_yi_hi,
                     xmin = MAU_pred_yi_lo,
                     y = seq_len(n_studies),
                     color = "MA-U")) +
  
  geom_errorbarh(data = dr_pred_interval,
                 aes(xmax = MAC_pred_yi_hi,
                     xmin = MAC_pred_yi_lo,
                     y = seq_len(n_studies),
                     color = "MA-C")) +
  
  geom_point(data = data.frame(x = NA_real_, y = NA_real_),
             aes(x = x, y = y, color = "MA-U")) +
  geom_point(data = data.frame(x = NA_real_, y = NA_real_),
             aes(x = x, y = y, color = "MA-C")) +
  
  scale_y_continuous(breaks = seq_len(n_studies), labels = study_labels) +
  scale_color_manual(
    name   = NULL,
    values = c("REP effect size" = "black",
               "MA-U" = "#D55E00",
               "MA-C" = "#0072B2"),
    breaks = c("REP effect size", "MA-U", "MA-C"),
    guide  = legend_guide
  ) +
  shared_theme

ggsave("outputs/MB4_PI_MA.jpeg", p_ma, height = 12, width = 10, dpi = 500)


########### Plot 2: MR-U and MR-C ###########
p_mr <- ggplot() +
  geom_point(data = dr_pred_est,
             aes(x = yi, y = seq_len(n_studies),
                 color = "REP effect size"),
             size = 1) +
  
  geom_point(data = dr_pred_est,
             aes(x = pred.MRU, y = seq_len(n_studies),
                 color = "MR-U"),
             shape = 17, size = 1.5) +
  
  geom_point(data = dr_pred_est,
             aes(x = pred.MRC, y = seq_len(n_studies),
                 color = "MR-C"),
             shape = 15, size = 1.5) +
  
  geom_errorbarh(data = dr_pred_interval,
                 aes(xmax = MRU_pred_yi_hi,
                     xmin = MRU_pred_yi_lo,
                     y = seq_len(n_studies),
                     color = "MR-U")) +
  
  geom_errorbarh(data = dr_pred_interval,
                 aes(xmax = MRC_pred_yi_hi,
                     xmin = MRC_pred_yi_lo,
                     y = seq_len(n_studies),
                     color = "MR-C")) +
  
  scale_y_continuous(breaks = seq_len(n_studies), labels = study_labels) +
  scale_color_manual(
    name   = NULL,
    values = c("REP effect size" = "black",
               "MR-U" = "#E69F00",
               "MR-C" = "#009E73"),
    breaks = c("REP effect size", "MR-U", "MR-C"),
    guide  = legend_guide
  ) +
  shared_theme

ggsave("outputs/MB4_PI_MR.jpeg", p_mr, height = 12, width = 10, dpi = 500)
