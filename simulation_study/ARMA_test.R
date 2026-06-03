# =============================================================================
# ARMA(2,2) — single train/test split, horizon = 5, draws = 2000
# Matches the tourism-style accuracy framework
# =============================================================================

library(forecast)
library(dplyr)
library(scoringRules)

# -----------------------------------------------------------------------------
# 0.  CONFIG
# -----------------------------------------------------------------------------
TRAIN_N  <- 100L
HORIZON  <- 5L
N_DRAWS  <- 2000L
set.seed(2025)

# -----------------------------------------------------------------------------
# 1.  LOAD DATA
# -----------------------------------------------------------------------------
synthetic_data <- read.csv("synthetic_hierarchy.csv")

# Expects columns: Total, A, B, AA, AB, BA, BB
stopifnot(all(c("Total","A","B","AA","AB","BA","BB") %in% names(synthetic_data)))

series_names <- c("Total", "A", "B", "AA", "AB", "BA", "BB")

train <- synthetic_data[1:TRAIN_N, series_names]
test  <- synthetic_data[(TRAIN_N + 1):(TRAIN_N + HORIZON), series_names]

# Index helpers (matching tourism IDX pattern)
IDX <- list(
  total = "Total",
  A     = "A",
  B     = "B",
  bot   = c("AA", "AB", "BA", "BB")
)

# -----------------------------------------------------------------------------
# 2.  FIT ARMA(2,2) AND GENERATE DRAWS  [draws x horizon]
# -----------------------------------------------------------------------------
fit_and_draw <- function(y_train, h = HORIZON, n_draws = N_DRAWS) {
  fit    <- Arima(y_train, order = c(2, 0, 2))
  sim_mat <- matrix(NA_real_, nrow = n_draws, ncol = h)
  for (d in seq_len(n_draws)) {
    sim_mat[d, ] <- simulate(fit, nsim = h, future = TRUE)
  }
  sim_mat   # [n_draws x h]
}

cat("Fitting ARMA(2,2) models...\n")
fc <- lapply(series_names, function(s) {
  cat("  ", s, "\n")
  fit_and_draw(train[[s]])
})
names(fc) <- series_names

# Convenience handles matching tourism variable names
fc_total <- fc[["Total"]]                          # [draws x h]
fc_A     <- fc[["A"]]
fc_B     <- fc[["B"]]
fc_bot   <- array(                                 # [draws x h x 4]
  data = c(fc[["AA"]], fc[["AB"]], fc[["BA"]], fc[["BB"]]),
  dim  = c(N_DRAWS, HORIZON, 4)
)

# -----------------------------------------------------------------------------
# 3.  ACCURACY HELPERS  (verbatim from tourism script)
# -----------------------------------------------------------------------------
JSD_metric <- function(p_draws, q_draws, n_grid = 1000) {
  lo    <- min(c(p_draws, q_draws)) - 3 * sd(c(p_draws, q_draws))
  hi    <- max(c(p_draws, q_draws)) + 3 * sd(c(p_draws, q_draws))
  p_kde <- pmax(density(p_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  q_kde <- pmax(density(q_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  m_kde <- 0.5 * (p_kde + q_kde)
  dx    <- (hi - lo) / (n_grid - 1)
  jsd   <- 0.5 * (sum(p_kde * log(p_kde / m_kde)) +
                    sum(q_kde * log(q_kde / m_kde))) * dx
  sqrt(max(0, min(log(2), jsd)))
}

compute_energy_score <- function(draws_mat, actuals, sub_idx = NULL) {
  if (!is.null(sub_idx)) draws_mat <- draws_mat[sub_idx, , drop = FALSE]
  col_sds    <- apply(draws_mat, 2, sd)
  draws_sc   <- sweep(draws_mat, 2, col_sds, "/")
  actuals_sc <- as.numeric(actuals) / col_sds
  diffs1     <- sweep(draws_sc, 2, actuals_sc, "-")
  t1         <- mean(sqrt(rowSums(diffs1^2)))
  t2         <- mean(dist(draws_sc, method = "euclidean"))
  t1 - 0.5 * t2
}

n_es_subsample <- 1000L
sub_idx <- sample.int(N_DRAWS, size = min(n_es_subsample, N_DRAWS), replace = FALSE)

bot_names <- c("AA", "AB", "BA", "BB")

make_row <- function(draws, actual, series, level) {
  data.frame(
    series     = series,
    level      = level,
    actual     = actual,
    fc_mean    = mean(draws),
    ae         = abs(mean(draws) - actual),
    se         = (mean(draws) - actual)^2,
    crps       = scoringRules::crps_sample(actual, draws),
    logS       = scoringRules::logs_sample(actual, draws),
    covered_90 = as.integer(actual >= quantile(draws, 0.05) &
                              actual <= quantile(draws, 0.95)),
    covered_95 = as.integer(actual >= quantile(draws, 0.025) &
                              actual <= quantile(draws, 0.975)),
    pit        = mean(draws < actual)
  )
}

# -----------------------------------------------------------------------------
# 4.  EVALUATE OVER HORIZONS
# -----------------------------------------------------------------------------
roll <- 1L   # single split — roll = 1 throughout

horizon_results <- lapply(seq_len(HORIZON), function(h) {
  
  actual_bot   <- as.numeric(test[h, IDX$bot])
  actual_A     <- as.numeric(test[h, IDX$A])
  actual_B     <- as.numeric(test[h, IDX$B])
  actual_total <- as.numeric(test[h, IDX$total])
  
  univariate_df <- bind_rows(
    lapply(1:4, function(j)
      make_row(fc_bot[, h, j], actual_bot[j], bot_names[j], "bottom")),
    make_row(fc_A[, h],     actual_A,     "A",     "mid"),
    make_row(fc_B[, h],     actual_B,     "B",     "mid"),
    make_row(fc_total[, h], actual_total, "Total", "top")
  ) %>% mutate(roll = roll, horizon = h, .before = series)
  
  # Multivariate coherence metrics
  bot_mat      <- fc_bot[, h, ]
  sum_bot      <- rowSums(bot_mat)
  sum_mid      <- fc_A[, h] + fc_B[, h]
  draws_total  <- fc_total[, h]
  
  joint_mat     <- cbind(bot_mat, fc_A[, h], fc_B[, h], draws_total)
  actuals_joint <- c(actual_bot, actual_A, actual_B, actual_total)
  es_val        <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
  
  multivariate_df <- data.frame(
    roll             = roll,
    horizon          = h,
    es               = es_val,
    mai_bot_mid      = mean(sum_bot)  - mean(sum_mid),
    mai_mid_top      = mean(sum_mid)  - mean(draws_total),
    mai_bot_top      = mean(sum_bot)  - mean(draws_total),
    jsd_bot_mid      = JSD_metric(sum_bot,  sum_mid),
    jsd_mid_top      = JSD_metric(sum_mid,  draws_total),
    jsd_bot_top      = JSD_metric(sum_bot,  draws_total)
  )
  
  list(univariate = univariate_df, multivariate = multivariate_df)
})

# -----------------------------------------------------------------------------
# 5.  COLLECT AND DISPLAY
# -----------------------------------------------------------------------------
roll_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
roll_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))

cat("\n===== UNIVARIATE RESULTS =====\n")
print(roll_univariate, digits = 4)

cat("\n===== MULTIVARIATE / COHERENCE RESULTS =====\n")
print(roll_multivariate, digits = 4)
# 
# # Optional: write to CSV
# write.csv(roll_univariate,   "results_univariate.csv",   row.names = FALSE)
# write.csv(roll_multivariate, "results_multivariate.csv", row.names = FALSE)
# cat("\nResults written to results_univariate.csv and results_multivariate.csv\n")