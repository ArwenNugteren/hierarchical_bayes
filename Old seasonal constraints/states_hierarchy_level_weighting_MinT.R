library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

# ── Data prep ─────────────────────────────────────────────────────────────────

load("tourism.Rdata")

training_end <- 96
m            <- 12L
horizon      <- 12

state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)

series_cols <- setdiff(names(tourism), c("Date", "Total",
                                         names(tourism)[grepl("^Total", names(tourism))]))

Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

Y_mean <- mean(Y_raw_state[1:training_end, ])
Y_sd   <- sd(Y_raw_state[1:training_end, ])

Y_state        <- (Y_raw_state - Y_mean) / Y_sd
Y_total_scaled <- (Y_raw_total - Y_mean * n_state) / (Y_sd * n_state)

lambda <- 7

# ── Stan model ────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> H;
  matrix[T_train, n_state] Y_state;
  vector[T_train]           Y_total;
  real<lower=0> lambda;
}
parameters {
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_state_raw;
  vector<lower=0>[n_state]          sigma_state;
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_total_raw;
  real<lower=0>          sigma_total;
  vector[n_state]    l0_state;
  matrix[n_state, m] s0_state;
  real      l0_total;
  vector[m] s0_total;
}
transformed parameters {
  vector[n_state] gamma_state;
  real gamma_total;
  for (n in 1:n_state)
    gamma_state[n] = gamma_state_raw[n] * (1 - alpha_state[n]);
  gamma_total = gamma_total_raw * (1 - alpha_total);

  array[n_state] vector[T_train + 1] l_s;
  array[n_state] vector[T_train + m] s_s;
  for (n in 1:n_state) {
    for (j in 1:m) s_s[n][j] = s0_state[n, j];
    l_s[n][1] = l0_state[n];
    for (i in 1:T_train) {
      real s_lag   = s_s[n][i];
      l_s[n][i+1] = alpha_state[n] * (Y_state[i, n] - s_lag)
                  + (1 - alpha_state[n]) * l_s[n][i];
      s_s[n][i+m] = gamma_state[n] * (Y_state[i, n] - l_s[n][i])
                  + (1 - gamma_state[n]) * s_lag;
    }
  }

  vector[T_train + 1] l_t;
  vector[T_train + m] s_t;
  for (j in 1:m) s_t[j] = s0_total[j];
  l_t[1] = l0_total;
  for (i in 1:T_train) {
    real s_lag = s_t[i];
    l_t[i+1]  = alpha_total * (Y_total[i] - s_lag)
              + (1 - alpha_total) * l_t[i];
    s_t[i+m]  = gamma_total * (Y_total[i] - l_t[i])
              + (1 - gamma_total) * s_lag;
  }
}
model {
  mu_alpha     ~ uniform(0, 1);
  mu_gamma_raw ~ uniform(0, 1);
  alpha_state     ~ normal(mu_alpha, 0.1);
  gamma_state_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_state     ~ cauchy(0, 0.3);
  alpha_total     ~ normal(mu_alpha, 0.1);
  gamma_total_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 0.3);
  l0_state ~ normal(0, 1);
  l0_total ~ normal(0, 1);
  for (n in 1:n_state) {
    s0_state[n] ~ normal(0, 1);
    sum(s0_state[n]) ~ normal(0, 0.01 * m);
  }
  s0_total ~ normal(0, 1);
  sum(s0_total) ~ normal(0, 0.01 * m);
  for (n in 1:n_state)
    for (i in 1:T_train)
      target += normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
  for (i in 1:T_train)
    target += lambda * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
}
generated quantities {
  matrix[H, n_state] Y_fc_state;
  vector[H]          Y_fc_total;

  for (n in 1:n_state) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_s[n][T_train + 1];
    for (h in 1:H) sf[h] = s_s[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_state[h, n] = normal_rng(lf[h] + sf[h], sigma_state[n]);
      lf[h+1] = alpha_state[n] * (Y_fc_state[h, n] - sf[h])
              + (1 - alpha_state[n]) * lf[h];
    }
    Y_fc_state[H, n] = normal_rng(lf[H] + sf[H], sigma_state[n]);
  }
  {
    vector[H] lf; vector[H] sf;
    lf[1] = l_t[T_train + 1];
    for (h in 1:H) sf[h] = s_t[T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_total[h] = normal_rng(lf[h] + sf[h], sigma_total);
      lf[h+1] = alpha_total * (Y_fc_total[h] - sf[h])
              + (1 - alpha_total) * lf[h];
    }
    Y_fc_total[H] = normal_rng(lf[H] + sf[H], sigma_total);
  }
}
"

data_list <- list(
  T_train = training_end,
  m       = m,
  n_state = n_state,
  H       = horizon + 1,
  Y_state = Y_state[1:training_end, ],
  Y_total = Y_total_scaled[1:training_end],
  lambda  = lambda
)

stanDso <- stan_model(model_code = model_string)

save_seed  <- 2025
start_time <- Sys.time()
stanFit <- sampling(
  object = stanDso,
  data   = data_list,
  seed   = save_seed,
  chains = 10,
  iter   = 2000,
  warmup = 1000,
  cores  = 10
)
end_time <- Sys.time()
print(end_time - start_time)

cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")






# ── Extract and back-transform draws ──────────────────────────────────────────

fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

# Back-transform to original units
# fc_state_orig: [draws x horizon x n_state]
# fc_total_orig: [draws x horizon]
fc_state_orig <- fc_state_draws[, 1:horizon, ] * Y_sd + Y_mean
fc_total_orig <- fc_total_draws[, 1:horizon]   * Y_sd * n_state + Y_mean * n_state

actuals_state <- Y_raw_state[(training_end + 1):(training_end + horizon), ]
actuals_total <- Y_raw_total[(training_end + 1):(training_end + horizon)]

# ── MinT reconciliation ───────────────────────────────────────────────────────
# Hierarchy: Total (1 series) aggregates n_state bottom series
# Summing matrix S: [n_state+1 x n_state]
#   top row: all ones (Total = sum of states)
#   bottom:  identity (states map to themselves)

n_all <- n_state + 1L   # total number of series in hierarchy
S <- rbind(rep(1, n_state), diag(n_state))   # [(n_state+1) x n_state]

mint_reconcile <- function(base_means, W) {
  # base_means: named vector of length n_all, order: Total, state1..stateN
  # W:          covariance matrix of base forecast errors [n_all x n_all]
  # Returns reconciled means for all series
  P <- solve(t(S) %*% solve(W) %*% S) %*% t(S) %*% solve(W)  # [n_state x n_all]
  G <- S %*% P                                                  # [n_all x n_all]
  as.vector(G %*% base_means)
}

# Estimate W (shrinkage estimator) from in-sample residuals across draws
# Use posterior mean parameters to get a single residual matrix for W estimation
# Residuals: [T_train x n_all] — total first, then states
resid_total  <- Y_raw_total[1:training_end] -
  (colMeans(rstan::extract(stanFit, "l_t")$l_t) +
     colMeans(rstan::extract(stanFit, "s_t")$s_t))[1:training_end] * Y_sd * n_state +
  Y_mean * n_state

# Build residual matrix from posterior means of fitted values
l_s_mean <- apply(rstan::extract(stanFit, "l_s")$l_s, c(2,3), mean)  # [n_state x T+1]
s_s_mean <- apply(rstan::extract(stanFit, "s_s")$s_s, c(2,3), mean)  # [n_state x T+m]
l_t_mean <- colMeans(rstan::extract(stanFit, "l_t")$l_t)              # [T+1]
s_t_mean <- colMeans(rstan::extract(stanFit, "s_t")$s_t)              # [T+m]

fitted_states <- t(l_s_mean[, 1:training_end] + s_s_mean[, 1:training_end]) * Y_sd + Y_mean
fitted_total  <- (l_t_mean[1:training_end] + s_t_mean[1:training_end]) * Y_sd * n_state +
  Y_mean * n_state

resid_mat <- cbind(
  Y_raw_total[1:training_end] - fitted_total,
  Y_raw_state[1:training_end, ] - fitted_states
)  # [T_train x n_all], columns: Total, state1..stateN

# Ledoit-Wolf shrinkage covariance (diagonal target)
lw_shrink <- function(E) {
  T  <- nrow(E)
  n  <- ncol(E)
  S  <- cov(E)
  mu <- mean(diag(S))
  # shrinkage intensity via Ledoit-Wolf formula
  num <- sum(apply(E, 1, function(e) sum((tcrossprod(e) - S)^2))) / T^2
  den <- sum((S - diag(mu, n))^2)
  rho <- min(1, num / den)
  (1 - rho) * S + rho * diag(mu, n)
}

W_hat <- lw_shrink(resid_mat)   # [n_all x n_all]

# ── Apply MinT per horizon using posterior mean forecasts and base variances ─────────────────────
base_sds <- bind_rows(
  data.frame(
    series  = "Total",
    horizon = 1:horizon,
    fc_sd   = apply(fc_total_orig, 2, sd)
  ),
  do.call(rbind, lapply(1:n_state, function(s)
    data.frame(
      series  = state_codes[s],
      horizon = 1:horizon,
      fc_sd   = apply(fc_state_orig[, , s], 2, sd)
    )
  ))
)

mint_results <- lapply(1:horizon, function(h) {
  base_means <- c(mean(fc_total_orig[, h]),
                  colMeans(fc_state_orig[, h, ]))
  reconciled <- mint_reconcile(base_means, W_hat)
  data.frame(
    horizon      = h,
    series       = c("Total", state_codes),
    fc_mean_base = base_means,
    fc_mean_mint = reconciled
  )
}) %>%
  bind_rows() %>%
  left_join(base_sds, by = c("series", "horizon"))

# ── Accuracy: base forecasts ───────────────────────────────────────────────────

horizon_results_base <- lapply(1:horizon, function(h) {
  series_results <- lapply(1:n_state, function(s) {
    draws  <- fc_state_orig[, h, s]
    actual <- actuals_state[h, s]
    q05    <- quantile(draws, 0.05)
    q95    <- quantile(draws, 0.95)
    data.frame(
      horizon = h, series = state_codes[s], actual = actual,
      fc_mean = mean(draws),
      ae      = abs(mean(draws) - actual),
      se      = (mean(draws) - actual)^2,
      crps    = scoringRules::crps_sample(actual, draws),
      logS    = scoringRules::logs_sample(actual, draws),
      covered = as.integer(actual >= q05 & actual <= q95)
    )
  })
  draws_t  <- fc_total_orig[, h]
  actual_t <- actuals_total[h]
  bind_rows(series_results,
            data.frame(
              horizon = h, series = "Total", actual = actual_t,
              fc_mean = mean(draws_t),
              ae      = abs(mean(draws_t) - actual_t),
              se      = (mean(draws_t) - actual_t)^2,
              crps    = scoringRules::crps_sample(actual_t, draws_t),
              logS    = scoringRules::logs_sample(actual_t, draws_t),
              covered = as.integer(actual_t >= quantile(draws_t, 0.05) &
                                     actual_t <= quantile(draws_t, 0.95))
            ))
}) %>% bind_rows()

# ── Accuracy: MinT forecasts ───────────────────────────────────────────────────
# Note: MinT operates on posterior means only, so no draws available for
# CRPS/LogS/Coverage — point forecast metrics only

actuals_df <- bind_rows(
  data.frame(series = "Total",    horizon = 1:horizon, actual = actuals_total),
  do.call(rbind, lapply(1:n_state, function(s)
    data.frame(series  = state_codes[s],
               horizon = 1:horizon,
               actual  = actuals_state[, s])
  ))
)

horizon_results_mint <- mint_results %>%
  left_join(actuals_df, by = c("series", "horizon")) %>%
  mutate(
    ae      = abs(fc_mean_mint - actual),
    se      = (fc_mean_mint - actual)^2,
    # Shifted normal: same SD as base, new reconciled mean
    crps    = scoringRules::crps_norm(actual, mean = fc_mean_mint, sd = fc_sd),
    logS    = scoringRules::logs_norm(actual, mean = fc_mean_mint, sd = fc_sd),
    covered = as.integer(actual >= qnorm(0.05, fc_mean_mint, fc_sd) &
                           actual <= qnorm(0.95, fc_mean_mint, fc_sd))
  )

accuracy_mint <- horizon_results_mint %>%
  group_by(series) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    LogS     = mean(logS),
    Coverage = mean(covered),
    .groups  = "drop"
  ) %>%
  bind_rows(horizon_results_mint %>%
              summarise(
                series   = "Overall",
                MAE      = mean(ae),
                RMSE     = sqrt(mean(se)),
                CRPS     = mean(crps),
                LogS     = mean(logS),
                Coverage = mean(covered)
              ))

# ── Summary tables ─────────────────────────────────────────────────────────────

accuracy_base <- horizon_results_base %>%
  group_by(series) %>%
  summarise(MAE = mean(ae), RMSE = sqrt(mean(se)),
            CRPS = mean(crps), LogS = mean(logS),
            Coverage = mean(covered), .groups = "drop") %>%
  bind_rows(horizon_results_base %>%
              summarise(series = "Overall", MAE = mean(ae), RMSE = sqrt(mean(se)),
                        CRPS = mean(crps), LogS = mean(logS), Coverage = mean(covered)))

# ── Coherence ─────────────────────────────────────────────────────────────────

coherence_base <- lapply(1:horizon, function(h) {
  gap <- mean(fc_total_orig[, h]) - sum(colMeans(fc_state_orig[, h, ]))
  data.frame(horizon = h, MAE_coherence = abs(gap), RMSE_coherence = sqrt(gap^2))
}) %>% bind_rows()

coherence_mint <- mint_results %>%
  group_by(horizon) %>%
  summarise(
    gap = fc_mean_mint[series == "Total"] - sum(fc_mean_mint[series != "Total"]),
    .groups = "drop"
  ) %>%
  mutate(MAE_coherence = abs(gap), RMSE_coherence = sqrt(gap^2))

print(accuracy_base)
print(accuracy_mint)
print(summarise(coherence_base,
        MAE_coherence = mean(MAE_coherence), RMSE_coherence = sqrt(mean(RMSE_coherence^2))))
print(summarise(coherence_mint,
        MAE_coherence = mean(MAE_coherence), RMSE_coherence = sqrt(mean(RMSE_coherence^2))))
