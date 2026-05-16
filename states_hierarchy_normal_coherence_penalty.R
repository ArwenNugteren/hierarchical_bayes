library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

# ── Data prep ─────────────────────────────────────────────────────────────────

load("tourism.Rdata")

# times
training_end <- 96
m            <- 12L
horizon      <- 12

# State series A-G (single letter codes)
state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)

# Number of bottom-level series per state (for consistent scaling)
# Count series matching "^X" where X is the state letter
series_cols <- setdiff(names(tourism), c("Date", "Total",
                                         names(tourism)[grepl("^Total", names(tourism))]))

regions_per_state <- sapply(state_codes, function(s) {
  sum(grepl(paste0("^", s, "[A-Z]{2}(Hol|Vis|Bus|Oth)$"), names(tourism)))
})
regions_per_total <- sum(regions_per_state)

Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)


# Scale from bottom-level region data, consistent with earlier approach
# Use the state-level training data to define scale
Y_mean <- mean(Y_raw_state[1:training_end, ])
Y_sd   <- sd(Y_raw_state[1:training_end, ])

# States: standard scaling
Y_state <- (Y_raw_state - Y_mean) / Y_sd

# Total: scaled by n_states equivalent (sum of states)
Y_total_scaled <- (Y_raw_total - Y_mean * n_state) / Y_sd

lambda <- 1
sigma_coh <- 0.1  # variance

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
  real<lower=0> sigma_coh;
}
parameters {
  // ── Hyperpriors on alpha and gamma ────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;

  // ── State-level smoothing parameters ─────────────────────────────────────
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_state_raw;
  vector<lower=0>[n_state]          sigma_state;

  // ── Total-level smoothing parameters (drawn from same hyperprior) ─────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_total_raw;
  real<lower=0>          sigma_total;

  // ── Initial values  ──────────────────────────────────────────
  real           l0_total;
  vector[n_state] l0_state;

  matrix[n_state, m-1] s0_state_free;
  vector[m-1]          s0_total_free;
}
transformed parameters {
  // constraining gamma to be smaller than 1 - alpha
  vector[n_state] gamma_state;
  real gamma_total;
  for (n in 1:n_state){
    gamma_state[n] = gamma_state_raw[n] * (1 - alpha_state[n]);
  }
  gamma_total = gamma_total_raw * (1-alpha_total);
  
  // initial value sum-to-zero constraint
  matrix[n_state, m] s0_state;
  vector[m]          s0_total;
  
  for (n in 1:n_state) {
    s0_state[n, 1:(m-1)] = s0_state_free[n];
    s0_state[n, m]        = -sum(s0_state_free[n]);
  }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);
  
  
  // ── State recursions ──────────────────────────────────────────────────────
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

  // ── Total recursion ───────────────────────────────────────────────────────
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
  // ── Hyperpriors ───────────────────────────────────────────────────────────
  mu_alpha ~ uniform(0,1);
  mu_gamma_raw ~ uniform(0,1);

  // ── State-level priors drawn from hyperprior ──────────────────────────────
  alpha_state ~ normal(mu_alpha, 0.1);
  gamma_state_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_state ~ cauchy(0, 0.3);

  // ── Total-level priors drawn from same hyperprior ─────────────────────────
  alpha_total ~ normal(mu_alpha, 0.1);
  gamma_total_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_total ~ cauchy(0, 0.3);

  // ── Initial values ──────────────────────────────────────────────────
  l0_state ~ normal(0, 1);
  l0_total ~ normal(0, n_state);
  
  to_vector(s0_state_free) ~ normal(0, 1); 
  s0_total_free ~ normal(0, n_state);
  
  // -- constrain initial values
  for (k in 1:m){
    real diff = (l0_total + s0_total[k]) - (sum(l0_state) + sum(s0_state[,k]));
    target += normal_lpdf(diff | 0, sigma_coh);
  }

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_state) {
    for (i in 1:T_train){
      target += normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
    }
  }
  for (i in 1:T_train) {
    target += lambda * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
  }
    
  // coherence penalty
  for (i in 1:T_train) {
    real fitted_total = l_t[i] + s_t[i];
    real fitted_sum   = 0;
    for (n in 1:n_state)
      fitted_sum += l_s[n][i] + s_s[n][i];
  // add penalty based on normal
    target += normal_lpdf(fitted_total - fitted_sum | 0, sigma_coh);
  }
}
generated quantities {
  matrix[H, n_state] Y_fc_state;
  vector[H]          Y_fc_total;

  // ── Diagnostics ───────────────────────────────────────────────────────────
  vector[T_train] lpdf_state_sum;   // sum across states per timestep
  vector[T_train] lpdf_total;       // total-level log likelihood per timestep
  vector[T_train] coherence_gap;    // fitted_total - fitted_sum
  vector[T_train] lpdf_coherence;   // coherence penalty per timestep
  vector[T_train] target_contrib;   // total target += per timestep

  for (i in 1:T_train) {
    // State log likelihoods summed across states
    real state_ll = 0;
    for (n in 1:n_state)
      state_ll += normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
    lpdf_state_sum[i] = state_ll;

    // Total log likelihood
    lpdf_total[i] = lambda * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);

    // Coherence gap and penalty
    real fitted_total = l_t[i] + s_t[i];
    real fitted_sum   = 0;
    for (n in 1:n_state)
      fitted_sum += l_s[n][i] + s_s[n][i];
    coherence_gap[i]    = fitted_total - fitted_sum;
    lpdf_coherence[i]   = normal_lpdf(coherence_gap[i] | 0, sigma_coh);

    // Total contribution at this timestep
    target_contrib[i] = lpdf_state_sum[i] + lpdf_total[i] + lpdf_coherence[i];
  }

  // ── Forecasts (unchanged) ─────────────────────────────────────────────────
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

# ── Data list ─────────────────────────────────────────────────────────────────

data_list <- list(
  T_train  = training_end,
  m        = m,
  n_state  = n_state,
  H        = horizon + 1,
  Y_state  = Y_state[1:training_end, ],
  Y_total  = Y_total_scaled[1:training_end], 
  lambda = lambda, 
  sigma_coh = sigma_coh
)

stanDso <- stan_model(model_code = model_string)

save_seed <- 2025
start_time <- Sys.time()
stanFit <- sampling(
  object = stanDso,
  data   = data_list,
  seed   = save_seed,
  chains = 10,
  iter   = 2000,
  warmup = 1000,
  cores  = 10
  # control = list(adapt_delta = 0.95)
)
end_time <- Sys.time()
print(end_time - start_time)

# Extract posterior means across draws for each timestep
diag_vars <- c("lpdf_state_sum", "lpdf_total", "coherence_gap", 
               "lpdf_coherence", "target_contrib")

diag_summary <- lapply(diag_vars, function(v) {
  draws <- rstan::extract(stanFit, pars = v)[[v]]  # [draws x T_train]
  data.frame(
    timestep = 1:training_end,
    variable = v,
    mean     = colMeans(draws),
    q05      = apply(draws, 2, quantile, 0.05),
    q95      = apply(draws, 2, quantile, 0.95)
  )
}) %>% bind_rows()

# View as wide table for easy comparison
diag_wide <- diag_summary %>%
  select(timestep, variable, mean) %>%
  tidyr::pivot_wider(names_from = variable, values_from = mean)

View(diag_wide)


# cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
# cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")
# 
# 
# View(summary(stanFit)$summary)

all_pars <- names(stanFit)
coh_pars <- all_pars[grepl("sigma_coh", all_pars)]
sigma_coh_mean <- summary(stanFit)$summary[coh_pars, "mean"]
print(sigma_coh_mean)
# View(summary(stanFit_VB)$summary)

all_pars <- names(stanFit)
sigma_pars <- all_pars[grepl("sigma", all_pars)]
sigma_means <- summary(stanFit)$summary[sigma_pars, "mean"]
sigma_means

# # hyperparameters
# pairs(stanFit, pars = c("mu_alpha", "mu_gamma_raw"))
# 
# # lower level parameters
# # get all parameter names from the fit
# all_pars <- names(stanFit)
# 
# # grab all alpha-related parameters
# alpha_pars <- all_pars[grepl("alpha_state", all_pars)]
# pairs(stanFit, pars = alpha_pars)
# 
# alpha_means <- summary(stanFit)$summary[alpha_pars, "mean"]
# print(alpha_means)
# 
# # grab all gamma-related parameters
# gamma_pars <- all_pars[grepl("gamma_state", all_pars)]
# pairs(stanFit, pars = gamma_pars)
# 
# gamma_means <- summary(stanFit)$summary[gamma_pars, "mean"]
# print(gamma_means)
# 
# # grab all l0-related parameters
# l0_pars <- all_pars[grepl("l0", all_pars)]
# 
# pairs(stanFit, pars = l0_pars)




# ── accuracy calculation ─────────────────────────────────────────────────
fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

fc_state_orig  <- fc_state_draws[, 1:horizon, ] * Y_sd + Y_mean
fc_total_orig  <- fc_total_draws[, 1:horizon]   * Y_sd + Y_mean * n_state

actuals_state  <- Y_raw_state[(training_end + 1):(training_end + horizon), ]  # [horizon x n_state]
actuals_total <- Y_raw_total[(training_end + 1):(training_end + horizon)]    # [horizon]


# ── Compute per-horizon accuracy ───────────────────────────────────────────────

horizon_results <- lapply(1:horizon, function(h) {
  
  series_results <- lapply(1:n_state, function(s) {
    draws  <- fc_state_orig[, h, s]
    actual <- actuals_state[h, s]
    q05    <- quantile(draws, 0.05)
    q95    <- quantile(draws, 0.95)
    data.frame(
      origin  = training_end,
      horizon = h,
      series  = state_codes[s],
      actual  = actual,
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
  q05_t    <- quantile(draws_t, 0.05)
  q95_t    <- quantile(draws_t, 0.95)
  
  bind_rows(
    series_results,
    data.frame(
      origin  = training_end,
      horizon = h,
      series  = "Total",
      actual  = actual_t,
      fc_mean = mean(draws_t),
      ae      = abs(mean(draws_t) - actual_t),
      se      = (mean(draws_t) - actual_t)^2,
      crps    = scoringRules::crps_sample(actual_t, draws_t),
      logS    = scoringRules::logs_sample(actual_t, draws_t),
      covered = as.integer(actual_t >= q05_t & actual_t <= q95_t)
    )
  )
})

accuracy <- bind_rows(horizon_results) %>%
  group_by(series) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    LogS     = mean(logS),
    Coverage = mean(covered),
    .groups  = "drop")
# ── Coherence accuracy ─────────────────────────────────────────────────────────
# fc_state_orig: [draws x horizon x n_state], fc_total_orig: [draws x horizon]

coherence_results <- lapply(1:horizon, function(h) {
  total_mean     <- mean(fc_total_orig[, h])
  state_sum_mean <- sum(apply(fc_state_orig[, h, ], 2, mean))  # mean per state, then sum
  gap            <- total_mean - state_sum_mean
  
  data.frame(
    horizon        = h,
    MAE_coherence  = abs(gap),
    RMSE_coherence = sqrt(gap^2)
  )
})

coherence_by_horizon <- bind_rows(coherence_results)

coherence_overall <- coherence_by_horizon %>%
  summarise(
    MAE_coherence  = mean(MAE_coherence),
    RMSE_coherence = sqrt(mean(RMSE_coherence^2))  # correct aggregation for RMSE
  )

# ── Forecast vs actual means by horizon ───────────────────────────────────────

fc_vs_actual <- bind_rows(horizon_results) %>%
  group_by(horizon, series) %>%
  summarise(
    fc_mean = mean(fc_mean),
    actual  = mean(actual),
    .groups = "drop"
  ) %>%
  mutate(gap = fc_mean - actual)
View(fc_vs_actual)

print(accuracy)
print(coherence_overall)

