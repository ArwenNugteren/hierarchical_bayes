library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

load("tourism.Rdata")

# read job name and directories
args     <- commandArgs(trailingOnly = TRUE)
job_name <- args[1]

dir_accuracy <- file.path(job_name, "accuracy")
dir_stanfits <- file.path(job_name, "stanfits")

# ── Setup ──────────────────────────────────────────────────────────────────────

state_codes  <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state      <- length(state_codes)
Y_raw_state  <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_total  <- tourism %>% pull(Total)

m           <- 12L
window_size <- 96
horizon     <- 12
test_set    <- 121
origins     <- 96:(96+test_set-1)   # 121 windows


# ── Get this task's origin from PBS array index ────────────────────────────────

task_id <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
origin  <- origins[task_id]

# ── Training window ────────────────────────────────────────────────────────────

train_idx <- (origin - window_size + 1):origin
test_idx  <- (origin + 1):(origin + horizon)

Y_mean_w  <- mean(Y_raw_state[train_idx, ])
Y_sd_w    <- sd(Y_raw_state[train_idx, ])

Y_state_w <- (Y_raw_state[train_idx, ] - Y_mean_w) / Y_sd_w
Y_total_w <- (Y_raw_total[train_idx]   - Y_mean_w * n_state) / (Y_sd_w * n_state)

# hierarchical weighting
lambda <- 1

# coherence penalty weighting
kappa <- 1 

# ── Stan model ─────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> H;
  matrix[T_train, n_state] Y_state;
  vector[T_train]           Y_total;
  
  real<lower=0> lambda; 
  real<lower=0> kappa;
}
parameters {
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma;
  //real<lower=0> sigma_alpha;
  //real<lower=0> sigma_gamma;
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_state;
  vector<lower=0>[n_state]          sigma_state;
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_total;
  real<lower=0>          sigma_total;
  vector[n_state] l0_state;
  matrix[n_state, m] s0_state;
  real           l0_total;
  vector[m]      s0_total;
}
transformed parameters {
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
  mu_alpha    ~ uniform(0, 1);
  mu_gamma    ~ uniform(0, 1);
  //sigma_alpha ~ normal(0, 0.3);
  //sigma_gamma ~ normal(0, 0.3);
  alpha_state ~ normal(mu_alpha, 0.1);
  gamma_state ~ normal(mu_gamma, 0.1);
  sigma_state ~ cauchy(0, 1);
  alpha_total ~ normal(mu_alpha, 0.1);
  gamma_total ~ normal(mu_gamma, 0.1);
  sigma_total ~ cauchy(0, 1);
  l0_state ~ normal(0, 1);
  l0_total ~ normal(0, 1);
  for (n in 1:n_state) {
    s0_state[n] ~ normal(0, 1);
    sum(s0_state[n]) ~ normal(0, 0.01 * m);
  }
  s0_total ~ normal(0, 1);
  sum(s0_total) ~ normal(0, 0.01 * m);
  
  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_state) {
    for (i in 1:T_train){
      //Y_state[i, n] ~ normal(l_s[n][i] + s_s[n][i], sigma_state[n]);
      target += normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
    }
  }
  for (i in 1:T_train){
    target += lambda * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
  }
  
  // coherence penalty
  for (i in 1:T_train) {
    real fitted_total = l_t[i] + s_t[i];
    real fitted_sum   = 0;
    for (n in 1:n_state)
      fitted_sum += l_s[n][i] + s_s[n][i];
    // add penalty based on the square of the normalised difference
    target += -kappa * square((fitted_total - fitted_sum)/sigma_total);
  }
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

stanDso <- stan_model(model_code = model_string)

stanFit <- sampling(
  object  = stanDso,
  data    = list(
    T_train = window_size,
    m       = m,
    n_state = n_state,
    H       = horizon + 1,
    Y_state = Y_state_w,
    Y_total = Y_total_w, 
    lambda = lambda, 
    kappa = kappa
  ),
  seed    = 2025,
  chains  = 10,
  iter    = 2000,
  warmup  = 1000,
  cores   = 10,
  refresh = 0
)

# save the stanFits (in case I need to check them)
saveRDS(stanFit,   file.path(dir_stanfits, sprintf("stanfit_timestep_%03d.rds", origin)))

# ── Extract and back-transform ─────────────────────────────────────────────────

fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

fc_state_orig  <- fc_state_draws[, 1:horizon, ] * Y_sd_w + Y_mean_w
fc_total_orig  <- fc_total_draws[, 1:horizon]   * Y_sd_w * n_state + Y_mean_w * n_state

actuals_state  <- Y_raw_state[test_idx, ]
actuals_total  <- Y_raw_total[test_idx]

# ── Compute per-horizon accuracy ───────────────────────────────────────────────

horizon_results <- lapply(1:horizon, function(h) {
  
  series_results <- lapply(1:n_state, function(s) {
    draws  <- fc_state_orig[, h, s]
    actual <- actuals_state[h, s]
    q05    <- quantile(draws, 0.05)
    q95    <- quantile(draws, 0.95)
    data.frame(
      origin  = origin,
      horizon = h,
      series  = state_codes[s],
      actual  = actual,
      fc_mean = mean(draws),
      ae      = abs(mean(draws) - actual),
      se      = (mean(draws) - actual)^2,
      crps    = scoringRules::crps_sample(actual, draws),
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
      origin  = origin,
      horizon = h,
      series  = "Total",
      actual  = actual_t,
      fc_mean = mean(draws_t),
      ae      = abs(mean(draws_t) - actual_t),
      se      = (mean(draws_t) - actual_t)^2,
      crps    = scoringRules::crps_sample(actual_t, draws_t),
      covered = as.integer(actual_t >= q05_t & actual_t <= q95_t)
    )
  )
})

# save results (only accuracy results at this point)
origin_df <- bind_rows(horizon_results)
saveRDS(origin_df, file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin)))