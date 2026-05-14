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

state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)

zone_codes  <- names(tourism)[grepl("^[A-G][A-Z]$", names(tourism))]
n_zone      <- length(zone_codes)

zone_to_state <- match(substr(zone_codes, 1, 1), state_codes)

zones_per_state <- as.integer(table(zone_to_state)[order(as.integer(names(table(zone_to_state))))])

Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_zone  <- tourism %>% select(all_of(zone_codes))  %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

m           <- 12L
window_size <- 96
horizon     <- 12
test_set    <- 121
origins     <- 96:(96 + test_set - 1)

# ── Get this task's origin from PBS array index ────────────────────────────────

task_id <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
origin  <- origins[task_id]

# ── Training window ────────────────────────────────────────────────────────────

train_idx <- (origin - window_size + 1):origin
test_idx  <- (origin + 1):(origin + horizon)

Y_mean_w <- mean(Y_raw_state[train_idx, ])
Y_sd_w   <- sd(Y_raw_state[train_idx, ])

Y_state_w <- (Y_raw_state[train_idx, ] - Y_mean_w) / Y_sd_w
Y_zone_w  <- sweep(Y_raw_zone[train_idx, ], 2, Y_mean_w * zones_per_state[zone_to_state], "-")
Y_zone_w  <- sweep(Y_zone_w,               2, Y_sd_w   * zones_per_state[zone_to_state], "/")
Y_total_w <- (Y_raw_total[train_idx] - Y_mean_w * n_state) / (Y_sd_w * n_state)


# hierarchical weighting
lambda_total <- 2


# ── Stan model ─────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> n_zone;
  int<lower=1> H;
  matrix[T_train, n_state] Y_state;
  matrix[T_train, n_zone]  Y_zone;
  vector[T_train]           Y_total;
  array[n_zone] int<lower=1, upper=n_state> zone_to_state;
  
  real<lower=0> lambda_total;
}
parameters {
  // ── Top-level hyperpriors ─────────────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;   // hyperprior on gamma fraction

  // ── State-level hyperpriors ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_state] mu_alpha_state;
  vector<lower=0, upper=1>[n_state] mu_gamma_raw_state;

  // ── State-level parameters ────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_raw_state;  // fraction of (1 - alpha)
  vector<lower=0>[n_state]          sigma_state;

  // ── Zone-level parameters ─────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zone]  alpha_zone;
  vector<lower=0, upper=1>[n_zone]  gamma_raw_zone;   // fraction of (1 - alpha)
  vector<lower=0>[n_zone]           sigma_zone;

  // ── Total-level parameters ────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;

  // ── Initial states ────────────────────────────────────────────────────────
  vector[n_state]    l0_state;
  matrix[n_state, m] s0_state;
  vector[n_zone]     l0_zone;
  matrix[n_zone, m]  s0_zone;
  real               l0_total;
  vector[m]          s0_total;
}
transformed parameters {
  // ── Recover actual gamma values ───────────────────────────────────────────
  vector[n_state] gamma_state;
  vector[n_zone]  gamma_zone;
  real            gamma_total;

  for (n in 1:n_state)
    gamma_state[n] = gamma_raw_state[n] * (1 - alpha_state[n]);
  for (n in 1:n_zone)
    gamma_zone[n]  = gamma_raw_zone[n]  * (1 - alpha_zone[n]);
  gamma_total      = gamma_raw_total    * (1 - alpha_total);

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

  // ── Zone recursions ───────────────────────────────────────────────────────
  array[n_zone] vector[T_train + 1] l_z;
  array[n_zone] vector[T_train + m] s_z;

  for (n in 1:n_zone) {
    for (j in 1:m) s_z[n][j] = s0_zone[n, j];
    l_z[n][1] = l0_zone[n];
    for (i in 1:T_train) {
      real s_lag   = s_z[n][i];
      l_z[n][i+1] = alpha_zone[n] * (Y_zone[i, n] - s_lag)
                  + (1 - alpha_zone[n]) * l_z[n][i];
      s_z[n][i+m] = gamma_zone[n] * (Y_zone[i, n] - l_z[n][i])
                  + (1 - gamma_zone[n]) * s_lag;
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
  // ── Top-level hyperpriors ─────────────────────────────────────────────────
  mu_alpha     ~ uniform(0, 1);
  mu_gamma_raw ~ uniform(0, 1);

  // ── State-level hyperpriors drawn from top ────────────────────────────────
  mu_alpha_state     ~ normal(mu_alpha,     0.1);
  mu_gamma_raw_state ~ normal(mu_gamma_raw, 0.1);

  // ── State-level parameters drawn from top hyperprior ─────────────────────
  alpha_state     ~ normal(mu_alpha,     0.1);
  gamma_raw_state ~ normal(mu_gamma_raw, 0.1);
  sigma_state     ~ cauchy(0, 0.3);

  // ── Zone-level parameters drawn from parent state hyperprior ─────────────
  for (n in 1:n_zone) {
    int s = zone_to_state[n];
    alpha_zone[n]     ~ normal(mu_alpha_state[s],     0.1);
    gamma_raw_zone[n] ~ normal(mu_gamma_raw_state[s], 0.1);
  }
  sigma_zone ~ cauchy(0, 0.3);

  // ── Total-level parameters drawn from top hyperprior ─────────────────────
  alpha_total     ~ normal(mu_alpha,     0.1);
  gamma_raw_total ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 0.3);

  // ── Initial state priors ──────────────────────────────────────────────────
  l0_state ~ normal(0, 1);
  l0_zone  ~ normal(0, 1);
  l0_total ~ normal(0, 1);

  for (n in 1:n_state) {
    s0_state[n] ~ normal(0, 1);
    sum(s0_state[n]) ~ normal(0, 0.01 * m);
  }
  for (n in 1:n_zone) {
    s0_zone[n] ~ normal(0, 1);
    sum(s0_zone[n]) ~ normal(0, 0.01 * m);
  }
  s0_total ~ normal(0, 1);
  sum(s0_total) ~ normal(0, 0.01 * m);

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_state)
    for (i in 1:T_train)
      //Y_state[i, n] ~ normal(l_s[n][i] + s_s[n][i], sigma_state[n]);
      target += normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);

  for (n in 1:n_zone)
    for (i in 1:T_train)
      //Y_zone[i, n] ~ normal(l_z[n][i] + s_z[n][i], sigma_zone[n]);
      target += normal_lpdf(Y_zone[i,n] | l_z[n][i] + s_z[n][i], sigma_zone[n]);

  for (i in 1:T_train)
    //Y_total[i] ~ normal(l_t[i] + s_t[i], sigma_total);
    target += lambda_total * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
}
generated quantities {
  matrix[H, n_state] Y_fc_state;
  matrix[H, n_zone]  Y_fc_zone;
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

  for (n in 1:n_zone) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_z[n][T_train + 1];
    for (h in 1:H) sf[h] = s_z[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_zone[h, n] = normal_rng(lf[h] + sf[h], sigma_zone[n]);
      lf[h+1] = alpha_zone[n] * (Y_fc_zone[h, n] - sf[h])
              + (1 - alpha_zone[n]) * lf[h];
    }
    Y_fc_zone[H, n] = normal_rng(lf[H] + sf[H], sigma_zone[n]);
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
  object = stanDso,
  data   = list(
    T_train       = window_size,
    m             = m,
    n_state       = n_state,
    n_zone        = n_zone,
    H             = horizon + 1,
    Y_state       = Y_state_w,
    Y_zone        = Y_zone_w,
    Y_total       = Y_total_w,
    zone_to_state = zone_to_state, 
    lambda_total = lambda_total
  ),
  seed    = 2025,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  refresh = 0
)

saveRDS(stanFit, file.path(dir_stanfits, sprintf("stanfit_timestep_%03d.rds", origin)))

# ── Extract and back-transform ─────────────────────────────────────────────────

fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_zone_draws  <- rstan::extract(stanFit, pars = "Y_fc_zone")$Y_fc_zone
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

fc_state_orig <- fc_state_draws[, 1:horizon, ] * Y_sd_w + Y_mean_w
fc_zone_orig  <- sweep(fc_zone_draws[, 1:horizon, ], 3,
                       Y_sd_w  * zones_per_state[zone_to_state], "*")
fc_zone_orig  <- sweep(fc_zone_orig, 3,
                       Y_mean_w * zones_per_state[zone_to_state], "+")
fc_total_orig <- fc_total_draws[, 1:horizon] * Y_sd_w * n_state + Y_mean_w * n_state

actuals_state <- Y_raw_state[test_idx, ]
actuals_zone  <- Y_raw_zone[test_idx, ]
actuals_total <- Y_raw_total[test_idx]

# ── Compute per-horizon accuracy ───────────────────────────────────────────────

horizon_results <- lapply(1:horizon, function(h) {
  
  state_results <- lapply(1:n_state, function(s) {
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
  
  zone_results <- lapply(1:n_zone, function(z) {
    draws  <- fc_zone_orig[, h, z]
    actual <- actuals_zone[h, z]
    q05    <- quantile(draws, 0.05)
    q95    <- quantile(draws, 0.95)
    data.frame(
      origin  = origin,
      horizon = h,
      series  = zone_codes[z],
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
    state_results,
    zone_results,
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

origin_df <- bind_rows(horizon_results)
saveRDS(origin_df, file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin)))