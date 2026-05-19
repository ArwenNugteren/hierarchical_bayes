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


# ── Scaling ────────────────────────────────────────────────────────────────────
# Scale defined from state-level training data as before
Y_mean <- mean(Y_raw_zone[train_idx, ])
Y_sd   <- sd(Y_raw_zone[train_idx, ])

# Zones: default scaling
Y_zone <- (Y_raw_zone - Y_mean) / Y_sd

# States: scaled by zones_per_state
Y_state <- sweep(Y_raw_state, 2, zones_per_state, function(x, n) {
  (x - (Y_mean * n)) / Y_sd
})

# Total: scaled by n_zones
Y_total_scaled <- (Y_raw_total - Y_mean * n_zone) / Y_sd

# hierarchical weighting
lambda_total <- 1
lambda_state <- 1
sigma_coh <- 0.01 * n_zones

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
  real<lower=0> lambda_state;
  real<lower=0> sigma_coh; 
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

// ── Initial values ────────────────────────────────────────────────────────
  vector[n_state]        l0_state;
  matrix[n_state, m-1]   s0_state_free;
  vector[n_zone]         l0_zone;
  matrix[n_zone, m-1]    s0_zone_free;
  real                   l0_total;
  vector[m-1]            s0_total_free;
}
transformed parameters {
  // ── Reconstruct full seasonals with sum-to-zero enforced ──────────────────
  matrix[n_state, m] s0_state;
  matrix[n_zone, m]  s0_zone;
  vector[m]          s0_total;

  for (n in 1:n_state) {
    s0_state[n, 1:(m-1)] = s0_state_free[n];
    s0_state[n, m]        = -sum(s0_state_free[n]);
  }
  for (n in 1:n_zone) {
    s0_zone[n, 1:(m-1)] = s0_zone_free[n];
    s0_zone[n, m]        = -sum(s0_zone_free[n]);
  }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);
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

  // ── Initial priors ──────────────────────────────────────────────────
  // Zones: unit scale (bottom level)
  l0_zone  ~ normal(0, 1);
  to_vector(s0_zone_free) ~ normal(0, 1);

  // States: scaled by number of zones per state
  for (n in 1:n_state) {
    int n_zones_in_state = 0;
    for (z in 1:n_zone)
      if (zone_to_state[z] == n) n_zones_in_state += 1;
    l0_state[n]       ~ normal(0, n_zones_in_state);
    s0_state_free[n]  ~ normal(0, n_zones_in_state);
  }

  // Total: scaled by total number of zones
  l0_total      ~ normal(0, n_zone);
  s0_total_free ~ normal(0, n_zone);
  
  // -- constrain initial values
  for (k in 1:m){
    real diff_total = (l0_total + s0_total[k]) - (sum(l0_state) + sum(s0_state[,k]));
    real diff_middle = (sum(l0_state) + sum(s0_state[,k]) - (sum(l0_zone) + sum(s0_zone[,k])));
    
    target += normal_lpdf(diff_total | 0, sigma_coh);
    target += normal_lpdf(diff_middle | 0, sigma_coh);
  }
  

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_zone)
    for (i in 1:T_train){
      target += normal_lpdf(Y_zone[i,n] | l_z[n][i] + s_z[n][i], sigma_zone[n]);
    }
  for (n in 1:n_state)
    for (i in 1:T_train){
      target += lambda_state * normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
    }
  for (i in 1:T_train){
    target += lambda_total * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
  }
  
  // ── Coherence penalty ─────────────────────────────────────────────────────
  for (i in 1:T_train) {
    real fitted_total  = l_t[i] + s_t[i];
    real fitted_states = 0;
    real fitted_zones  = 0;
  
    for (n in 1:n_state)
      fitted_states += l_s[n][i] + s_s[n][i];
    for (n in 1:n_zones)
      fitted_zones += l_z[n][i] + s_z[n][i];
  
    // Total vs sum of states
    target += normal_lpdf(fitted_total  - fitted_states | 0, sigma_coh);
    // Sum of states vs sum of zones
    target += normal_lpdf(fitted_states - fitted_zones  | 0, sigma_coh);
  }
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
    lambda_total = lambda_total, 
    lambda_state = lambda_state,
    sigma_coh    = sigma_coh
  ),
  seed    = 2025,
  chains  = 10,
  iter    = 2000,
  warmup  = 1000,
  cores   = 10,
  refresh = 0
)

saveRDS(stanFit, file.path(dir_stanfits, sprintf("stanfit_timestep_%03d.rds", origin)))






# ── Extract and back-transform ─────────────────────────────────────────────────
n_es_subsample <- 1000L   # draws used for ES pairwise term; increase for publication

# ── Scaling constants (zone-level, matching estimation script) ─────────────────
Y_mean <- mean(Y_raw_zone[train_idx, ])
Y_sd   <- sd(Y_raw_zone[train_idx, ])

# ── Helpers ────────────────────────────────────────────────────────────────────

# Energy Score (standardised, subsampled pairwise term)
# Sub-setting is done before calling this function; standardisation is internal.
compute_energy_score <- function(draws_mat, actuals, sub_idx = NULL) {
  if (!is.null(sub_idx))
    draws_mat <- draws_mat[sub_idx, , drop = FALSE]
  col_sds    <- apply(draws_mat, 2, sd)
  draws_sc   <- sweep(draws_mat, 2, col_sds, "/")
  actuals_sc <- as.numeric(actuals) / col_sds
  diffs1     <- sweep(draws_sc, 2, actuals_sc, "-")
  t1         <- mean(sqrt(rowSums(diffs1^2)))
  pw         <- dist(draws_sc, method = "euclidean")
  t2         <- mean(pw)
  t1 - 0.5 * t2
}

# KL divergence helper: returns all three directions
# p_draws, q_draws: numeric vectors (1-D distributions)
kl_triple <- function(p_draws, q_draws) {
  pm <- matrix(p_draws, ncol = 1)
  qm <- matrix(q_draws, ncol = 1)
  pq <- tryCatch(kldest::kld_est_nn(pm, qm), error = function(e) NA_real_)
  qp <- tryCatch(kldest::kld_est_nn(qm, pm), error = function(e) NA_real_)
  list(
    kl_pq  = pq,
    kl_qp  = qp,
    kl_sym = 0.5 * (pq + qp)
  )
}

# ── extract MCMC samples ─────────────────────────────────────────────
fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_zone_draws  <- rstan::extract(stanFit, pars = "Y_fc_zone")$Y_fc_zone
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

n_draws <- nrow(fc_zone_draws)

# ── Back-transform to original scale ──────────────────────────────────────────
# Zones: Y_zone_scaled = (Y_raw_zone - Y_mean) / Y_sd
fc_zone_orig <- fc_zone_draws[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean

# States: Y_state_scaled = (Y_raw_state - Y_mean * zones_per_state) / Y_sd
# So: Y_raw_state = Y_state_scaled * Y_sd + Y_mean * zones_per_state
fc_state_orig <- array(NA_real_, dim = c(n_draws, horizon, n_state))
for (s in 1:n_state) {
  fc_state_orig[, , s] <- fc_state_draws[, 1:horizon, s] * Y_sd +
    Y_mean * zones_per_state[s]
}

# Total: Y_total_scaled = (Y_raw_total - Y_mean * n_zone) / Y_sd
fc_total_orig <- fc_total_draws[, 1:horizon, drop = FALSE] * Y_sd +
  Y_mean * n_zone

# ── Actuals ────────────────────────────────────────────────────────────────────

actuals_state <- Y_raw_state[test_idx, ]   # horizon x n_state
actuals_zone  <- Y_raw_zone[test_idx, ]    # horizon x n_zone
actuals_total <- Y_raw_total[test_idx]     # horizon

# ── Fixed subsample index for ES ──────────────────────────────────────────────

set.seed(origin) # mix it up a little
sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)

# ── Per-horizon accuracy ───────────────────────────────────────────────────────

horizon_results <- lapply(1:horizon, function(h) {
  
  # ── Univariate: zones ──────────────────────────────────────────────────────
  zone_rows <- lapply(1:n_zone, function(z) {
    draws  <- fc_zone_orig[, h, z]
    actual <- actuals_zone[h, z]
    data.frame(
      origin  = origin,
      horizon = h,
      series  = zone_codes[z],
      level   = "Zone",
      actual  = actual,
      fc_mean = mean(draws),
      ae      = abs(mean(draws) - actual),
      se      = (mean(draws) - actual)^2,
      pit     = mean(draws < actual)
    )
  })
  
  # ── Univariate: states ─────────────────────────────────────────────────────
  state_rows <- lapply(1:n_state, function(s) {
    draws  <- fc_state_orig[, h, s]
    actual <- actuals_state[h, s]
    data.frame(
      origin  = origin,
      horizon = h,
      series  = state_codes[s],
      level   = "State",
      actual  = actual,
      fc_mean = mean(draws),
      ae      = abs(mean(draws) - actual),
      se      = (mean(draws) - actual)^2,
      pit     = mean(draws < actual)
    )
  })
  
  # ── Univariate: total ──────────────────────────────────────────────────────
  draws_t  <- fc_total_orig[, h]
  actual_t <- actuals_total[h]
  
  total_row <- data.frame(
    origin  = origin,
    horizon = h,
    series  = "Total",
    level   = "Total",
    actual  = actual_t,
    fc_mean = mean(draws_t),
    ae      = abs(mean(draws_t) - actual_t),
    se      = (mean(draws_t) - actual_t)^2,
    pit     = mean(draws_t < actual_t)
  )
  
  univariate_df <- bind_rows(c(zone_rows, state_rows, list(total_row)))
  
  # ── Draw matrices for this horizon ────────────────────────────────────────
  zone_mat  <- fc_zone_orig[, h, ]    # n_draws x n_zone
  state_mat <- fc_state_orig[, h, ]   # n_draws x n_state
  
  # ── Energy Score: all series jointly ──────────────────────────────────────
  # joint vector: zones | states | total  (n_zone + n_state + 1 columns)
  joint_mat     <- cbind(zone_mat, state_mat, draws_t)
  actuals_joint <- c(actuals_zone[h, ], actuals_state[h, ], actual_t)
  es_val        <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
  
  # ── Incoherence ───────────────────────────────────────────────────────────
  sum_zones_draws  <- rowSums(zone_mat)    # n_draws: sum of all zones
  sum_states_draws <- rowSums(state_mat)   # n_draws: sum of all states
  
  # Total vs sum of states
  incoherence_total_state <- mean(draws_t)       - mean(sum_states_draws)
  # Sum of states vs sum of zones
  incoherence_state_zone  <- mean(sum_states_draws) - mean(sum_zones_draws)
  # Total vs sum of zones (derived, but reported separately)
  incoherence_total_zone  <- mean(draws_t)       - mean(sum_zones_draws)
  
  # ── KL divergence: aggregate level ────────────────────────────────────────
  # KL(total || sum_states): p = total, q = sum_states
  kl_total_state <- kl_triple(draws_t,          sum_states_draws)
  # KL(sum_states || sum_zones): p = sum_states, q = sum_zones
  kl_state_zone  <- kl_triple(sum_states_draws, sum_zones_draws)
  # KL(total || sum_zones): p = total, q = sum_zones
  kl_total_zone  <- kl_triple(draws_t,          sum_zones_draws)
  
  # ── KL divergence: per-state (state vs its own zones) ─────────────────────
  kl_per_state <- lapply(1:n_state, function(s) {
    state_draws_s <- fc_state_orig[, h, s]
    # sum of zones belonging to state s
    zone_idx_s    <- which(zone_to_state == s)
    sum_zones_s   <- rowSums(fc_zone_orig[, h, zone_idx_s, drop = FALSE])
    kl            <- kl_triple(state_draws_s, sum_zones_s)
    data.frame(
      state       = state_codes[s],
      kl_pq       = kl$kl_pq,    # KL(state || sum_zones_s)
      kl_qp       = kl$kl_qp,    # KL(sum_zones_s || state)
      kl_sym      = kl$kl_sym
    )
  })
  kl_per_state_df <- bind_rows(kl_per_state)
  
  # ── Assemble multivariate row ──────────────────────────────────────────────
  multivariate_df <- data.frame(
    origin  = origin,
    horizon = h,
    es      = es_val,
    # Incoherence
    incoherence_total_state = incoherence_total_state,
    incoherence_state_zone  = incoherence_state_zone,
    incoherence_total_zone  = incoherence_total_zone,
    # Aggregate KL: total vs sum of states
    kl_ts_pq  = kl_total_state$kl_pq,
    kl_ts_qp  = kl_total_state$kl_qp,
    kl_ts_sym = kl_total_state$kl_sym,
    # Aggregate KL: sum of states vs sum of zones
    kl_sz_pq  = kl_state_zone$kl_pq,
    kl_sz_qp  = kl_state_zone$kl_qp,
    kl_sz_sym = kl_state_zone$kl_sym,
    # Aggregate KL: total vs sum of zones
    kl_tz_pq  = kl_total_zone$kl_pq,
    kl_tz_qp  = kl_total_zone$kl_qp,
    kl_tz_sym = kl_total_zone$kl_sym
  )
  
  list(
    univariate      = univariate_df,
    multivariate    = multivariate_df,
    kl_per_state    = cbind(data.frame(origin = origin, horizon = h),
                            kl_per_state_df)
  )
})

# ── Collect and save ───────────────────────────────────────────────────────────

origin_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
origin_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))
origin_kl_state     <- bind_rows(lapply(horizon_results, `[[`, "kl_per_state"))

saveRDS(
  list(
    univariate   = origin_univariate,
    multivariate = origin_multivariate,
    kl_per_state = origin_kl_state
  ),
  file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin))
)
