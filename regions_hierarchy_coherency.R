library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

load("tourism.Rdata")

training_end <- 96
m            <- 12L
horizon      <- 12

# ── Series identification ──────────────────────────────────────────────────────

state_codes  <- names(tourism)[grepl("^[A-G]$",          names(tourism))]
zone_codes   <- names(tourism)[grepl("^[A-G][A-Z]$",     names(tourism))]
region_codes <- names(tourism)[grepl("^[A-G][A-Z][A-Z]$", names(tourism))]

n_states  <- length(state_codes)
n_zones   <- length(zone_codes)
n_regions <- length(region_codes)

# ── Hierarchy mappings ─────────────────────────────────────────────────────────

zone_to_state   <- match(substr(zone_codes,   1, 1), state_codes)
region_to_zone  <- match(substr(region_codes, 1, 2), zone_codes)
region_to_state <- match(substr(region_codes, 1, 1), state_codes)

zones_per_state    <- as.integer(table(zone_to_state)[order(as.integer(names(table(zone_to_state))))])
regions_per_zone   <- as.integer(table(region_to_zone)[order(as.integer(names(table(region_to_zone))))])
regions_per_state  <- as.integer(table(region_to_state)[order(as.integer(names(table(region_to_state))))])

# ── Raw data ───────────────────────────────────────────────────────────────────

Y_raw_state  <- tourism %>% select(all_of(state_codes))  %>% as.matrix()
Y_raw_zone   <- tourism %>% select(all_of(zone_codes))   %>% as.matrix()
Y_raw_region <- tourism %>% select(all_of(region_codes)) %>% as.matrix()
Y_raw_total  <- tourism %>% pull(Total)

# ── Scaling (region-based) ─────────────────────────────────────────────────────

Y_mean <- mean(Y_raw_region[1:training_end, ])
Y_sd   <- sd(Y_raw_region[1:training_end, ])

Y_region <- (Y_raw_region - Y_mean) / Y_sd

Y_zone <- sweep(Y_raw_zone, 2, Y_mean * regions_per_zone, "-")
Y_zone <- sweep(Y_zone,     2, Y_sd,                      "/")

Y_state <- sweep(Y_raw_state, 2, Y_mean * regions_per_state, "-")
Y_state <- sweep(Y_state,     2, Y_sd,                       "/")

Y_total <- (Y_raw_total - Y_mean * n_regions) / Y_sd
lambda_total <- 1
lambda_state <- 1
lambda_zone <- 1
sigma_coh <- 0.01 * n_regions

# ── Stan model ─────────────────────────────────────────────────────────────────
model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_states;
  int<lower=1> n_zones;
  int<lower=1> n_regions;
  int<lower=1> H;

  matrix[T_train, n_states]  Y_state;
  matrix[T_train, n_zones]   Y_zone;
  matrix[T_train, n_regions] Y_region;
  vector[T_train]             Y_total;

  array[n_zones]   int<lower=1, upper=n_states> zone_to_state;
  array[n_regions] int<lower=1, upper=n_zones>  region_to_zone;

  array[n_zones]  int<lower=1> regions_per_zone;
  array[n_states] int<lower=1> regions_per_state;

  real<lower=0> lambda_total;
  real<lower=0> lambda_state;
  real<lower=0> lambda_zone;
  real<lower=0> sigma_coh;
}
parameters {
  // ── Top-level hyperpriors ─────────────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;

  // ── State-level hyperpriors ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_states] mu_alpha_state;
  vector<lower=0, upper=1>[n_states] mu_gamma_raw_state;

  // ── Zone-level hyperpriors ────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zones] mu_alpha_zone;
  vector<lower=0, upper=1>[n_zones] mu_gamma_raw_zone;

  // ── State-level parameters ────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_states] alpha_state;
  vector<lower=0, upper=1>[n_states] gamma_raw_state;
  vector<lower=0>[n_states]          sigma_state;

  // ── Zone-level parameters ─────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zones] alpha_zone;
  vector<lower=0, upper=1>[n_zones] gamma_raw_zone;
  vector<lower=0>[n_zones]          sigma_zone;

  // ── Region-level parameters ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_regions] alpha_region;
  vector<lower=0, upper=1>[n_regions] gamma_raw_region;
  vector<lower=0>[n_regions]          sigma_region;

  // ── Total-level parameters ────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;

  // ── Initial states (free seasonal components, m-1) ────────────────────────
  vector[n_states]          l0_state;
  matrix[n_states,  m - 1]  s0_state_free;
  vector[n_zones]           l0_zone;
  matrix[n_zones,   m - 1]  s0_zone_free;
  vector[n_regions]         l0_region;
  matrix[n_regions, m - 1]  s0_region_free;
  real                      l0_total;
  vector[m - 1]             s0_total_free;
}
transformed parameters {
  // ── Reconstruct full seasonals with sum-to-zero enforced ──────────────────
  matrix[n_states,  m] s0_state;
  matrix[n_zones,   m] s0_zone;
  matrix[n_regions, m] s0_region;
  vector[m]            s0_total;

  for (n in 1:n_states) {
    s0_state[n, 1:(m-1)] = s0_state_free[n];
    s0_state[n, m]        = -sum(s0_state_free[n]);
  }
  for (n in 1:n_zones) {
    s0_zone[n, 1:(m-1)] = s0_zone_free[n];
    s0_zone[n, m]        = -sum(s0_zone_free[n]);
  }
  for (n in 1:n_regions) {
    s0_region[n, 1:(m-1)] = s0_region_free[n];
    s0_region[n, m]        = -sum(s0_region_free[n]);
  }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);

  // ── Recover actual gamma values ───────────────────────────────────────────
  vector[n_states]  gamma_state;
  vector[n_zones]   gamma_zone;
  vector[n_regions] gamma_region;
  real              gamma_total;

  for (n in 1:n_states)
    gamma_state[n]  = gamma_raw_state[n]  * (1 - alpha_state[n]);
  for (n in 1:n_zones)
    gamma_zone[n]   = gamma_raw_zone[n]   * (1 - alpha_zone[n]);
  for (n in 1:n_regions)
    gamma_region[n] = gamma_raw_region[n] * (1 - alpha_region[n]);
  gamma_total       = gamma_raw_total     * (1 - alpha_total);

  // ── State recursions ──────────────────────────────────────────────────────
  array[n_states] vector[T_train + 1] l_s;
  array[n_states] vector[T_train + m] s_s;

  for (n in 1:n_states) {
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
  array[n_zones] vector[T_train + 1] l_z;
  array[n_zones] vector[T_train + m] s_z;

  for (n in 1:n_zones) {
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

  // ── Region recursions ─────────────────────────────────────────────────────
  array[n_regions] vector[T_train + 1] l_r;
  array[n_regions] vector[T_train + m] s_r;

  for (n in 1:n_regions) {
    for (j in 1:m) s_r[n][j] = s0_region[n, j];
    l_r[n][1] = l0_region[n];
    for (i in 1:T_train) {
      real s_lag   = s_r[n][i];
      l_r[n][i+1] = alpha_region[n] * (Y_region[i, n] - s_lag)
                  + (1 - alpha_region[n]) * l_r[n][i];
      s_r[n][i+m] = gamma_region[n] * (Y_region[i, n] - l_r[n][i])
                  + (1 - gamma_region[n]) * s_lag;
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

  // ── State-level hyperpriors ───────────────────────────────────────────────
  mu_alpha_state     ~ normal(mu_alpha,     0.1);
  mu_gamma_raw_state ~ normal(mu_gamma_raw, 0.1);

  // ── Zone-level hyperpriors ────────────────────────────────────────────────
  for (n in 1:n_zones) {
    int s = zone_to_state[n];
    mu_alpha_zone[n]     ~ normal(mu_alpha_state[s],     0.1);
    mu_gamma_raw_zone[n] ~ normal(mu_gamma_raw_state[s], 0.1);
  }

  // ── Total parameters ──────────────────────────────────────────────────────
  alpha_total     ~ normal(mu_alpha,     0.1);
  gamma_raw_total ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 1);

  // ── State parameters ──────────────────────────────────────────────────────
  alpha_state     ~ normal(mu_alpha,     0.1);
  gamma_raw_state ~ normal(mu_gamma_raw, 0.1);
  sigma_state     ~ cauchy(0, 0.5);

  // ── Zone parameters ───────────────────────────────────────────────────────
  for (n in 1:n_zones) {
    int s = zone_to_state[n];
    alpha_zone[n]     ~ normal(mu_alpha_state[s],     0.1);
    gamma_raw_zone[n] ~ normal(mu_gamma_raw_state[s], 0.1);
  }
  sigma_zone ~ cauchy(0, 0.3);

  // ── Region parameters ─────────────────────────────────────────────────────
  for (n in 1:n_regions) {
    int z = region_to_zone[n];
    alpha_region[n]     ~ normal(mu_alpha_zone[z],     0.1);
    gamma_raw_region[n] ~ normal(mu_gamma_raw_zone[z], 0.1);
  }
  sigma_region ~ cauchy(0, 0.3);

  // ── Initial value priors ──────────────────────────────────────────────────
  l0_region ~ normal(0, 1);
  to_vector(s0_region_free) ~ normal(0, 1);

  for (n in 1:n_zones) {
    l0_zone[n]      ~ normal(0, regions_per_zone[n]);
    s0_zone_free[n] ~ normal(0, regions_per_zone[n]);
  }
  for (n in 1:n_states) {
    l0_state[n]      ~ normal(0, regions_per_state[n]);
    s0_state_free[n] ~ normal(0, regions_per_state[n]);
  }
  l0_total      ~ normal(0, n_regions);
  s0_total_free ~ normal(0, n_regions);

  // ── Initial-value coherence penalties ─────────────────────────────────────
  for (k in 1:m) {
    real init_total   = l0_total        + s0_total[k];
    real init_states  = sum(l0_state)   + sum(s0_state[, k]);
    real init_zones   = sum(l0_zone)    + sum(s0_zone[, k]);
    real init_regions = sum(l0_region)  + sum(s0_region[, k]);

    // Total vs sum of states
    target += normal_lpdf(init_total   - init_states  | 0, sigma_coh);
    // Sum of states vs sum of zones
    target += normal_lpdf(init_states  - init_zones   | 0, sigma_coh);
    // Sum of zones vs sum of regions
    target += normal_lpdf(init_zones   - init_regions | 0, sigma_coh);
  }

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_regions)
    for (i in 1:T_train)
      target += normal_lpdf(Y_region[i, n] | l_r[n][i] + s_r[n][i], sigma_region[n]);

  for (n in 1:n_zones)
    for (i in 1:T_train)
      target += lambda_zone * normal_lpdf(Y_zone[i, n] | l_z[n][i] + s_z[n][i], sigma_zone[n]);

  for (n in 1:n_states)
    for (i in 1:T_train)
      target += lambda_state * normal_lpdf(Y_state[i, n] | l_s[n][i] + s_s[n][i], sigma_state[n]);

  for (i in 1:T_train)
    target += lambda_total * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);

  // ── Fitted-value coherence penalties ──────────────────────────────────────
  for (i in 1:T_train) {
    real fitted_total   = l_t[i] + s_t[i];
    real fitted_states  = 0;
    real fitted_zones   = 0;
    real fitted_regions = 0;

    for (n in 1:n_states)
      fitted_states  += l_s[n][i] + s_s[n][i];
    for (n in 1:n_zones)
      fitted_zones   += l_z[n][i] + s_z[n][i];
    for (n in 1:n_regions)
      fitted_regions += l_r[n][i] + s_r[n][i];

    // Total vs sum of states
    target += normal_lpdf(fitted_total  - fitted_states  | 0, sigma_coh);
    // Sum of states vs sum of zones
    target += normal_lpdf(fitted_states - fitted_zones   | 0, sigma_coh);
    // Sum of zones vs sum of regions
    target += normal_lpdf(fitted_zones  - fitted_regions | 0, sigma_coh);
  }
}
generated quantities {
  matrix[H, n_states]  Y_fc_state;
  matrix[H, n_zones]   Y_fc_zone;
  matrix[H, n_regions] Y_fc_region;
  vector[H]            Y_fc_total;

  for (n in 1:n_states) {
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

  for (n in 1:n_zones) {
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

  for (n in 1:n_regions) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_r[n][T_train + 1];
    for (h in 1:H) sf[h] = s_r[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_region[h, n] = normal_rng(lf[h] + sf[h], sigma_region[n]);
      lf[h+1] = alpha_region[n] * (Y_fc_region[h, n] - sf[h])
              + (1 - alpha_region[n]) * lf[h];
    }
    Y_fc_region[H, n] = normal_rng(lf[H] + sf[H], sigma_region[n]);
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

# ── Data list ──────────────────────────────────────────────────────────────────

data_list <- list(
  T_train          = training_end,
  m                = m,
  n_states         = n_states,
  n_zones          = n_zones,
  n_regions        = n_regions,
  H                = horizon + 1,
  Y_state          = Y_state[1:training_end, ],
  Y_zone           = Y_zone[1:training_end, ],
  Y_region         = Y_region[1:training_end, ],
  Y_total          = Y_total[1:training_end],
  zone_to_state    = zone_to_state,
  region_to_zone   = region_to_zone,
  regions_per_zone  = regions_per_zone,
  regions_per_state = regions_per_state, 
  lambda_total      = lambda_total,
  lambda_state      = lambda_state, 
  lambda_zone       = lambda_zone, 
  sigma_coh         = sigma_coh
)

stanDso <- stan_model(model_code = model_string)

save_seed  <- 2025
start_time <- Sys.time()
stanFit <- sampling(
  object = stanDso,
  data   = data_list,
  seed   = save_seed,
  chains = 10,
  iter   = 200,
  warmup = 100,
  cores  = 10
)
end_time <- Sys.time()
print(end_time - start_time)

cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")





# ── Accuracy ───────────────────────────────────────────────────────────────────

Y_mean <- mean(Y_raw_region[1:training_end, ])
Y_sd   <- sd(Y_raw_region[1:training_end, ])

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

JSD_pair <- function(p_draws, q_draws, n_grid = 1000) {
  lo    <- min(c(p_draws, q_draws)) - 3 * sd(c(p_draws, q_draws))
  hi    <- max(c(p_draws, q_draws)) + 3 * sd(c(p_draws, q_draws))
  p_kde <- pmax(density(p_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  q_kde <- pmax(density(q_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  m_kde <- 0.5 * (p_kde + q_kde)
  dx    <- (hi - lo) / (n_grid - 1)
  jsd   <- 0.5 * (sum(p_kde * log(p_kde / m_kde)) +
                    sum(q_kde * log(q_kde / m_kde))) * dx
  max(0, min(log(2), jsd))
}

fc_region_draws <- rstan::extract(stanFit, pars = "Y_fc_region")$Y_fc_region
fc_zone_draws   <- rstan::extract(stanFit, pars = "Y_fc_zone")$Y_fc_zone
fc_state_draws  <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_total_draws  <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

n_draws <- nrow(fc_region_draws)

fc_region_orig <- fc_region_draws[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean

fc_zone_orig <- array(NA_real_, dim = c(n_draws, horizon, n_zones))
for (z in 1:n_zones)
  fc_zone_orig[, , z] <- fc_zone_draws[, 1:horizon, z] * Y_sd +
  Y_mean * regions_per_zone[z]

fc_state_orig <- array(NA_real_, dim = c(n_draws, horizon, n_states))
for (s in 1:n_states)
  fc_state_orig[, , s] <- fc_state_draws[, 1:horizon, s] * Y_sd +
  Y_mean * regions_per_state[s]

fc_total_orig <- fc_total_draws[, 1:horizon, drop = FALSE] * Y_sd +
  Y_mean * n_regions

actuals_region <- Y_raw_region[(training_end + 1):(training_end + horizon), ]
actuals_zone   <- Y_raw_zone[  (training_end + 1):(training_end + horizon), ]
actuals_state  <- Y_raw_state[ (training_end + 1):(training_end + horizon), ]
actuals_total  <- Y_raw_total[ (training_end + 1):(training_end + horizon)]

n_es_subsample <- 1000L
set.seed(2025)
sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)

horizon_results <- lapply(1:horizon, function(h) {
  
  make_rows <- function(draws_arr, actuals_mat, codes, level) {
    lapply(seq_len(ncol(actuals_mat)), function(i) {
      draws  <- draws_arr[, h, i]
      actual <- actuals_mat[h, i]
      data.frame(
        horizon    = h,
        series     = codes[i],
        level      = level,
        actual     = actual,
        fc_mean    = mean(draws),
        ae         = abs(mean(draws) - actual),
        se         = (mean(draws) - actual)^2,
        crps       = scoringRules::crps_sample(actual, draws),
        logS       = scoringRules::logs_sample(actual, draws),
        covered_90 = as.integer(actual >= quantile(draws, 0.05)  &
                                  actual <= quantile(draws, 0.95)),
        covered_95 = as.integer(actual >= quantile(draws, 0.025) &
                                  actual <= quantile(draws, 0.975)),
        pit        = mean(draws < actual)
      )
    })
  }
  
  draws_t  <- fc_total_orig[, h]
  actual_t <- actuals_total[h]
  total_row <- data.frame(
    horizon    = h,
    series     = "Total",
    level      = "Total",
    actual     = actual_t,
    fc_mean    = mean(draws_t),
    ae         = abs(mean(draws_t) - actual_t),
    se         = (mean(draws_t) - actual_t)^2,
    crps       = scoringRules::crps_sample(actual_t, draws_t),
    logS       = scoringRules::logs_sample(actual_t, draws_t),
    covered_90 = as.integer(actual_t >= quantile(draws_t, 0.05)  &
                              actual_t <= quantile(draws_t, 0.95)),
    covered_95 = as.integer(actual_t >= quantile(draws_t, 0.025) &
                              actual_t <= quantile(draws_t, 0.975)),
    pit        = mean(draws_t < actual_t)
  )
  
  univariate_df <- bind_rows(c(
    make_rows(fc_region_orig, actuals_region, region_codes, "Region"),
    make_rows(fc_zone_orig,   actuals_zone,   zone_codes,   "Zone"),
    make_rows(fc_state_orig,  actuals_state,  state_codes,  "State"),
    list(total_row)
  ))
  
  region_mat <- fc_region_orig[, h, ]
  zone_mat   <- fc_zone_orig[,   h, ]
  state_mat  <- fc_state_orig[,  h, ]
  
  sum_regions_draws <- rowSums(region_mat)
  sum_zones_draws   <- rowSums(zone_mat)
  sum_states_draws  <- rowSums(state_mat)
  
  joint_mat     <- cbind(region_mat, zone_mat, state_mat, draws_t)
  actuals_joint <- c(actuals_region[h, ], actuals_zone[h, ],
                     actuals_state[h, ],  actual_t)
  es_val <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
  
  mai_region_zone  <- mean(sum_regions_draws) - mean(sum_zones_draws)
  mai_zone_state   <- mean(sum_zones_draws)   - mean(sum_states_draws)
  mai_state_total  <- mean(sum_states_draws)  - mean(draws_t)
  mai_region_total <- mean(sum_regions_draws) - mean(draws_t)
  
  jsd_region_zone  <- JSD_pair(sum_regions_draws, sum_zones_draws)
  jsd_zone_state   <- JSD_pair(sum_zones_draws,   sum_states_draws)
  jsd_state_total  <- JSD_pair(sum_states_draws,  draws_t)
  jsd_region_total <- JSD_pair(sum_regions_draws, draws_t)
  
  multivariate_df <- data.frame(
    horizon          = h,
    es               = es_val,
    mai_region_zone  = mai_region_zone,
    mai_zone_state   = mai_zone_state,
    mai_state_total  = mai_state_total,
    mai_region_total = mai_region_total,
    jsd_region_zone  = jsd_region_zone,
    jsd_zone_state   = jsd_zone_state,
    jsd_state_total  = jsd_state_total,
    jsd_region_total = jsd_region_total
  )
  
  list(univariate = univariate_df, multivariate = multivariate_df)
})

origin_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
origin_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))

make_level_row <- function(lev) {
  origin_univariate %>%
    filter(level == lev) %>%
    summarise(
      level      = lev,
      MAE        = mean(ae),
      RMSE       = sqrt(mean(se)),
      CRPS       = mean(crps),
      LogS       = mean(logS),
      Cover_90   = mean(covered_90),
      Cover_95   = mean(covered_95),
      ES               = NA_real_,
      MAI_region_zone  = NA_real_,
      MAI_zone_state   = NA_real_,
      MAI_state_total  = NA_real_,
      MAI_region_total = NA_real_,
      JSD_region_zone  = NA_real_,
      JSD_zone_state   = NA_real_,
      JSD_state_total  = NA_real_,
      JSD_region_total = NA_real_
    )
}

overall_row <- bind_cols(
  origin_univariate %>%
    summarise(
      level    = "Overall",
      MAE      = mean(ae),
      RMSE     = sqrt(mean(se)),
      CRPS     = mean(crps),
      LogS     = mean(logS),
      Cover_90 = mean(covered_90),
      Cover_95 = mean(covered_95)
    ),
  origin_multivariate %>%
    summarise(
      ES               = mean(es),
      MAI_region_zone  = mean(abs(mai_region_zone)),
      MAI_zone_state   = mean(abs(mai_zone_state)),
      MAI_state_total  = mean(abs(mai_state_total)),
      MAI_region_total = mean(abs(mai_region_total)),
      JSD_region_zone  = mean(jsd_region_zone,  na.rm = TRUE),
      JSD_zone_state   = mean(jsd_zone_state,   na.rm = TRUE),
      JSD_state_total  = mean(jsd_state_total,  na.rm = TRUE),
      JSD_region_total = mean(jsd_region_total, na.rm = TRUE)
    )
)

summary_table <- bind_rows(
  make_level_row("Region"),
  make_level_row("Zone"),
  make_level_row("State"),
  make_level_row("Total"),
  overall_row
)

print(summary_table)

save(stanFit, summary_table, origin_univariate, origin_multivariate,
     file = "region_hierarchy_test2.RData")