library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

load("tourism.Rdata")

training_end <- 96
m            <- 12L
horizon      <- 12

# ── Series identification ──────────────────────────────────────────────────────

state_codes  <- names(tourism)[grepl("^[A-G]$",      names(tourism))]
zone_codes   <- names(tourism)[grepl("^[A-G][A-Z]$", names(tourism))]
region_codes <- names(tourism)[grepl("^[A-G][A-Z][A-Z]$", names(tourism))]

n_state  <- length(state_codes)
n_zone   <- length(zone_codes)
n_region <- length(region_codes)

# Hierarchy mappings
zone_to_state    <- match(substr(zone_codes,   1, 1), state_codes)
region_to_zone   <- match(substr(region_codes, 1, 2), zone_codes)
region_to_state  <- match(substr(region_codes, 1, 1), state_codes)

# Counts for scaling
zones_per_state    <- as.integer(table(zone_to_state)[order(as.integer(names(table(zone_to_state))))])
regions_per_zone   <- as.integer(table(region_to_zone)[order(as.integer(names(table(region_to_zone))))])
regions_per_state  <- as.integer(table(region_to_state)[order(as.integer(names(table(region_to_state))))])
n_regions_total    <- n_region

# ── Raw data ───────────────────────────────────────────────────────────────────

Y_raw_state  <- tourism %>% select(all_of(state_codes))  %>% as.matrix()
Y_raw_zone   <- tourism %>% select(all_of(zone_codes))   %>% as.matrix()
Y_raw_region <- tourism %>% select(all_of(region_codes)) %>% as.matrix()
Y_raw_total  <- tourism %>% pull(Total)

# ── Scaling ────────────────────────────────────────────────────────────────────

Y_mean <- mean(Y_raw_state[1:training_end, ])
Y_sd   <- sd(Y_raw_state[1:training_end, ])

Y_state  <- (Y_raw_state - Y_mean) / Y_sd

Y_zone   <- sweep(Y_raw_zone,   2, Y_mean * zones_per_state[zone_to_state],     "-")
Y_zone   <- sweep(Y_zone,       2, Y_sd   * zones_per_state[zone_to_state],     "/")

Y_region <- sweep(Y_raw_region, 2, Y_mean * regions_per_zone[region_to_zone],   "-")
Y_region <- sweep(Y_region,     2, Y_sd   * regions_per_zone[region_to_zone],   "/")

Y_total_scaled <- (Y_raw_total - Y_mean * n_state) / (Y_sd * n_state)

# ── Stan model ─────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> n_zone;
  int<lower=1> n_region;
  int<lower=1> H;

  matrix[T_train, n_state]  Y_state;
  matrix[T_train, n_zone]   Y_zone;
  matrix[T_train, n_region] Y_region;
  vector[T_train]            Y_total;

  array[n_zone]   int<lower=1, upper=n_state> zone_to_state;
  array[n_region] int<lower=1, upper=n_zone>  region_to_zone;
}
parameters {
  // ── Top-level hyperpriors ─────────────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;

  // ── State-level hyperpriors ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_state] mu_alpha_state;
  vector<lower=0, upper=1>[n_state] mu_gamma_raw_state;

  // ── Zone-level hyperpriors (one per zone, drawn from parent state) ─────────
  vector<lower=0, upper=1>[n_zone] mu_alpha_zone;
  vector<lower=0, upper=1>[n_zone] mu_gamma_raw_zone;

  // ── State-level parameters ────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_raw_state;
  vector<lower=0>[n_state]          sigma_state;

  // ── Zone-level parameters ─────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zone]  alpha_zone;
  vector<lower=0, upper=1>[n_zone]  gamma_raw_zone;
  vector<lower=0>[n_zone]           sigma_zone;

  // ── Region-level parameters ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_region] alpha_region;
  vector<lower=0, upper=1>[n_region] gamma_raw_region;
  vector<lower=0>[n_region]          sigma_region;

  // ── Total-level parameters ────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;

  // ── Initial states ────────────────────────────────────────────────────────
  vector[n_state]     l0_state;
  matrix[n_state, m]  s0_state;
  vector[n_zone]      l0_zone;
  matrix[n_zone, m]   s0_zone;
  vector[n_region]    l0_region;
  matrix[n_region, m] s0_region;
  real                l0_total;
  vector[m]           s0_total;
}
transformed parameters {
  // ── Recover actual gamma values ───────────────────────────────────────────
  vector[n_state]  gamma_state;
  vector[n_zone]   gamma_zone;
  vector[n_region] gamma_region;
  real             gamma_total;

  for (n in 1:n_state)
    gamma_state[n]  = gamma_raw_state[n]  * (1 - alpha_state[n]);
  for (n in 1:n_zone)
    gamma_zone[n]   = gamma_raw_zone[n]   * (1 - alpha_zone[n]);
  for (n in 1:n_region)
    gamma_region[n] = gamma_raw_region[n] * (1 - alpha_region[n]);
  gamma_total       = gamma_raw_total     * (1 - alpha_total);

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

  // ── Region recursions ─────────────────────────────────────────────────────
  array[n_region] vector[T_train + 1] l_r;
  array[n_region] vector[T_train + m] s_r;

  for (n in 1:n_region) {
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

  // ── State-level hyperpriors drawn from top ────────────────────────────────
  mu_alpha_state     ~ normal(mu_alpha,     0.1);
  mu_gamma_raw_state ~ normal(mu_gamma_raw, 0.1);

  // ── Zone-level hyperpriors drawn from parent state ────────────────────────
  for (n in 1:n_zone) {
    int s = zone_to_state[n];
    mu_alpha_zone[n]     ~ normal(mu_alpha_state[s],     0.1);
    mu_gamma_raw_zone[n] ~ normal(mu_gamma_raw_state[s], 0.1);
  }

  // ── State series drawn from top hyperprior ────────────────────────────────
  alpha_state     ~ normal(mu_alpha,     0.1);
  gamma_raw_state ~ normal(mu_gamma_raw, 0.1);
  sigma_state     ~ cauchy(0, 0.3);

  // ── Zone series drawn from parent state hyperprior ────────────────────────
  for (n in 1:n_zone) {
    int s = zone_to_state[n];
    alpha_zone[n]     ~ normal(mu_alpha_state[s],     0.1);
    gamma_raw_zone[n] ~ normal(mu_gamma_raw_state[s], 0.1);
  }
  sigma_zone ~ cauchy(0, 0.3);

  // ── Region series drawn from parent zone hyperprior ───────────────────────
  for (n in 1:n_region) {
    int z = region_to_zone[n];
    alpha_region[n]     ~ normal(mu_alpha_zone[z],     0.1);
    gamma_raw_region[n] ~ normal(mu_gamma_raw_zone[z], 0.1);
  }
  sigma_region ~ cauchy(0, 0.3);

  // ── Total drawn from top hyperprior ───────────────────────────────────────
  alpha_total     ~ normal(mu_alpha,     0.1);
  gamma_raw_total ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 0.3);

  // ── Initial state priors ──────────────────────────────────────────────────
  l0_state  ~ normal(0, 1);
  l0_zone   ~ normal(0, 1);
  l0_region ~ normal(0, 1);
  l0_total  ~ normal(0, 1);

  for (n in 1:n_state) {
    s0_state[n]  ~ normal(0, 1);
    sum(s0_state[n])  ~ normal(0, 0.01 * m);
  }
  for (n in 1:n_zone) {
    s0_zone[n]   ~ normal(0, 1);
    sum(s0_zone[n])   ~ normal(0, 0.01 * m);
  }
  for (n in 1:n_region) {
    s0_region[n] ~ normal(0, 1);
    sum(s0_region[n]) ~ normal(0, 0.01 * m);
  }
  s0_total ~ normal(0, 1);
  sum(s0_total) ~ normal(0, 0.01 * m);

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_state)
    for (i in 1:T_train)
      Y_state[i, n] ~ normal(l_s[n][i] + s_s[n][i], sigma_state[n]);

  for (n in 1:n_zone)
    for (i in 1:T_train)
      Y_zone[i, n] ~ normal(l_z[n][i] + s_z[n][i], sigma_zone[n]);

  for (n in 1:n_region)
    for (i in 1:T_train)
      Y_region[i, n] ~ normal(l_r[n][i] + s_r[n][i], sigma_region[n]);

  for (i in 1:T_train)
    Y_total[i] ~ normal(l_t[i] + s_t[i], sigma_total);
}
generated quantities {
  matrix[H, n_state]  Y_fc_state;
  matrix[H, n_zone]   Y_fc_zone;
  matrix[H, n_region] Y_fc_region;
  vector[H]           Y_fc_total;

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

  for (n in 1:n_region) {
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
  T_train        = training_end,
  m              = m,
  n_state        = n_state,
  n_zone         = n_zone,
  n_region       = n_region,
  H              = horizon + 1,
  Y_state        = Y_state[1:training_end, ],
  Y_zone         = Y_zone[1:training_end, ],
  Y_region       = Y_region[1:training_end, ],
  Y_total        = Y_total_scaled[1:training_end],
  zone_to_state  = zone_to_state,
  region_to_zone = region_to_zone
)

stanDso <- stan_model(model_code = model_string)

save_seed <- 2025
start_time <- Sys.time()
stanFit <- sampling(
  object = stanDso,
  data   = data_list,
  seed   = save_seed,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  cores  = 4
)
end_time <- Sys.time()
print(end_time - start_time)

cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")

# ── Accuracy ───────────────────────────────────────────────────────────────────

fc_state_draws  <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_zone_draws   <- rstan::extract(stanFit, pars = "Y_fc_zone")$Y_fc_zone
fc_region_draws <- rstan::extract(stanFit, pars = "Y_fc_region")$Y_fc_region
fc_total_draws  <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

fc_state_means  <- apply(fc_state_draws,  c(2, 3), mean)[1:horizon, ]
fc_zone_means   <- apply(fc_zone_draws,   c(2, 3), mean)[1:horizon, ]
fc_region_means <- apply(fc_region_draws, c(2, 3), mean)[1:horizon, ]
fc_total_means  <- colMeans(fc_total_draws)[1:horizon]

fc_state_orig  <- fc_state_means * Y_sd + Y_mean
fc_zone_orig   <- sweep(fc_zone_means,   2, Y_sd * zones_per_state[zone_to_state],   "*")
fc_zone_orig   <- sweep(fc_zone_orig,    2, Y_mean * zones_per_state[zone_to_state], "+")
fc_region_orig <- sweep(fc_region_means, 2, Y_sd * regions_per_zone[region_to_zone],   "*")
fc_region_orig <- sweep(fc_region_orig,  2, Y_mean * regions_per_zone[region_to_zone], "+")
fc_total_orig  <- fc_total_means * Y_sd * n_state + Y_mean * n_state

actuals_state  <- Y_raw_state[ (training_end + 1):(training_end + horizon), ]
actuals_zone   <- Y_raw_zone[  (training_end + 1):(training_end + horizon), ]
actuals_region <- Y_raw_region[(training_end + 1):(training_end + horizon), ]
actuals_total  <- Y_raw_total[ (training_end + 1):(training_end + horizon)]

accuracy_state <- lapply(1:n_state, function(s) {
  errors <- fc_state_orig[, s] - actuals_state[, s]
  data.frame(series = state_codes[s], horizon = 1:horizon,
             actual = actuals_state[, s], forecast = fc_state_orig[, s],
             ae = abs(errors), se = errors^2)
}) %>% bind_rows()

accuracy_zone <- lapply(1:n_zone, function(z) {
  errors <- fc_zone_orig[, z] - actuals_zone[, z]
  data.frame(series = zone_codes[z], horizon = 1:horizon,
             actual = actuals_zone[, z], forecast = fc_zone_orig[, z],
             ae = abs(errors), se = errors^2)
}) %>% bind_rows()

accuracy_region <- lapply(1:n_region, function(r) {
  errors <- fc_region_orig[, r] - actuals_region[, r]
  data.frame(series = region_codes[r], horizon = 1:horizon,
             actual = actuals_region[, r], forecast = fc_region_orig[, r],
             ae = abs(errors), se = errors^2)
}) %>% bind_rows()

total_errors <- fc_total_orig - actuals_total
accuracy_total <- data.frame(
  series = "Total", horizon = 1:horizon,
  actual = actuals_total, forecast = fc_total_orig,
  ae = abs(total_errors), se = total_errors^2
)

all_accuracy <- bind_rows(accuracy_state, accuracy_zone, accuracy_region, accuracy_total)

summary_accuracy <- all_accuracy %>%
  mutate(level = case_when(
    series == "Total"   ~ "Total",
    nchar(series) == 1  ~ "State",
    nchar(series) == 2  ~ "Zone",
    nchar(series) == 3  ~ "Region"
  )) %>%
  group_by(level, series) %>%
  summarise(MAE = mean(ae), RMSE = sqrt(mean(se)), .groups = "drop") %>%
  bind_rows(
    all_accuracy %>%
      mutate(level = case_when(
        series == "Total"   ~ "Total",
        nchar(series) == 1  ~ "State",
        nchar(series) == 2  ~ "Zone",
        nchar(series) == 3  ~ "Region"
      )) %>%
      group_by(level) %>%
      summarise(series = "Overall (level)",
                MAE = mean(ae), RMSE = sqrt(mean(se)), .groups = "drop")
  ) %>%
  bind_rows(
    all_accuracy %>%
      summarise(level = "All", series = "Overall",
                MAE = mean(ae), RMSE = sqrt(mean(se)))
  )

print(summary_accuracy)

save(summary_accuracy, stanFit,
     file = "region_hierarchy_test.RData")