library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

load("tourism.Rdata")

training_end <- 96
m            <- 12L
horizon      <- 12

# ── Series identification ──────────────────────────────────────────────────────

state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)

zone_codes  <- names(tourism)[grepl("^[A-G][A-Z]$", names(tourism))]
n_zones      <- length(zone_codes)

# Map each zone to its parent state index
zone_to_state <- match(substr(zone_codes, 1, 1), state_codes)

# Number of zones per state (for scaling)
zones_per_state <- as.integer(table(zone_to_state)[order(as.integer(names(table(zone_to_state))))])
zones_per_total <- n_zones

# ── Raw data ───────────────────────────────────────────────────────────────────

Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_zone  <- tourism %>% select(all_of(zone_codes))  %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

# ── Scaling ────────────────────────────────────────────────────────────────────
# Scale defined from state-level training data as before

Y_mean <- mean(Y_raw_zone[1:training_end, ])
Y_sd   <- sd(Y_raw_zone[1:training_end, ])

# Zones: default scaling
Y_zone <- (Y_raw_zone - Y_mean) / Y_sd

# States: scaled by zones_per_state
Y_state <- sweep(Y_raw_state, 2, zones_per_state, function(x, n) {
  (x - (Y_mean * n)) / Y_sd
})

# Total: scaled by n_zones
Y_total_scaled <- (Y_raw_total - Y_mean * n_zones) / Y_sd


# coherency and level penalties
lambda_total <- 1
lambda_state <- 1
sigma_coh <- 0.01 * n_zones

# ── Stan model ─────────────────────────────────────────────────────────────────
model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> n_zones;
  int<lower=1> H;
  matrix[T_train, n_state] Y_state;
  matrix[T_train, n_zones]  Y_zone;
  vector[T_train]           Y_total;
  array[n_zones] int<lower=1, upper=n_state> zone_to_state;
  
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
  vector<lower=0, upper=1>[n_zones]  alpha_zone;
  vector<lower=0, upper=1>[n_zones]  gamma_raw_zone;   // fraction of (1 - alpha)
  vector<lower=0>[n_zones]           sigma_zone;

  // ── Total-level parameters ────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;

// ── Initial states ────────────────────────────────────────────────────────
  vector[n_state]        l0_state;
  matrix[n_state, m-1]   s0_state_free;
  vector[n_zones]         l0_zone;
  matrix[n_zones, m-1]    s0_zone_free;
  real                   l0_total;
  vector[m-1]            s0_total_free;
}
transformed parameters {
  // ── Reconstruct full seasonals with sum-to-zero enforced ──────────────────
  matrix[n_state, m] s0_state;
  matrix[n_zones, m]  s0_zone;
  vector[m]          s0_total;

  for (n in 1:n_state) {
    s0_state[n, 1:(m-1)] = s0_state_free[n];
    s0_state[n, m]        = -sum(s0_state_free[n]);
  }
  for (n in 1:n_zones) {
    s0_zone[n, 1:(m-1)] = s0_zone_free[n];
    s0_zone[n, m]        = -sum(s0_zone_free[n]);
  }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);
  
  
  // ── Recover actual gamma values ───────────────────────────────────────────
  vector[n_state] gamma_state;
  vector[n_zones]  gamma_zone;
  real            gamma_total;

  for (n in 1:n_state)
    gamma_state[n] = gamma_raw_state[n] * (1 - alpha_state[n]);
  for (n in 1:n_zones)
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
  for (n in 1:n_zones) {
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
    for (z in 1:n_zones)
      if (zone_to_state[z] == n) n_zones_in_state += 1;
    l0_state[n]       ~ normal(0, n_zones_in_state);
    s0_state_free[n]  ~ normal(0, n_zones_in_state);
  }

  // Total: scaled by total number of zones
  l0_total      ~ normal(0, n_zones);
  s0_total_free ~ normal(0, n_zones);
  
  // -- constrain initial values
  for (k in 1:m){
    real diff_total = (l0_total + s0_total[k]) - (sum(l0_state) + sum(s0_state[,k]));
    real diff_middle = (sum(l0_state) + sum(s0_state[,k]) - (sum(l0_zone) + sum(s0_zone[,k])));
    
    target += normal_lpdf(diff_total | 0, sigma_coh);
    target += normal_lpdf(diff_middle | 0, sigma_coh);
  }

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_zones)
    for (i in 1:T_train)
      //Y_zone[i, n] ~ normal(l_z[n][i] + s_z[n][i], sigma_zone[n]);
      target += normal_lpdf(Y_zone[i,n] | l_z[n][i] + s_z[n][i], sigma_zone[n]);

  for (n in 1:n_state)
    for (i in 1:T_train)
      //Y_state[i, n] ~ normal(l_s[n][i] + s_s[n][i], sigma_state[n]);
      target += lambda_state * normal_lpdf(Y_state[i,n] | l_s[n][i] + s_s[n][i], sigma_state[n]);

  for (i in 1:T_train)
    //Y_total[i] ~ normal(l_t[i] + s_t[i], sigma_total);
    target += lambda_total * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);

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
  matrix[H, n_zones]  Y_fc_zone;
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
  T_train      = training_end,
  m            = m,
  n_state      = n_state,
  n_zones       = n_zones,
  H            = horizon + 1,
  Y_state      = Y_state[1:training_end, ],
  Y_zone       = Y_zone[1:training_end, ],
  Y_total      = Y_total_scaled[1:training_end],
  zone_to_state = zone_to_state, 
  lambda_total = lambda_total, 
  lambda_state = lambda_state, 
  sigma_coh    = sigma_coh
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
)
end_time <- Sys.time()
print(end_time - start_time)

# cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
# cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")
# 
# View(summary(stanFit)$summary)
# 
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
# gamma_pars <- all_pars[grepl("^gamma_state", all_pars)]
# pairs(stanFit, pars = gamma_pars)
# 
# gamma_means <- summary(stanFit)$summary[gamma_pars, "mean"]
# print(gamma_means)
# 
# # grab all l0-related parameters
# l0_pars <- all_pars[grepl("l0", all_pars)]
# 
# pairs(stanFit, pars = l0_pars)



# ── Accuracy ───────────────────────────────────────────────────────────────────
fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_zone_draws  <- rstan::extract(stanFit, pars = "Y_fc_zone")$Y_fc_zone
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

fc_state_means <- apply(fc_state_draws, c(2, 3), mean)[1:horizon, ]
fc_zone_means  <- apply(fc_zone_draws,  c(2, 3), mean)[1:horizon, ]
fc_total_means <- colMeans(fc_total_draws)[1:horizon]

# back-transform
fc_zone_orig  <- fc_zone_means * Y_sd + Y_mean
fc_state_orig <- sweep(fc_state_means, 2, zones_per_state, function(x, n) {
  (x * Y_sd) + (Y_mean * n)
})
fc_total_orig <- fc_total_means * Y_sd + (Y_mean * n_zones)


actuals_state <- Y_raw_state[(training_end + 1):(training_end + horizon), ]
actuals_zone  <- Y_raw_zone[ (training_end + 1):(training_end + horizon), ]
actuals_total <- Y_raw_total[(training_end + 1):(training_end + horizon)]

accuracy_state <- lapply(1:n_state, function(s) {
  errors <- fc_state_orig[, s] - actuals_state[, s]
  data.frame(series = state_codes[s], horizon = 1:horizon,
             actual = actuals_state[, s], forecast = fc_state_orig[, s],
             ae = abs(errors), se = errors^2)
}) %>% bind_rows()

accuracy_zone <- lapply(1:n_zones, function(z) {
  errors <- fc_zone_orig[, z] - actuals_zone[, z]
  data.frame(series = zone_codes[z], horizon = 1:horizon,
             actual = actuals_zone[, z], forecast = fc_zone_orig[, z],
             ae = abs(errors), se = errors^2)
}) %>% bind_rows()

total_errors <- fc_total_orig - actuals_total
accuracy_total <- data.frame(
  series = "Total", horizon = 1:horizon,
  actual = actuals_total, forecast = fc_total_orig,
  ae = abs(total_errors), se = total_errors^2
)

all_accuracy <- bind_rows(accuracy_state, accuracy_zone, accuracy_total)

summary_accuracy <- all_accuracy %>%
  group_by(series) %>%
  summarise(MAE = mean(ae), RMSE = sqrt(mean(se)), .groups = "drop") %>%
  bind_rows(
    all_accuracy %>%
      summarise(series = "Overall", MAE = mean(ae), RMSE = sqrt(mean(se)))
  )

print(summary_accuracy)
