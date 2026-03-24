# ── Data prep ─────────────────────────────────────────────────────────────────

load("tourism.Rdata")

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
Y_total_scaled <- (Y_raw_total - Y_mean * n_state) / (Y_sd * n_state)

training_end <- 96
m            <- 12L
horizon      <- 12

# ── Stan model ────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> H;

  matrix[T_train, n_state] Y_state;
  vector[T_train]           Y_total;
}
parameters {
  // ── Hyperpriors on alpha and gamma ────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_gamma;

  // ── State-level smoothing parameters ─────────────────────────────────────
  vector<lower=0, upper=1>[n_state] alpha_state;
  vector<lower=0, upper=1>[n_state] gamma_state;
  vector<lower=0>[n_state]          sigma_state;

  // ── Total-level smoothing parameters (drawn from same hyperprior) ─────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_total;
  real<lower=0>          sigma_total;

  // ── Initial states — state level ──────────────────────────────────────────
  vector[n_state] l0_state;
  matrix[n_state, m] s0_state;

  // ── Initial states — total level ──────────────────────────────────────────
  real           l0_total;
  vector[m]      s0_total;
}
transformed parameters {
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
  mu_alpha ~ uniform(0, 1);
  mu_gamma ~ uniform(0, 1);
  sigma_alpha ~ normal(0, 0.3);
  sigma_gamma ~ normal(0, 0.3);

  // ── State-level priors drawn from hyperprior ──────────────────────────────
  alpha_state ~ normal(mu_alpha, sigma_alpha);
  gamma_state ~ normal(mu_gamma, sigma_gamma);
  sigma_state ~ cauchy(0, 1);

  // ── Total-level priors drawn from same hyperprior ─────────────────────────
  alpha_total ~ normal(mu_alpha, sigma_alpha);
  gamma_total ~ normal(mu_gamma, sigma_gamma);
  sigma_total ~ cauchy(0, 1);

  // ── Initial state priors ──────────────────────────────────────────────────
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
    for (i in 1:T_train)
      Y_state[i, n] ~ normal(l_s[n][i] + s_s[n][i], sigma_state[n]);
  }
  for (i in 1:T_train)
    Y_total[i] ~ normal(l_t[i] + s_t[i], sigma_total);
}
generated quantities {
  matrix[H, n_state] Y_fc_state;
  vector[H]          Y_fc_total;

  // ── State forecasts ───────────────────────────────────────────────────────
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

  // ── Total forecasts ───────────────────────────────────────────────────────
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

  // ── Log likelihood ────────────────────────────────────────────────────────
  matrix[T_train, n_state] log_lik_state;
  vector[T_train]           log_lik_total;
  for (n in 1:n_state)
    for (i in 1:T_train)
      log_lik_state[i, n] = normal_lpdf(Y_state[i, n] | l_s[n][i] + s_s[n][i], sigma_state[n]);
  for (i in 1:T_train)
    log_lik_total[i] = normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
}
"

# ── Data list ─────────────────────────────────────────────────────────────────

data_list <- list(
  T_train  = training_end,
  m        = m,
  n_state  = n_state,
  H        = horizon + 1,
  Y_state  = Y_state[1:training_end, ],
  Y_total  = Y_total_scaled[1:training_end]
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

View(summary(stanFit)$summary)

# hyperparameters
pairs(stanFit, pars = c("mu_alpha", "mu_gamma"))

# lower level parameters
# get all parameter names from the fit
all_pars <- names(stanFit)

# grab all alpha-related parameters
alpha_pars <- all_pars[grepl("alpha", all_pars)]

pairs(stanFit, pars = alpha_pars)

# grab all gamma-related parameters
gamma_pars <- all_pars[grepl("gamma", all_pars)]

pairs(stanFit, pars = gamma_pars)

# grab all gamma-related parameters
l0_pars <- all_pars[grepl("l0", all_pars)]

pairs(stanFit, pars = l0_pars)


# model fit evaluations
k <- 4 + 3*n_state + 3 + (n_state + 1) + (n_state + 1)*m
n <- training_end

# Extract log likelihoods
log_lik_state <- loo::extract_log_lik(stanFit, parameter_name = "log_lik_state")
log_lik_total <- loo::extract_log_lik(stanFit, parameter_name = "log_lik_total")

# Sum pointwise posterior means across all series and time points
mean_log_lik <- sum(colMeans(log_lik_state)) + sum(colMeans(log_lik_total))

# Jacobian adjustment to get back to original scale
# One adjustment per series (n_state state series + 1 total)
mean_log_lik_orig <- mean_log_lik - (n_state + 1) * training_end * log(Y_sd)

n <- training_end * (n_state + 1)  # 96 * 8 = 768

aic  <- -2 * mean_log_lik_orig + 2 * k
aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1)
bic  <- -2 * mean_log_lik_orig + k * log(n)

cat(sprintf("n:    %d\n",   n))
cat(sprintf("k:    %d\n",   k))
cat(sprintf("AIC:  %.4f\n", aic))
cat(sprintf("AICc: %.4f\n", aicc))
cat(sprintf("BIC:  %.4f\n", bic))


## accuracy calculation 
fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total
# fc_state_draws is [n_draws x H x n_state]
# fc_total_draws is [n_draws x H]

# Take posterior means, drop the last horizon (H = horizon + 1)
fc_state_means <- apply(fc_state_draws, c(2, 3), mean)[1:horizon, ]  # [horizon x n_state]
fc_total_means <- colMeans(fc_total_draws)[1:horizon]                 # [horizon]

# ── Back-transform to original scale ──────────────────────────────────────────

fc_state_orig <- fc_state_means * Y_sd + Y_mean
fc_total_orig <- fc_total_means * Y_sd * n_state + Y_mean * n_state

# ── Actuals ────────────────────────────────────────────────────────────────────

actuals_state <- Y_raw_state[(training_end + 1):(training_end + horizon), ]  # [horizon x n_state]
actuals_total <- Y_raw_total[(training_end + 1):(training_end + horizon)]    # [horizon]

# ── Per-series accuracy ────────────────────────────────────────────────────────

state_accuracy <- lapply(1:n_state, function(s) {
  errors <- fc_state_orig[, s] - actuals_state[, s]
  data.frame(
    series  = state_codes[s],
    horizon = 1:horizon,
    actual  = actuals_state[, s],
    forecast = fc_state_orig[, s],
    error   = errors,
    mae     = mean(abs(errors)),
    rmse    = sqrt(mean(errors^2))
  )
}) %>% bind_rows()

# Total series accuracy
total_errors <- fc_total_orig - actuals_total
total_accuracy <- data.frame(
  series   = "Total",
  horizon  = 1:horizon,
  actual   = actuals_total,
  forecast = fc_total_orig,
  error    = total_errors,
  mae      = mean(abs(total_errors)),
  rmse     = sqrt(mean(total_errors^2))
)

# ── Summary table ──────────────────────────────────────────────────────────────

summary_accuracy <- bind_rows(state_accuracy, total_accuracy) %>%
  group_by(series) %>%
  summarise(MAE = first(mae), RMSE = first(rmse), .groups = "drop") %>%
  bind_rows(
    bind_rows(state_accuracy, total_accuracy) %>%
      summarise(series = "Overall", 
                MAE  = mean(abs(error)), 
                RMSE = sqrt(mean(error^2)))
  )

print(summary_accuracy)
