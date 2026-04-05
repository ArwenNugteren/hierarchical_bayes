library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)
library(fable)
library(tsibble)
library(feasts)

# ── Data prep ─────────────────────────────────────────────────────────────────

load("tourism.Rdata")

training_end <- 96
m            <- 12L
horizon      <- 12

state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)

Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

Y_mean <- mean(Y_raw_state[1:training_end, ])
Y_sd   <- sd(Y_raw_state[1:training_end, ])

Y_state        <- (Y_raw_state - Y_mean) / Y_sd
Y_total_scaled <- (Y_raw_total - Y_mean * n_state) / (Y_sd * n_state)

# ── Fit fable ETS to extract initial values ───────────────────────────────────

# Helper: fit ETS(A,N,A) to a single numeric vector (training portion, scaled)
# and return list(l0, s0) in Stan's parameterisation.
#
# Stan's recursion initialises:
#   s_s[n][j] = s0_state[n, j]  for j in 1:m        (the oldest m seasonals)
#   l_s[n][1] = l0_state[n]
#
# fable stores seasonals as s_t for t = 1..T in components().
# At the start of the recursion (before observing any data) we need
# the m seasonal values that will be "lagged into" periods 1..m,
# i.e. the initial seasonal vector s_{1-m}, s_{2-m}, ..., s_0.
# These correspond to the *first* m rows of the components() seasonal column
# (before the recursion updates them). We use the fitted components at
# time 1..m as a good proxy (they are very close to the true initials for ETS).

fit_ets_inits <- function(y_scaled, m) {
  ts_df <- tsibble(
    time  = 1:training_end,
    value = y_scaled[1:training_end],
    index = time
  )
  
  fit <- ts_df %>%
    model(ets = ETS(value ~ error("A") + trend("N") + season("A", period = m)))
  
  comp <- components(fit)
  
  # Initial seasonal states are the first m rows (negative time indices)
  # fable orders them time = -(m-1), -(m-2), ..., 0  i.e. oldest to most recent
  s0 <- comp$season[1:m]   # already oldest-first, matches Stan's s0[1..m]
  
  # l0: the level at t=0 is the first non-NA level value (at time = 1, post first update)
  # Back-calculate: level[1] = alpha*(y[1] - s0_lag) + (1-alpha)*l0
  # Simpler and robust: just grab the first non-NA level as a close proxy,
  # OR extract from the fitted parameter list
  alpha <- fit %>% tidy() %>% filter(term == "alpha") %>% pull(estimate)
  level_t1 <- comp$level[comp$time == 1]
  s_lag_1  <- comp$season[comp$time == 0]   # seasonal that was lagged into t=1
  y_t1     <- comp$value[comp$time == 1]
  
  # Invert the ETS level update to recover l0:
  # level[1] = alpha*(y[1] - s_lag) + (1 - alpha)*l0
  # => l0 = (level[1] - alpha*(y[1] - s_lag)) / (1 - alpha)
  l0 <- (level_t1 - alpha * (y_t1 - s_lag_1)) / (1 - alpha)
  
  list(l0 = unname(l0), s0 = unname(s0))
}

# Fit for each state (scaled)
state_inits <- lapply(seq_len(n_state), function(n) {
  fit_ets_inits(Y_state[, n], m)
})

# Fit for total (scaled)
total_inits <- fit_ets_inits(Y_total_scaled, m)

# ── Build Stan init function ──────────────────────────────────────────────────

# Assemble l0_state (length n_state) and s0_state (n_state x m matrix)
l0_state_init <- sapply(state_inits, `[[`, "l0")
s0_state_init <- do.call(rbind, lapply(state_inits, `[[`, "s0"))  # n_state x m

l0_total_init <- total_inits$l0
s0_total_init <- total_inits$s0

make_init <- function() {
  list(
    # Hyperpriors — neutral starting points
    mu_alpha     = 0.3,
    mu_gamma_raw = 0.2,
    
    # State smoothing params — near the hyperprior means
    alpha_state     = rep(0.3, n_state),
    gamma_state_raw = rep(0.2, n_state),
    sigma_state     = rep(0.5, n_state),
    
    # Total smoothing params
    alpha_total     = 0.3,
    gamma_total_raw = 0.2,
    sigma_total     = 0.5,
    
    # ETS-derived initial states
    l0_state = l0_state_init,          # vector[n_state]
    s0_state = s0_state_init,          # matrix[n_state, m]  (row = state)
    l0_total = l0_total_init,          # scalar
    s0_total = s0_total_init           # vector[m]
  )
}

# ── Stan model ────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_state;
  int<lower=1> H;

  matrix[T_train, n_state] Y_state;
  vector[T_train]           Y_total;

  // ── ETS-derived initial states (fixed, not estimated) ─────────────────────
  vector[n_state]    l0_state;
  matrix[n_state, m] s0_state;
  real               l0_total;
  vector[m]          s0_total;
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
}
transformed parameters {
  // constraining gamma to be smaller than 1 - alpha
  vector[n_state] gamma_state;
  real gamma_total;
  for (n in 1:n_state){
    gamma_state[n] = gamma_state_raw[n] * (1 - alpha_state[n]);
  }
  gamma_total = gamma_total_raw * (1-alpha_total);
  
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
  mu_alpha ~ uniform(0,1);
  mu_gamma_raw ~ uniform(0,1);

  alpha_state ~ normal(mu_alpha, 0.1);
  gamma_state_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_state ~ cauchy(0, 0.3);

  alpha_total ~ normal(mu_alpha, 0.1);
  gamma_total_raw ~ normal(mu_gamma_raw, 0.1);
  sigma_total ~ cauchy(0, 0.3);

  // all l0/s0 priors and sum-to-zero constraints removed

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
}
"

# ── Data list ─────────────────────────────────────────────────────────────────
l0_state_init <- sapply(state_inits, `[[`, "l0")
s0_state_init <- do.call(rbind, lapply(state_inits, `[[`, "s0"))  # n_state x m

data_list <- list(
  T_train  = training_end,
  m        = m,
  n_state  = n_state,
  H        = horizon + 1,
  Y_state  = Y_state[1:training_end, ],
  Y_total  = Y_total_scaled[1:training_end],
  
  # Fixed initial states from fable
  l0_state = l0_state_init,
  s0_state = s0_state_init,
  l0_total = total_inits$l0,
  s0_total = total_inits$s0
)
# ── Compile and sample ────────────────────────────────────────────────────────

stan_model <- stan_model(model_code = model_string)

stanFit <- sampling(
  stan_model,
  data    = data_list,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,         # can reduce warmup since inits are good
  cores   = 4
  # control = list(adapt_delta = 0.95, max_treedepth = 12)
)


cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")


View(summary(stanFit)$summary)

# hyperparameters
pairs(stanFit, pars = c("mu_alpha", "mu_gamma_raw"))

# lower level parameters
# get all parameter names from the fit
all_pars <- names(stanFit)

# grab all alpha-related parameters
alpha_pars <- all_pars[grepl("^alpha_state", all_pars)]
pairs(stanFit, pars = alpha_pars)

alpha_means <- summary(stanFit)$summary[alpha_pars, "mean"]
print(alpha_means)

# grab all gamma-related parameters
gamma_pars <- all_pars[
  grepl("^gamma_state", all_pars) & 
    !grepl("gamma_state_raw", all_pars)
]
pairs(stanFit, pars = gamma_pars)

gamma_means <- summary(stanFit)$summary[gamma_pars, "mean"]
print(gamma_means)

# # grab all l0-related parameters
# l0_pars <- all_pars[grepl("l0", all_pars)]
# 
# pairs(stanFit, pars = l0_pars)



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



# After sampling, extract and plot the implied seasonal components
# Compare the magnitude of l vs s contributions to fitted values
fit_summary <- summary(stanFit)$summary
l_contrib <- fit_summary[grep("^l_s\\[", rownames(fit_summary)), "mean"]
s_contrib <- fit_summary[grep("^s_s\\[", rownames(fit_summary)), "mean"]

row_index_l <- as.integer(sub("l_s\\[(.*),.*\\]", "\\1", names(l_contrib)))
row_index_s <- as.integer(sub("s_s\\[(.*),.*\\]", "\\1", names(s_contrib)))


summary_by_row_l <- data.frame(
  value = l_contrib,
  i = row_index_l
) %>%
  group_by(i) %>%
  summarise(
    mean_contrib = mean(value),
    .groups = "drop"
  )

summary_by_row_s <- data.frame(
  value = s_contrib,
  i = row_index_s
) %>%
  group_by(i) %>%
  summarise(
    mean_contrib = mean(value),
    .groups = "drop"
  )



