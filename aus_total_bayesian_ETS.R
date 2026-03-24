library(fable)
library(dplyr)
library(stringr)
library(rstan)
library(posterior)
library(loo)
library(bayesplot)
library(ggplot2)

load("tourism.Rdata")

# start with australian total
aus_tourism <- tourism %>%
  select(Total)

## scale data
Y_mean <- mean(aus_tourism$Total)
Y_sd   <- sd(aus_tourism$Total)
Y_scaled <- (aus_tourism$Total - Y_mean) / Y_sd

training_start  <- 1                       
training_end    <- 96
T_full          <- length(Y_scaled)
n_forecast      <- T_full - training_end

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  vector[T_train] Y;
}
parameters {
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=1> gamma;
  real<lower=0>         sigma;
  real                  l_0;
  vector[m]             s_0;
}
transformed parameters {
  vector[T_train + m] s;
  vector[T_train + 1] l;

  for (j in 1:m) s[j] = s_0[j];
  l[1] = l_0;

  for (i in 1:T_train) {
    real s_lag = s[i];
    l[i+1] = alpha * (Y[i] - s_lag) + (1 - alpha) * l[i];
    s[i+m] = gamma * (Y[i] - l[i]) + (1 - gamma) * s_lag;
  }
}
model {
  alpha ~ uniform(0, 1);
  gamma ~ uniform(0, 1);
  sigma ~ cauchy(0, 1);
  l_0   ~ normal(Y[1], 1);
  
  s_0 ~ normal(0, 1);
  sum(s_0) ~ normal(0, 0.01 * m); // enforcing sum to 0 constraint for identifiability
  
  for (i in 1:T_train)
    Y[i] ~ normal(l[i] + s[i], sigma);
}
generated quantities {
  vector[13] l_forecast;
  vector[13] s_forecast;
  vector[13] Y_forecast;

  l_forecast[1] = l[T_train + 1];
  for (h in 1:13) s_forecast[h] = s[T_train + h - m];

  for (h in 1:12) {
    Y_forecast[h]   = normal_rng(l_forecast[h] + s_forecast[h], sigma);
    l_forecast[h+1] = alpha * (Y_forecast[h] - s_forecast[h])
                    + (1 - alpha) * l_forecast[h];
  }
  Y_forecast[13] = normal_rng(l_forecast[13] + s_forecast[13], sigma);

  real Y_pred_1  = Y_forecast[1];
  real Y_pred_2  = Y_forecast[2];
  real Y_pred_3  = Y_forecast[3];
  real Y_pred_6  = Y_forecast[6];
  real Y_pred_12 = Y_forecast[12];
  
  // calculate a log likelihood for AIC/BIC calculations
  vector[T_train] log_lik;
  for (i in 1:T_train) {
    log_lik[i] = normal_lpdf(Y[i] | l[i] + s[i], sigma);
  }
  // do some posterior predictive checks
  vector[T_train] Y_rep;
  for (i in 1:T_train) {
    Y_rep[i] = normal_rng(l[i] + s[i], sigma);
  }
}
"

data_list <- list(
  T_train = training_end,
  m       = 12L,
  Y       = Y_scaled[1:training_end]
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

alpha_post <- rstan::extract(stanFit, pars = "alpha")$alpha
hist(alpha_post)

gamma_post <- rstan::extract(stanFit, pars = "gamma")$gamma
hist(gamma_post)

sigma_post <- rstan::extract(stanFit, pars = "sigma")$sigma

cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")

pairs(stanFit,    pars = c("alpha", "gamma", "sigma"))
traceplot(stanFit, pars = c("alpha", "gamma"))

View(summary(stanFit)$summary)

pairs(stanFit, pars = c("alpha", "gamma", "sigma", "l_0"))
pairs(stanFit, pars = c("s_0[1]", "s_0[2]", "s_0[3]", "s_0[4]", "s_0[5]", "s_0[6]", 
                           "s_0[7]", "s_0[8]", "s_0[9]", "s_0[10]", "s_0[11]", "s_0[12]"))

# posterior predictive checks on in-sample training data
Y_rep_draws <- rstan::extract(stanFit, pars = "Y_rep")$Y_rep
# Y_rep_draws is [n_draws x T_train]
Y_obs <- Y_scaled[1:training_end]

# ── 1. Density overlay ────────────────────────────────────────────────────────
# Simulated densities (grey) overlaid on observed density (dark)
# Good for checking overall distributional shape
ppc_dens_overlay(Y_obs, Y_rep_draws[1:200, ])

# ── 2. Mean and interval bands over time ──────────────────────────────────────
# Checks whether the model captures the time series dynamics
ppc_ribbon(Y_obs, Y_rep_draws[1:200, ], x = 1:training_end)

# ── 3. Test statistics ────────────────────────────────────────────────────────
# Compare summary statistics of replicated vs observed data
# Ideally the observed value (dark line) sits well within the simulated distribution
ppc_stat(Y_obs, Y_rep_draws, stat = "mean")
ppc_stat(Y_obs, Y_rep_draws, stat = "sd")
ppc_stat(Y_obs, Y_rep_draws, stat = "max")
ppc_stat(Y_obs, Y_rep_draws, stat = "min")

# ── 4. Seasonal check ─────────────────────────────────────────────────────────
# Group by month to check if seasonal structure is captured
month_idx <- rep(1:12, length.out = training_end)
ppc_stat_grouped(Y_obs, Y_rep_draws, group = month_idx, stat = "mean")


# predictive checks over the forecast values
# Extract forecast draws [n_draws x 13], take first 12 horizons
Y_fc_draws <- rstan::extract(stanFit, pars = "Y_forecast")$Y_forecast[, 1:12]

# Back-transform to original scale
Y_fc_draws_orig <- Y_fc_draws * Y_sd + Y_mean

# Actual out-of-sample values
actuals <- aus_tourism$Total[(training_end + 1):(training_end + 12)]

# Month labels — continuing from wherever training ended
# training_end = 96 months from start, so forecast starts at month 97
last_train_month <- ((training_end - 1) %% 12) + 1   # which calendar month is obs 96
fc_months <- ((last_train_month:(last_train_month + 11)) %% 12) + 1

# Build a long data frame of draws for plotting
draws_df <- as.data.frame(Y_fc_draws_orig)
colnames(draws_df) <- paste0("h", 1:12)
draws_df$draw <- 1:nrow(draws_df)

draws_long <- draws_df %>%
  tidyr::pivot_longer(-draw, names_to = "horizon", values_to = "value") %>%
  mutate(horizon = as.integer(sub("h", "", horizon)),
         month   = fc_months[horizon])

actuals_df <- data.frame(
  horizon = 1:12,
  month   = fc_months,
  actual  = actuals
)

ggplot(draws_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", alpha = 0.4) +
  geom_vline(data = actuals_df, aes(xintercept = actual),
             colour = "navy", linewidth = 1) +
  facet_wrap(~ horizon, scales = "free",
             labeller = labeller(horizon = setNames(
               paste0("h=", 1:12, "  (month ", fc_months, ")"),
               1:12
             ))) +
  labs(x = NULL, y = NULL,
       title = "Forecast posterior predictive vs actuals",
       subtitle = "Dark line = observed value") +
  theme_minimal()



#### accuracy evaluation 
# Extract posterior means of forecasts
fc_draws <- rstan::extract(stanFit, pars = "Y_forecast")$Y_forecast
# fc_draws is [n_draws x 13], we want horizons 1:12
fc_means <- colMeans(fc_draws)[1:12]

# Back-transform to original scale
fc_means_orig <- fc_means * Y_sd + Y_mean

# Actual values: the 12 observations immediately after training_end
actuals <- aus_tourism$Total[(training_end + 1):(training_end + 12)]

# Compute errors
errors    <- fc_means_orig - actuals
ae        <- abs(errors)
se        <- errors^2

mae  <- mean(ae)
rmse <- sqrt(mean(se))

# Per-horizon accuracy table
accuracy_df <- data.frame(
  horizon  = 1:12,
  actual   = actuals,
  forecast = fc_means_orig,
  error    = errors,
  abs_error = ae,
  sq_error  = se
)

cat(sprintf("MAE:  %.4f\n", mae))
cat(sprintf("RMSE: %.4f\n", rmse))
print(accuracy_df)



#### calculating AIC, AICc and BIC on the mean of the log-likelihood
log_lik      <- loo::extract_log_lik(stanFit, parameter_name = "log_lik")
mean_log_lik_original <- sum(colMeans(log_lik)) - training_end * log(Y_sd)


# Parameter count: alpha, gamma, sigma, l_0, s_raw[1:12] = 16
k <- 16
n <- training_end

aic  <- -2 * mean_log_lik_original + 2 * k
aicc <- aic + (2 * k^2 + 2 * k) / (n - k - 1)
bic  <- -2 * mean_log_lik_original + k * log(n)

cat(sprintf("AIC:  %.4f\n", aic))
cat(sprintf("AICc: %.4f\n", aicc))
cat(sprintf("BIC:  %.4f\n", bic))
