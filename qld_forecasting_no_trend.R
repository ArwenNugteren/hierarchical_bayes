library(fable)
library(dplyr)
library(stringr)
library(rstan)
library(posterior)
library(loo)

load("tourism.Rdata")

# start with Queensland (state C)
qld_tourism <- tourism %>%
  select(starts_with("C"))

brisbane_tourism <- qld_tourism %>%
  select(C)

## scale data
Y_mean <- mean(brisbane_tourism$C)
Y_sd   <- sd(brisbane_tourism$C)
Y_scaled <- (brisbane_tourism$C - Y_mean) / Y_sd

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
  vector[m]             s_raw;
}
transformed parameters {
  vector[m] s_0;
  s_0 = s_raw - mean(s_raw);

  vector[T_train + 1] l;
  vector[T_train + m] s;

  for (j in 1:m) s[j] = s_0[j];
  l[1] = l_0;

  for (i in 1:T_train) {
    real s_lag = s[i];
    l[i+1] = alpha * (Y[i] - s_lag) + (1 - alpha) * l[i];
    s[i+m] = gamma * (Y[i] - l[i]) + (1 - gamma) * s_lag;
  }
}
model {
  alpha ~ normal(0.5, 0.5);
  gamma ~ normal(0.5, 0.5);
  sigma ~ cauchy(0, 1);
  l_0   ~ normal(Y[1], 1);
  s_raw ~ normal(0, 1);

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


#### accuracy evaluation 
# Extract posterior means of forecasts
fc_draws <- rstan::extract(stanFit, pars = "Y_forecast")$Y_forecast
# fc_draws is [n_draws x 13], we want horizons 1:12
fc_means <- colMeans(fc_draws)[1:12]

# Back-transform to original scale
fc_means_orig <- fc_means * Y_sd + Y_mean

# Actual values: the 12 observations immediately after training_end
actuals <- brisbane_tourism$C[(training_end + 1):(training_end + 12)]

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
