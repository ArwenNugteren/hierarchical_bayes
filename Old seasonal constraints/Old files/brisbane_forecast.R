library(fable)
library(dplyr)
library(stringr)
library(rstan)
library(posterior)

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
T_full          <- length(Y_scaled)         # full dataset length
n_forecast      <- T_full - training_end    # number of one-step ahead forecasts


### stan model
model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;              // seasonal period (12 for monthly)
  vector[T_train] Y;
}
parameters {
  real<lower=0,upper=1> alpha;
  real<lower=0,upper=1> beta;
  real<lower=0,upper=1> gamma;
  real<lower=0>         sigma;
  real                  l_0;
  real                  b_0;
  vector[m]             s_raw;  // unconstrained seasonal seeds
}
transformed parameters {
  // Enforce sum-to-zero on initial seasonals
  vector[m] s_0;
  s_0 = s_raw - mean(s_raw);
  // State vectors (1-indexed: index 1 = initial state)
  vector[T_train + 1] l;
  vector[T_train + 1] b;
  // Seasonal indices: need m lags before training + T_train updates
  vector[T_train + m] s;
  // Seed the first m seasonal slots from s_0
  for (j in 1:m)
    s[j] = s_0[j];
  l[1] = l_0;
  b[1] = b_0;
  for (i in 1:T_train) {
    real s_lag = s[i];                            // s[t - m]
    l[i+1] = alpha * (Y[i] - s_lag)
             + (1 - alpha) * (l[i] + b[i]);
    b[i+1] = beta  * (l[i+1] - l[i])
             + (1 - beta)  * b[i];
    s[i+m] = gamma * (Y[i] - l[i] - b[i])
             + (1 - gamma) * s_lag;
  }
}
model {
  alpha ~ normal(0.5, 0.5);
  beta  ~ normal(0.5, 0.5);           // shrink toward zero — trend often slow-moving
  gamma ~ normal(0.5, 0.5);           // shrink toward zero — seasonality often stable
  sigma ~ cauchy(0, 1);
  l_0   ~ normal(Y[1], 1);
  b_0   ~ normal(0, 0.5);
  s_raw ~ normal(0, 1);
  for (i in 1:T_train) {
    Y[i] ~ normal(l[i] + b[i] + s[i], sigma);
  }
}
generated quantities {
  vector[13] l_forecast;
  vector[13] b_forecast;
  vector[13] s_forecast;
  vector[13] Y_forecast;
  l_forecast[1] = l[T_train + 1];
  b_forecast[1] = b[T_train + 1];
  // Carry forward the last m seasonal indices for future lags
  // h=1 needs s[T_train + 1], which sits at index T_train + 1 in s[]
  for (h in 1:13) {
    // seasonal lag: wrap back m steps into the estimated s vector
    int s_idx = T_train + h - m;     // index into s[] for the seasonal lag
    s_forecast[h] = s[s_idx];
  }
  // Simulate forward
  for (h in 1:12) {
    Y_forecast[h]     = normal_rng(l_forecast[h] + b_forecast[h] + s_forecast[h], sigma);
    l_forecast[h+1]   = alpha * (Y_forecast[h] - s_forecast[h])
                        + (1 - alpha) * (l_forecast[h] + b_forecast[h]);
    b_forecast[h+1]   = beta  * (l_forecast[h+1] - l_forecast[h])
                        + (1 - beta)  * b_forecast[h];
  }
  Y_forecast[13] = normal_rng(l_forecast[13] + b_forecast[13] + s_forecast[13], sigma);
  real Y_pred_1  = Y_forecast[1];
  real Y_pred_2  = Y_forecast[2];
  real Y_pred_3  = Y_forecast[3];
  real Y_pred_6  = Y_forecast[6];
  real Y_pred_12 = Y_forecast[12];
}
"
# Updated data list — add m = 12
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
  data = data_list, 
  seed = save_seed,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4
)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

alpha_post <- rstan::extract(stanFit, pars="alpha")$alpha
hist(alpha_post)

beta_post <- rstan::extract(stanFit, pars="beta")$beta
hist(beta_post)

gamma_post <- rstan::extract(stanFit, pars="gamma")$gamma
hist(gamma_post)

max(summary(stanFit)$summary[,"Rhat"])

min(summary(stanFit)$summary[,"n_eff"], na.rm=T)

View(summary(stanFit)$summary)

pairs(stanFit, pars=c("alpha", "beta", "gamma"))
traceplot(stanFit, pars=c("alpha", "beta", "gamma"))
