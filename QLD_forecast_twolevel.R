library(fable)
library(dplyr)
library(stringr)
library(rstan)
library(posterior)

load("tourism.Rdata")

# start with Queensland (state C)
qld_tourism <- tourism %>%
  select(starts_with("C"))

training_start  <- 1                       
training_end    <- 96
m <- 12
horizon <- 12

# ── Hierarchy construction (automatic) ───────────────────────────────────────

# All series with only letter codes (no numeric suffix)
region_codes <- names(qld_tourism)[grepl("^C[A-Z]*$", names(qld_tourism))]

# Classify by code length
top_series    <- region_codes[nchar(region_codes) == 1]  # "C"
zone_series   <- region_codes[nchar(region_codes) == 2]  # CA, CB, CC, CD
region_series <- region_codes[nchar(region_codes) == 3]  # CAA, CAB, ...

n_top    <- length(top_series)
n_zone   <- length(zone_series)
n_region <- length(region_series)

# Map each region to its parent zone index
region_to_zone <- match(substr(region_series, 1, 2), zone_series)

# ── Scaling ───────────────────────────────────────────────────────────────────

# Number of regions per zone and total (for consistent aggregate scaling)
regions_per_zone  <- as.integer(table(region_to_zone)[order(as.integer(names(table(region_to_zone))))])
regions_per_total <- n_region

Y_raw_region <- qld_tourism %>% select(all_of(region_series)) %>% as.matrix()
Y_raw_zone   <- qld_tourism %>% select(all_of(zone_series))   %>% as.matrix()
Y_raw_top    <- qld_tourism %>% select(all_of(top_series))    %>% as.matrix()

# Common scale estimated from region-level data (same as wind power approach)
Y_mean <- mean(Y_raw_region[1:training_end, ])
Y_sd   <- sd(Y_raw_region[1:training_end, ])

# Regions: standard scaling
Y_region <- (Y_raw_region - Y_mean) / Y_sd

# Zones: scaled by (agg - n_children * Y_mean) / (n_children * Y_sd)
# so that a zone of n regions has the same per-region-unit scale
Y_zone <- sweep(Y_raw_zone, 2, Y_mean * regions_per_zone, "-")
Y_zone <- sweep(Y_zone,     2, Y_sd   * regions_per_zone, "/")

# Top: same logic with total region count
Y_top_scaled <- (Y_raw_top - Y_mean * regions_per_total) / (Y_sd * regions_per_total)

# ── Stan model ────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_region;
  int<lower=1> H;

  matrix[T_train, n_region] Y_region;
}
parameters {
  // ── Shared hyperparameters ────────────────────────────────────────────────────
  real mu_alpha;
  real mu_beta;
  real mu_gamma;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;

  // ── Series-level smoothing parameters ─────────────────────────────────────
  vector[n_region] alpha_region;
  vector[n_region] beta_region;
  vector[n_region] gamma_region;
  vector<lower=0>[n_region] sigma_region;

  // ── Initial states ────────────────────────────────────────────────────────
  vector[n_region] l0_region;
  vector[n_region] b0_region;
  matrix[n_region, m] s_raw_region;
}
transformed parameters {
  array[n_region] vector[T_train + 1] l_r;
  array[n_region] vector[T_train + 1] b_r;
  array[n_region] vector[T_train + m] s_r;

  for (n in 1:n_region) {
    vector[m] s0 = (s_raw_region[n])' - mean((s_raw_region[n])');
    for (j in 1:m) s_r[n][j] = s0[j];
    l_r[n][1] = l0_region[n];
    b_r[n][1] = b0_region[n];
    for (i in 1:T_train) {
      real s_lag   = s_r[n][i];
      l_r[n][i+1] = alpha_region[n] * (Y_region[i, n] - s_lag)
                  + (1 - alpha_region[n]) * (l_r[n][i] + b_r[n][i]);
      b_r[n][i+1] = beta_region[n]  * (l_r[n][i+1] - l_r[n][i])
                  + (1 - beta_region[n])  * b_r[n][i];
      s_r[n][i+m] = gamma_region[n] * (Y_region[i, n] - l_r[n][i] - b_r[n][i])
                  + (1 - gamma_region[n]) * s_lag;
    }
  }
}
model {
  // ── Hyperpriors ───────────────────────────────────────────────────────────
  mu_alpha ~ normal(0.5, 0.5);
  mu_beta  ~ normal(0.5, 0.5);
  mu_gamma ~ normal(0.5, 0.5);
  sigma_alpha ~ normal(0, 0.3);
  sigma_beta  ~ normal(0, 0.3);
  sigma_gamma ~ normal(0, 0.3);

  // ── Region-level priors drawn from shared hyperprior ──────────────────────
  alpha_region ~ normal(mu_alpha, sigma_alpha);
  beta_region  ~ normal(mu_beta,  sigma_beta);
  gamma_region ~ normal(mu_gamma, sigma_gamma);
  sigma_region ~ cauchy(0, 1);

  // ── Initial state priors ──────────────────────────────────────────────────
  l0_region ~ normal(0, 2);
  b0_region ~ normal(0, 0.5);
  for (n in 1:n_region) s_raw_region[n] ~ normal(0, 1);

  // ── Likelihood ────────────────────────────────────────────────────────────
  for (n in 1:n_region) {
    for (i in 1:T_train)
      Y_region[i, n] ~ normal(l_r[n][i] + b_r[n][i] + s_r[n][i], sigma_region[n]);
  }
}
generated quantities {
  matrix[H, n_region] Y_fc_region;

  for (n in 1:n_region) {
    vector[H] lf;
    vector[H] bf;
    vector[H] sf;
    lf[1] = l_r[n][T_train + 1];
    bf[1] = b_r[n][T_train + 1];
    for (h in 1:H) sf[h] = s_r[n][T_train + h - m];
    for (h in 1:(H - 1)) {
      Y_fc_region[h, n] = normal_rng(lf[h] + bf[h] + sf[h], sigma_region[n]);
      lf[h+1] = alpha_region[n] * (Y_fc_region[h, n] - sf[h])
              + (1 - alpha_region[n]) * (lf[h] + bf[h]);
      bf[h+1] = beta_region[n] * (lf[h+1] - lf[h])
              + (1 - beta_region[n]) * bf[h];
    }
    Y_fc_region[H, n] = normal_rng(lf[H] + bf[H] + sf[H], sigma_region[n]);
  }
}
"

data_list <- list(
  T_train  = training_end,
  m        = m,
  n_region = n_region,
  H        = horizon + 1,
  Y_region = Y_region[1:training_end, ]
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