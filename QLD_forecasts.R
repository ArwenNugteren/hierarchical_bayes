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
  int<lower=1> n_zone;
  int H;

  matrix[T_train, n_region] Y_region;
  matrix[T_train, n_zone]   Y_zone;
  vector[T_train]            Y_top;

  array[n_region] int<lower=1, upper=n_zone> region_to_zone;
}
parameters {
  // ── Top-level hyperpriors ──────────────────────────────────────────────────
  real mu_alpha_top;
  real mu_beta_top;
  real mu_gamma_top;
  real<lower=0> sigma_alpha_top;
  real<lower=0> sigma_beta_top;
  real<lower=0> sigma_gamma_top;

  // ── Zone-level hyperpriors (one per zone, drawn from top) ─────────────────
  vector[n_zone] mu_alpha_zone;
  vector[n_zone] mu_beta_zone;
  vector[n_zone] mu_gamma_zone;
  real<lower=0> sigma_alpha_zone;   // shared spread at zone level
  real<lower=0> sigma_beta_zone;
  real<lower=0> sigma_gamma_zone;

  // ── Series-level smoothing parameters ─────────────────────────────────────
  // Region level — drawn from parent zone hyperprior
  vector[n_region] alpha_region;
  vector[n_region] beta_region;
  vector[n_region] gamma_region;
  vector<lower=0>[n_region] sigma_region;

  // Zone level — drawn from top hyperprior
  vector[n_zone] alpha_zone;
  vector[n_zone] beta_zone;
  vector[n_zone] gamma_zone;
  vector<lower=0>[n_zone] sigma_zone;

  // Top level — drawn from top hyperprior
  real alpha_top;
  real beta_top;
  real gamma_top;
  real<lower=0> sigma_top;

  // ── Initial states (per series) ───────────────────────────────────────────
  vector[n_region] l0_region;
  vector[n_region] b0_region;
  matrix[n_region, m] s_raw_region;

  vector[n_zone] l0_zone;
  vector[n_zone] b0_zone;
  matrix[n_zone, m] s_raw_zone;

  real l0_top;
  real b0_top;
  vector[m] s_raw_top;
}
transformed parameters {
  // ── State recursions — region ─────────────────────────────────────────────
  array[n_region] vector[T_train + 1] l_r;
  array[n_region] vector[T_train + 1] b_r;
  array[n_region] vector[T_train + m] s_r;

  for (n in 1:n_region) {
    vector[m] s0 = (s_raw_region[n])' - mean((s_raw_region[n])');
    for (j in 1:m) s_r[n][j] = s0[j];
    l_r[n][1] = l0_region[n];
    b_r[n][1] = b0_region[n];
    for (i in 1:T_train) {
      real s_lag     = s_r[n][i];
      l_r[n][i+1]   = alpha_region[n] * (Y_region[i, n] - s_lag)
                     + (1 - alpha_region[n]) * (l_r[n][i] + b_r[n][i]);
      b_r[n][i+1]   = beta_region[n]  * (l_r[n][i+1] - l_r[n][i])
                     + (1 - beta_region[n])  * b_r[n][i];
      s_r[n][i+m]   = gamma_region[n] * (Y_region[i, n] - l_r[n][i] - b_r[n][i])
                     + (1 - gamma_region[n]) * s_lag;
    }
  }

  // ── State recursions — zone ───────────────────────────────────────────────
  array[n_zone] vector[T_train + 1] l_z;
  array[n_zone] vector[T_train + 1] b_z;
  array[n_zone] vector[T_train + m] s_z;

  for (n in 1:n_zone) {
    vector[m] s0 = (s_raw_zone[n])' - mean((s_raw_zone[n])');
    for (j in 1:m) s_z[n][j] = s0[j];
    l_z[n][1] = l0_zone[n];
    b_z[n][1] = b0_zone[n];
    for (i in 1:T_train) {
      real s_lag     = s_z[n][i];
      l_z[n][i+1]   = alpha_zone[n] * (Y_zone[i, n] - s_lag)
                     + (1 - alpha_zone[n]) * (l_z[n][i] + b_z[n][i]);
      b_z[n][i+1]   = beta_zone[n]  * (l_z[n][i+1] - l_z[n][i])
                     + (1 - beta_zone[n])  * b_z[n][i];
      s_z[n][i+m]   = gamma_zone[n] * (Y_zone[i, n] - l_z[n][i] - b_z[n][i])
                     + (1 - gamma_zone[n]) * s_lag;
    }
  }

  // ── State recursions — top ────────────────────────────────────────────────
  vector[T_train + 1] l_t;
  vector[T_train + 1] b_t;
  vector[T_train + m] s_t;

  {
    vector[m] s0 = s_raw_top - mean(s_raw_top);
    for (j in 1:m) s_t[j] = s0[j];
    l_t[1] = l0_top;
    b_t[1] = b0_top;
    for (i in 1:T_train) {
      real s_lag  = s_t[i];
      l_t[i+1]   = alpha_top * (Y_top[i] - s_lag)
                 + (1 - alpha_top) * (l_t[i] + b_t[i]);
      b_t[i+1]   = beta_top  * (l_t[i+1] - l_t[i])
                 + (1 - beta_top)  * b_t[i];
      s_t[i+m]   = gamma_top * (Y_top[i] - l_t[i] - b_t[i])
                 + (1 - gamma_top) * s_lag;
    }
  }
}
model {
  // ── Top-level hyperpriors ─────────────────────────────────────────────────
  mu_alpha_top ~ normal(0.5, 0.5);
  mu_beta_top  ~ normal(0.5, 0.5);
  mu_gamma_top ~ normal(0.5, 0.5);
  sigma_alpha_top ~ normal(0, 0.3);
  sigma_beta_top  ~ normal(0, 0.3);
  sigma_gamma_top ~ normal(0, 0.3);

  // ── Zone hyperpriors drawn from top ──────────────────────────────────────
  mu_alpha_zone ~ normal(mu_alpha_top, sigma_alpha_top);
  mu_beta_zone  ~ normal(mu_beta_top,  sigma_beta_top);
  mu_gamma_zone ~ normal(mu_gamma_top, sigma_gamma_top);
  sigma_alpha_zone ~ normal(0, 0.3);
  sigma_beta_zone  ~ normal(0, 0.3);
  sigma_gamma_zone ~ normal(0, 0.3);

  // ── Region-level params drawn from parent zone ────────────────────────────
  for (n in 1:n_region) {
    int z = region_to_zone[n];
    alpha_region[n] ~ normal(mu_alpha_zone[z], sigma_alpha_zone);
    beta_region[n]  ~ normal(mu_beta_zone[z],  sigma_beta_zone);
    gamma_region[n] ~ normal(mu_gamma_zone[z], sigma_gamma_zone);
  }
  sigma_region ~ cauchy(0, 1);

  // ── Zone-level params drawn from top ─────────────────────────────────────
  alpha_zone ~ normal(mu_alpha_top, sigma_alpha_top);
  beta_zone  ~ normal(mu_beta_top,  sigma_beta_top);
  gamma_zone ~ normal(mu_gamma_top, sigma_gamma_top);
  sigma_zone ~ cauchy(0, 1);

  // ── Top-level series params drawn from top hyperprior ─────────────────────
  alpha_top ~ normal(mu_alpha_top, sigma_alpha_top);
  beta_top  ~ normal(mu_beta_top,  sigma_beta_top);
  gamma_top ~ normal(mu_gamma_top, sigma_gamma_top);
  sigma_top ~ cauchy(0, 1);

  // ── Initial state priors ──────────────────────────────────────────────────
  l0_region ~ normal(0, 2);
  b0_region ~ normal(0, 0.5);
  l0_zone   ~ normal(0, 2);
  b0_zone   ~ normal(0, 0.5);
  l0_top    ~ normal(0, 2);
  b0_top    ~ normal(0, 0.5);
  for (n in 1:n_region) s_raw_region[n] ~ normal(0, 1);
  for (n in 1:n_zone)   s_raw_zone[n]   ~ normal(0, 1);
  s_raw_top ~ normal(0, 1);

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_region) {
    for (i in 1:T_train)
      Y_region[i, n] ~ normal(l_r[n][i] + b_r[n][i] + s_r[n][i], sigma_region[n]);
  }
  for (n in 1:n_zone) {
    for (i in 1:T_train)
      Y_zone[i, n] ~ normal(l_z[n][i] + b_z[n][i] + s_z[n][i], sigma_zone[n]);
  }
  for (i in 1:T_train)
    Y_top[i] ~ normal(l_t[i] + b_t[i] + s_t[i], sigma_top);
}
generated quantities {
  matrix[H, n_region] Y_fc_region;
  matrix[H, n_zone]   Y_fc_zone;
  vector[H]           Y_fc_top;

  for (n in 1:n_region) {
    vector[H] lf; vector[H] bf; vector[H] sf;
    lf[1] = l_r[n][T_train + 1];
    bf[1] = b_r[n][T_train + 1];
    for (h in 1:H) sf[h] = s_r[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_region[h, n] = normal_rng(lf[h] + bf[h] + sf[h], sigma_region[n]);
      lf[h+1] = alpha_region[n] * (Y_fc_region[h, n] - sf[h])
              + (1 - alpha_region[n]) * (lf[h] + bf[h]);
      bf[h+1] = beta_region[n]  * (lf[h+1] - lf[h])
              + (1 - beta_region[n]) * bf[h];
    }
    Y_fc_region[H, n] = normal_rng(lf[H] + bf[H] + sf[H], sigma_region[n]);
  }

  for (n in 1:n_zone) {
    vector[H] lf; vector[H] bf; vector[H] sf;
    lf[1] = l_z[n][T_train + 1];
    bf[1] = b_z[n][T_train + 1];
    for (h in 1:H) sf[h] = s_z[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_zone[h, n] = normal_rng(lf[h] + bf[h] + sf[h], sigma_zone[n]);
      lf[h+1] = alpha_zone[n] * (Y_fc_zone[h, n] - sf[h])
              + (1 - alpha_zone[n]) * (lf[h] + bf[h]);
      bf[h+1] = beta_zone[n]  * (lf[h+1] - lf[h])
              + (1 - beta_zone[n]) * bf[h];
    }
    Y_fc_zone[H, n] = normal_rng(lf[H] + bf[H] + sf[H], sigma_zone[n]);
  }

  {
    vector[H] lf; vector[H] bf; vector[H] sf;
    lf[1] = l_t[T_train + 1];
    bf[1] = b_t[T_train + 1];
    for (h in 1:H) sf[h] = s_t[T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_top[h] = normal_rng(lf[h] + bf[h] + sf[h], sigma_top);
      lf[h+1] = alpha_top * (Y_fc_top[h] - sf[h])
              + (1 - alpha_top) * (lf[h] + bf[h]);
      bf[h+1] = beta_top  * (lf[h+1] - lf[h])
              + (1 - beta_top)  * bf[h];
    }
    Y_fc_top[H] = normal_rng(lf[H] + bf[H] + sf[H], sigma_top);
  }
}
"

# ── Data list ─────────────────────────────────────────────────────────────────

data_list <- list(
  T_train       = training_end,
  m             = m,
  n_region      = n_region,
  n_zone        = n_zone,
  Y_region      = Y_region[1:training_end, ],
  Y_zone        = Y_zone[1:training_end, ],
  Y_top         = Y_top_scaled[1:training_end, 1],
  region_to_zone = region_to_zone, 
  H = horizon + 1
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