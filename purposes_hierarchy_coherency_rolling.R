library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

load("tourism.Rdata")

# ── Read job name and directories ──────────────────────────────────────────────

args     <- commandArgs(trailingOnly = TRUE)
job_name <- args[1]

dir_accuracy <- file.path(job_name, "accuracy")
dir_stanfits <- file.path(job_name, "stanfits")

# ── Series identification ──────────────────────────────────────────────────────

all_names <- names(tourism)

rp_codes <- all_names[grepl("^[A-G][A-Z][A-Z](Hol|Vis|Bus|Oth)$", all_names)]
zp_codes <- all_names[grepl("^[A-G][A-Z](Hol|Vis|Bus|Oth)$",      all_names)]
sp_codes <- all_names[grepl("^[A-G](Hol|Vis|Bus|Oth)$",           all_names)]
pp_codes <- all_names[grepl("^Total(Hol|Vis|Bus|Oth)$",           all_names)]

n_rp <- length(rp_codes)
n_zp <- length(zp_codes)
n_sp <- length(sp_codes)
n_pp <- length(pp_codes)

# ── Hierarchy mappings ─────────────────────────────────────────────────────────

rp_to_zp <- match(paste0(substr(rp_codes, 1, 2),
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), zp_codes)
rp_to_sp <- match(paste0(substr(rp_codes, 1, 1),
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
rp_to_pp <- match(paste0("Total",
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)
zp_to_sp <- match(paste0(substr(zp_codes, 1, 1),
                         str_extract(zp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
zp_to_pp <- match(paste0("Total",
                         str_extract(zp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)
sp_to_pp <- match(paste0("Total",
                         str_extract(sp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)

rp_per_zp <- as.integer(table(rp_to_zp)[order(as.integer(names(table(rp_to_zp))))])
rp_per_sp <- as.integer(table(rp_to_sp)[order(as.integer(names(table(rp_to_sp))))])
rp_per_pp <- as.integer(table(rp_to_pp)[order(as.integer(names(table(rp_to_pp))))])

# ── Raw data ───────────────────────────────────────────────────────────────────

Y_raw_rp    <- tourism %>% select(all_of(rp_codes)) %>% as.matrix()
Y_raw_zp    <- tourism %>% select(all_of(zp_codes)) %>% as.matrix()
Y_raw_sp    <- tourism %>% select(all_of(sp_codes)) %>% as.matrix()
Y_raw_pp    <- tourism %>% select(all_of(pp_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

# ── Rolling window setup ───────────────────────────────────────────────────────

m           <- 12L
window_size <- 96
horizon     <- 12
test_set    <- 121
origins     <- 96:(96 + test_set - 1)

# ── Get this task's origin from PBS array index ────────────────────────────────

task_id <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
origin  <- origins[task_id]

# ── Training and test windows ──────────────────────────────────────────────────

train_idx <- (origin - window_size + 1):origin
test_idx  <- (origin + 1):(origin + horizon)

# ── Scaling (rp-based, from training window) ───────────────────────────────────

Y_mean <- mean(Y_raw_rp[train_idx, ])
Y_sd   <- sd(Y_raw_rp[train_idx, ])

Y_rp <- (Y_raw_rp - Y_mean) / Y_sd

Y_zp <- sweep(Y_raw_zp, 2, Y_mean * rp_per_zp, "-")
Y_zp <- sweep(Y_zp,     2, Y_sd,               "/")

Y_sp <- sweep(Y_raw_sp, 2, Y_mean * rp_per_sp, "-")
Y_sp <- sweep(Y_sp,     2, Y_sd,               "/")

Y_pp <- sweep(Y_raw_pp, 2, Y_mean * rp_per_pp, "-")
Y_pp <- sweep(Y_pp,     2, Y_sd,               "/")

Y_total <- (Y_raw_total - Y_mean * n_rp) / Y_sd


# weights and penalties
lambda_total <- 1
lambda_pp    <- 1
lambda_sp    <- 1
lambda_zp    <- 1
sigma_coh    <- 0.001 * n_rp

# ── Stan model ─────────────────────────────────────────────────────────────────

model_string <- "
data {
  int<lower=1> T_train;
  int<lower=1> m;
  int<lower=1> n_rp;
  int<lower=1> n_zp;
  int<lower=1> n_sp;
  int<lower=1> n_pp;
  int<lower=1> H;

  matrix[T_train, n_rp] Y_rp;
  matrix[T_train, n_zp] Y_zp;
  matrix[T_train, n_sp] Y_sp;
  matrix[T_train, n_pp] Y_pp;
  vector[T_train]        Y_total;

  array[n_rp] int<lower=1, upper=n_zp> rp_to_zp;
  array[n_rp] int<lower=1, upper=n_sp> rp_to_sp;
  array[n_rp] int<lower=1, upper=n_pp> rp_to_pp;
  array[n_zp] int<lower=1, upper=n_sp> zp_to_sp;
  array[n_zp] int<lower=1, upper=n_pp> zp_to_pp;
  array[n_sp] int<lower=1, upper=n_pp> sp_to_pp;

  array[n_zp] int<lower=1> rp_per_zp;
  array[n_sp] int<lower=1> rp_per_sp;
  array[n_pp] int<lower=1> rp_per_pp;

  real<lower=0> lambda_total;
  real<lower=0> lambda_pp;
  real<lower=0> lambda_sp;
  real<lower=0> lambda_zp;
  real<lower=0> sigma_coh;
}
parameters {
  // ── Global hyperpriors ────────────────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;

  // ── Purpose-level hyperpriors ─────────────────────────────────────────────
  vector<lower=0, upper=1>[n_pp] mu_alpha_pp;
  vector<lower=0, upper=1>[n_pp] mu_gamma_raw_pp;

  // ── State×purpose hyperpriors ─────────────────────────────────────────────
  vector<lower=0, upper=1>[n_sp] mu_alpha_sp;
  vector<lower=0, upper=1>[n_sp] mu_gamma_raw_sp;

  // ── Zone×purpose hyperpriors ──────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zp] mu_alpha_zp;
  vector<lower=0, upper=1>[n_zp] mu_gamma_raw_zp;

  // ── Purpose-total parameters ──────────────────────────────────────────────
  vector<lower=0, upper=1>[n_pp] alpha_pp;
  vector<lower=0, upper=1>[n_pp] gamma_raw_pp;
  vector<lower=0>[n_pp]          sigma_pp;

  // ── State×purpose parameters ──────────────────────────────────────────────
  vector<lower=0, upper=1>[n_sp] alpha_sp;
  vector<lower=0, upper=1>[n_sp] gamma_raw_sp;
  vector<lower=0>[n_sp]          sigma_sp;

  // ── Zone×purpose parameters ───────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zp] alpha_zp;
  vector<lower=0, upper=1>[n_zp] gamma_raw_zp;
  vector<lower=0>[n_zp]          sigma_zp;

  // ── Region×purpose parameters ─────────────────────────────────────────────
  vector<lower=0, upper=1>[n_rp] alpha_rp;
  vector<lower=0, upper=1>[n_rp] gamma_raw_rp;
  vector<lower=0>[n_rp]          sigma_rp;

  // ── Grand total parameters ────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;

  // ── Initial states ────────────────────────────────────────────────────────
  real                      l0_total;
  vector[m - 1]             s0_total_free;

  vector[n_pp]              l0_pp;
  matrix[n_pp, m - 1]      s0_pp_free;

  vector[n_sp]              l0_sp;
  matrix[n_sp, m - 1]      s0_sp_free;

  vector[n_zp]              l0_zp;
  matrix[n_zp, m - 1]      s0_zp_free;

  vector[n_rp]              l0_rp;
  matrix[n_rp, m - 1]      s0_rp_free;
}
transformed parameters {
  // ── Reconstruct full seasonals with sum-to-zero ───────────────────────────
  vector[m]            s0_total;
  matrix[n_pp, m]      s0_pp;
  matrix[n_sp, m]      s0_sp;
  matrix[n_zp, m]      s0_zp;
  matrix[n_rp, m]      s0_rp;

  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);

  for (n in 1:n_pp) {
    s0_pp[n, 1:(m-1)] = s0_pp_free[n];
    s0_pp[n, m]        = -sum(s0_pp_free[n]);
  }
  for (n in 1:n_sp) {
    s0_sp[n, 1:(m-1)] = s0_sp_free[n];
    s0_sp[n, m]        = -sum(s0_sp_free[n]);
  }
  for (n in 1:n_zp) {
    s0_zp[n, 1:(m-1)] = s0_zp_free[n];
    s0_zp[n, m]        = -sum(s0_zp_free[n]);
  }
  for (n in 1:n_rp) {
    s0_rp[n, 1:(m-1)] = s0_rp_free[n];
    s0_rp[n, m]        = -sum(s0_rp_free[n]);
  }

  // ── Recover actual gamma values ───────────────────────────────────────────
  real              gamma_total;
  vector[n_pp]      gamma_pp;
  vector[n_sp]      gamma_sp;
  vector[n_zp]      gamma_zp;
  vector[n_rp]      gamma_rp;

  gamma_total = gamma_raw_total * (1 - alpha_total);
  for (n in 1:n_pp) gamma_pp[n] = gamma_raw_pp[n] * (1 - alpha_pp[n]);
  for (n in 1:n_sp) gamma_sp[n] = gamma_raw_sp[n] * (1 - alpha_sp[n]);
  for (n in 1:n_zp) gamma_zp[n] = gamma_raw_zp[n] * (1 - alpha_zp[n]);
  for (n in 1:n_rp) gamma_rp[n] = gamma_raw_rp[n] * (1 - alpha_rp[n]);

  // ── Grand total recursion ─────────────────────────────────────────────────
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

  // ── Purpose-total recursions ──────────────────────────────────────────────
  array[n_pp] vector[T_train + 1] l_pp;
  array[n_pp] vector[T_train + m] s_pp;

  for (n in 1:n_pp) {
    for (j in 1:m) s_pp[n][j] = s0_pp[n, j];
    l_pp[n][1] = l0_pp[n];
    for (i in 1:T_train) {
      real s_lag    = s_pp[n][i];
      l_pp[n][i+1] = alpha_pp[n] * (Y_pp[i, n] - s_lag)
                   + (1 - alpha_pp[n]) * l_pp[n][i];
      s_pp[n][i+m] = gamma_pp[n] * (Y_pp[i, n] - l_pp[n][i])
                   + (1 - gamma_pp[n]) * s_lag;
    }
  }

  // ── State×purpose recursions ──────────────────────────────────────────────
  array[n_sp] vector[T_train + 1] l_sp;
  array[n_sp] vector[T_train + m] s_sp;

  for (n in 1:n_sp) {
    for (j in 1:m) s_sp[n][j] = s0_sp[n, j];
    l_sp[n][1] = l0_sp[n];
    for (i in 1:T_train) {
      real s_lag    = s_sp[n][i];
      l_sp[n][i+1] = alpha_sp[n] * (Y_sp[i, n] - s_lag)
                   + (1 - alpha_sp[n]) * l_sp[n][i];
      s_sp[n][i+m] = gamma_sp[n] * (Y_sp[i, n] - l_sp[n][i])
                   + (1 - gamma_sp[n]) * s_lag;
    }
  }

  // ── Zone×purpose recursions ───────────────────────────────────────────────
  array[n_zp] vector[T_train + 1] l_zp;
  array[n_zp] vector[T_train + m] s_zp;

  for (n in 1:n_zp) {
    for (j in 1:m) s_zp[n][j] = s0_zp[n, j];
    l_zp[n][1] = l0_zp[n];
    for (i in 1:T_train) {
      real s_lag    = s_zp[n][i];
      l_zp[n][i+1] = alpha_zp[n] * (Y_zp[i, n] - s_lag)
                   + (1 - alpha_zp[n]) * l_zp[n][i];
      s_zp[n][i+m] = gamma_zp[n] * (Y_zp[i, n] - l_zp[n][i])
                   + (1 - gamma_zp[n]) * s_lag;
    }
  }

  // ── Region×purpose recursions ─────────────────────────────────────────────
  array[n_rp] vector[T_train + 1] l_rp;
  array[n_rp] vector[T_train + m] s_rp;

  for (n in 1:n_rp) {
    for (j in 1:m) s_rp[n][j] = s0_rp[n, j];
    l_rp[n][1] = l0_rp[n];
    for (i in 1:T_train) {
      real s_lag    = s_rp[n][i];
      l_rp[n][i+1] = alpha_rp[n] * (Y_rp[i, n] - s_lag)
                   + (1 - alpha_rp[n]) * l_rp[n][i];
      s_rp[n][i+m] = gamma_rp[n] * (Y_rp[i, n] - l_rp[n][i])
                   + (1 - gamma_rp[n]) * s_lag;
    }
  }
}
model {
  // ── Global hyperpriors ────────────────────────────────────────────────────
  mu_alpha     ~ uniform(0, 1);
  mu_gamma_raw ~ uniform(0, 1);

  // ── Purpose-level hyperpriors ─────────────────────────────────────────────
  mu_alpha_pp     ~ normal(mu_alpha,     0.1);
  mu_gamma_raw_pp ~ normal(mu_gamma_raw, 0.1);

  // ── State×purpose hyperpriors ─────────────────────────────────────────────
  for (n in 1:n_sp) {
    int p = sp_to_pp[n];
    mu_alpha_sp[n]     ~ normal(mu_alpha_pp[p],     0.1);
    mu_gamma_raw_sp[n] ~ normal(mu_gamma_raw_pp[p], 0.1);
  }

  // ── Zone×purpose hyperpriors ──────────────────────────────────────────────
  for (n in 1:n_zp) {
    int s = zp_to_sp[n];
    mu_alpha_zp[n]     ~ normal(mu_alpha_sp[s],     0.1);
    mu_gamma_raw_zp[n] ~ normal(mu_gamma_raw_sp[s], 0.1);
  }

  // ── Grand total parameters ────────────────────────────────────────────────
  alpha_total     ~ normal(mu_alpha,     0.1);
  gamma_raw_total ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 1);

  // ── Purpose-total parameters ──────────────────────────────────────────────
  for (n in 1:n_pp) {
    alpha_pp[n]     ~ normal(mu_alpha_pp[n],     0.1);
    gamma_raw_pp[n] ~ normal(mu_gamma_raw_pp[n], 0.1);
  }
  sigma_pp ~ cauchy(0, 0.5);

  // ── State×purpose parameters ──────────────────────────────────────────────
  for (n in 1:n_sp) {
    int p = sp_to_pp[n];
    alpha_sp[n]     ~ normal(mu_alpha_pp[p],     0.1);
    gamma_raw_sp[n] ~ normal(mu_gamma_raw_pp[p], 0.1);
  }
  sigma_sp ~ cauchy(0, 0.5);

  // ── Zone×purpose parameters ───────────────────────────────────────────────
  for (n in 1:n_zp) {
    int s = zp_to_sp[n];
    alpha_zp[n]     ~ normal(mu_alpha_sp[s],     0.1);
    gamma_raw_zp[n] ~ normal(mu_gamma_raw_sp[s], 0.1);
  }
  sigma_zp ~ cauchy(0, 0.3);

  // ── Region×purpose parameters ─────────────────────────────────────────────
  for (n in 1:n_rp) {
    int z = rp_to_zp[n];
    alpha_rp[n]     ~ normal(mu_alpha_zp[z],     0.1);
    gamma_raw_rp[n] ~ normal(mu_gamma_raw_zp[z], 0.1);
  }
  sigma_rp ~ cauchy(0, 0.3);

  // ── Initial value priors ──────────────────────────────────────────────────
  l0_rp ~ normal(0, 1);
  to_vector(s0_rp_free) ~ normal(0, 1);

  for (n in 1:n_zp) {
    l0_zp[n]      ~ normal(0, rp_per_zp[n]);
    s0_zp_free[n] ~ normal(0, rp_per_zp[n]);
  }
  for (n in 1:n_sp) {
    l0_sp[n]      ~ normal(0, rp_per_sp[n]);
    s0_sp_free[n] ~ normal(0, rp_per_sp[n]);
  }
  for (n in 1:n_pp) {
    l0_pp[n]      ~ normal(0, rp_per_pp[n]);
    s0_pp_free[n] ~ normal(0, rp_per_pp[n]);
  }
  l0_total      ~ normal(0, n_rp);
  s0_total_free ~ normal(0, n_rp);

  // ── Initial-value coherence penalties ─────────────────────────────────────
  for (k in 1:m) {
    real init_total = l0_total       + s0_total[k];
    real init_pp    = sum(l0_pp)     + sum(s0_pp[, k]);
    real init_sp    = sum(l0_sp)     + sum(s0_sp[, k]);
    real init_zp    = sum(l0_zp)     + sum(s0_zp[, k]);
    real init_rp    = sum(l0_rp)     + sum(s0_rp[, k]);

    target += normal_lpdf(init_total - init_pp | 0, sigma_coh);
    target += normal_lpdf(init_pp    - init_sp | 0, sigma_coh);
    target += normal_lpdf(init_sp    - init_zp | 0, sigma_coh);
    target += normal_lpdf(init_zp    - init_rp | 0, sigma_coh);
  }

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (n in 1:n_rp)
    for (i in 1:T_train)
      target += normal_lpdf(Y_rp[i, n] | l_rp[n][i] + s_rp[n][i], sigma_rp[n]);

  for (n in 1:n_zp)
    for (i in 1:T_train)
      target += lambda_zp * normal_lpdf(Y_zp[i, n] | l_zp[n][i] + s_zp[n][i], sigma_zp[n]);

  for (n in 1:n_sp)
    for (i in 1:T_train)
      target += lambda_sp * normal_lpdf(Y_sp[i, n] | l_sp[n][i] + s_sp[n][i], sigma_sp[n]);

  for (n in 1:n_pp)
    for (i in 1:T_train)
      target += lambda_pp * normal_lpdf(Y_pp[i, n] | l_pp[n][i] + s_pp[n][i], sigma_pp[n]);

  for (i in 1:T_train)
    target += lambda_total * normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);

  // ── Fitted-value coherence penalties ──────────────────────────────────────
  for (i in 1:T_train) {
    real fitted_total = l_t[i] + s_t[i];
    real fitted_pp    = 0;
    real fitted_sp    = 0;
    real fitted_zp    = 0;
    real fitted_rp    = 0;

    for (n in 1:n_pp) fitted_pp += l_pp[n][i] + s_pp[n][i];
    for (n in 1:n_sp) fitted_sp += l_sp[n][i] + s_sp[n][i];
    for (n in 1:n_zp) fitted_zp += l_zp[n][i] + s_zp[n][i];
    for (n in 1:n_rp) fitted_rp += l_rp[n][i] + s_rp[n][i];

    target += normal_lpdf(fitted_total - fitted_pp | 0, sigma_coh);
    target += normal_lpdf(fitted_pp    - fitted_sp | 0, sigma_coh);
    target += normal_lpdf(fitted_sp    - fitted_zp | 0, sigma_coh);
    target += normal_lpdf(fitted_zp    - fitted_rp | 0, sigma_coh);
  }
}
generated quantities {
  matrix[H, n_rp] Y_fc_rp;
  matrix[H, n_zp] Y_fc_zp;
  matrix[H, n_sp] Y_fc_sp;
  matrix[H, n_pp] Y_fc_pp;
  vector[H]       Y_fc_total;

  for (n in 1:n_rp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_rp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_rp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_rp[h, n] = normal_rng(lf[h] + sf[h], sigma_rp[n]);
      lf[h+1] = alpha_rp[n] * (Y_fc_rp[h, n] - sf[h])
              + (1 - alpha_rp[n]) * lf[h];
    }
    Y_fc_rp[H, n] = normal_rng(lf[H] + sf[H], sigma_rp[n]);
  }

  for (n in 1:n_zp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_zp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_zp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_zp[h, n] = normal_rng(lf[h] + sf[h], sigma_zp[n]);
      lf[h+1] = alpha_zp[n] * (Y_fc_zp[h, n] - sf[h])
              + (1 - alpha_zp[n]) * lf[h];
    }
    Y_fc_zp[H, n] = normal_rng(lf[H] + sf[H], sigma_zp[n]);
  }

  for (n in 1:n_sp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_sp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_sp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_sp[h, n] = normal_rng(lf[h] + sf[h], sigma_sp[n]);
      lf[h+1] = alpha_sp[n] * (Y_fc_sp[h, n] - sf[h])
              + (1 - alpha_sp[n]) * lf[h];
    }
    Y_fc_sp[H, n] = normal_rng(lf[H] + sf[H], sigma_sp[n]);
  }

  for (n in 1:n_pp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_pp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_pp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_pp[h, n] = normal_rng(lf[h] + sf[h], sigma_pp[n]);
      lf[h+1] = alpha_pp[n] * (Y_fc_pp[h, n] - sf[h])
              + (1 - alpha_pp[n]) * lf[h];
    }
    Y_fc_pp[H, n] = normal_rng(lf[H] + sf[H], sigma_pp[n]);
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

# ── Compile and sample ─────────────────────────────────────────────────────────

stanDso <- stan_model(model_code = model_string)

stanFit <- sampling(
  object = stanDso,
  data   = list(
    T_train    = window_size,
    m          = m,
    n_rp       = n_rp,
    n_zp       = n_zp,
    n_sp       = n_sp,
    n_pp       = n_pp,
    H          = horizon + 1,
    Y_rp       = Y_rp[train_idx, ],
    Y_zp       = Y_zp[train_idx, ],
    Y_sp       = Y_sp[train_idx, ],
    Y_pp       = Y_pp[train_idx, ],
    Y_total    = Y_total[train_idx],
    rp_to_zp   = rp_to_zp,
    rp_to_sp   = rp_to_sp,
    rp_to_pp   = rp_to_pp,
    zp_to_sp   = zp_to_sp,
    zp_to_pp   = zp_to_pp,
    sp_to_pp   = sp_to_pp,
    rp_per_zp  = rp_per_zp,
    rp_per_sp  = rp_per_sp,
    rp_per_pp  = rp_per_pp,
    lambda_total = lambda_total,
    lambda_pp    = lambda_pp,
    lambda_sp    = lambda_sp,
    lambda_zp    = lambda_zp, 
    sigma_coh    = sigma_coh
  ),
  seed    = 2025,
  chains  = 10,
  iter    = 2000,
  warmup  = 1000,
  cores   = 10,
  refresh = 0
)

saveRDS(stanFit, file.path(dir_stanfits, sprintf("stanfit_timestep_%03d.rds", origin)))

# ── Helpers ────────────────────────────────────────────────────────────────────

compute_energy_score <- function(draws_mat, actuals, sub_idx = NULL) {
  if (!is.null(sub_idx)) draws_mat <- draws_mat[sub_idx, , drop = FALSE]
  col_sds    <- apply(draws_mat, 2, sd)
  draws_sc   <- sweep(draws_mat, 2, col_sds, "/")
  actuals_sc <- as.numeric(actuals) / col_sds
  diffs1     <- sweep(draws_sc, 2, actuals_sc, "-")
  t1         <- mean(sqrt(rowSums(diffs1^2)))
  pw         <- dist(draws_sc, method = "euclidean")
  t2         <- mean(pw)
  t1 - 0.5 * t2
}

JSD_metric <- function(p_draws, q_draws, n_grid = 1000) {
  lo    <- min(c(p_draws, q_draws)) - 3 * sd(c(p_draws, q_draws))
  hi    <- max(c(p_draws, q_draws)) + 3 * sd(c(p_draws, q_draws))
  p_kde <- pmax(density(p_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  q_kde <- pmax(density(q_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  m_kde <- 0.5 * (p_kde + q_kde)
  dx    <- (hi - lo) / (n_grid - 1)
  jsd   <- 0.5 * (sum(p_kde * log(p_kde / m_kde)) +
                    sum(q_kde * log(q_kde / m_kde))) * dx
  sqrt(max(0, min(log(2), jsd)))
}

# ── Extract draws ──────────────────────────────────────────────────────────────

fc_rp_draws    <- rstan::extract(stanFit, pars = "Y_fc_rp")$Y_fc_rp
fc_zp_draws    <- rstan::extract(stanFit, pars = "Y_fc_zp")$Y_fc_zp
fc_sp_draws    <- rstan::extract(stanFit, pars = "Y_fc_sp")$Y_fc_sp
fc_pp_draws    <- rstan::extract(stanFit, pars = "Y_fc_pp")$Y_fc_pp
fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total

n_draws <- nrow(fc_rp_draws)

# ── Back-transform to original scale ──────────────────────────────────────────

fc_rp_orig <- fc_rp_draws[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean

fc_zp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_zp))
for (z in 1:n_zp)
  fc_zp_orig[, , z] <- fc_zp_draws[, 1:horizon, z] * Y_sd + Y_mean * rp_per_zp[z]

fc_sp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_sp))
for (s in 1:n_sp)
  fc_sp_orig[, , s] <- fc_sp_draws[, 1:horizon, s] * Y_sd + Y_mean * rp_per_sp[s]

fc_pp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_pp))
for (p in 1:n_pp)
  fc_pp_orig[, , p] <- fc_pp_draws[, 1:horizon, p] * Y_sd + Y_mean * rp_per_pp[p]

fc_total_orig <- fc_total_draws[, 1:horizon, drop = FALSE] * Y_sd + Y_mean * n_rp

# ── Actuals ────────────────────────────────────────────────────────────────────

actuals_rp    <- Y_raw_rp[test_idx, ]
actuals_zp    <- Y_raw_zp[test_idx, ]
actuals_sp    <- Y_raw_sp[test_idx, ]
actuals_pp    <- Y_raw_pp[test_idx, ]
actuals_total <- Y_raw_total[test_idx]

# ── Subsample index for ES ─────────────────────────────────────────────────────

n_es_subsample <- 1000L
set.seed(origin)
sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)

# ── Per-horizon loop ───────────────────────────────────────────────────────────

horizon_results <- lapply(1:horizon, function(h) {
  
  make_rows <- function(draws_arr, actuals_mat, codes, level) {
    lapply(seq_len(ncol(actuals_mat)), function(i) {
      draws  <- draws_arr[, h, i]
      actual <- actuals_mat[h, i]
      data.frame(
        origin     = origin,
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
    origin     = origin,
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
    make_rows(fc_rp_orig, actuals_rp, rp_codes, "RegionPurpose"),
    make_rows(fc_zp_orig, actuals_zp, zp_codes, "ZonePurpose"),
    make_rows(fc_sp_orig, actuals_sp, sp_codes, "StatePurpose"),
    make_rows(fc_pp_orig, actuals_pp, pp_codes, "Purpose"),
    list(total_row)
  ))
  
  rp_mat <- fc_rp_orig[, h, ]
  zp_mat <- fc_zp_orig[, h, ]
  sp_mat <- fc_sp_orig[, h, ]
  pp_mat <- fc_pp_orig[, h, ]
  
  sum_rp_draws <- rowSums(rp_mat)
  sum_zp_draws <- rowSums(zp_mat)
  sum_sp_draws <- rowSums(sp_mat)
  sum_pp_draws <- rowSums(pp_mat)
  
  joint_mat     <- cbind(rp_mat, zp_mat, sp_mat, pp_mat, draws_t)
  actuals_joint <- c(actuals_rp[h, ], actuals_zp[h, ],
                     actuals_sp[h, ], actuals_pp[h, ], actual_t)
  es_val <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
  
  multivariate_df <- data.frame(
    origin           = origin,
    horizon          = h,
    es               = es_val,
    mai_rp_zp        = mean(sum_rp_draws) - mean(sum_zp_draws),
    mai_zp_sp        = mean(sum_zp_draws) - mean(sum_sp_draws),
    mai_sp_pp        = mean(sum_sp_draws) - mean(sum_pp_draws),
    mai_pp_total     = mean(sum_pp_draws) - mean(draws_t),
    mai_rp_total     = mean(sum_rp_draws) - mean(draws_t),
    jsd_rp_zp        = JSD_metric(sum_rp_draws, sum_zp_draws),
    jsd_zp_sp        = JSD_metric(sum_zp_draws, sum_sp_draws),
    jsd_sp_pp        = JSD_metric(sum_sp_draws, sum_pp_draws),
    jsd_pp_total     = JSD_metric(sum_pp_draws, draws_t),
    jsd_rp_total     = JSD_metric(sum_rp_draws, draws_t)
  )
  
  list(univariate = univariate_df, multivariate = multivariate_df)
})

# ── Collect and save ───────────────────────────────────────────────────────────

origin_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
origin_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))

saveRDS(
  list(
    univariate   = origin_univariate,
    multivariate = origin_multivariate
  ),
  file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin))
)