library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)
library(Matrix)

load("data/tourism.Rdata")

# ── Read job name and directories ──────────────────────────────────────────────

args     <- commandArgs(trailingOnly = TRUE)
job_name <- args[1]

dir_accuracy      <- file.path(job_name, "accuracy")
dir_stanfits      <- file.path(job_name, "stanfits")
dir_mint_samples  <- file.path(job_name, "MinT")
dir_mint_accuracy <- file.path(job_name, "accuracy_MinT")

for (d in c(dir_accuracy, dir_stanfits, dir_mint_samples, dir_mint_accuracy))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Series identification ──────────────────────────────────────────────────────

all_names <- names(tourism)

rp_codes <- all_names[grepl("^[A-G][A-Z][A-Z](Hol|Vis|Bus|Oth)$", all_names)]
zp_codes <- all_names[grepl("^[A-G][A-Z](Hol|Vis|Bus|Oth)$",      all_names)]
sp_codes <- all_names[grepl("^[A-G](Hol|Vis|Bus|Oth)$",           all_names)]
pp_codes <- all_names[grepl("^Total(Hol|Vis|Bus|Oth)$",           all_names)]

n_rp  <- length(rp_codes)
n_zp  <- length(zp_codes)
n_sp  <- length(sp_codes)
n_pp  <- length(pp_codes)
n_all <- n_rp + n_zp + n_sp + n_pp + 1L

# ── Hierarchy mappings (back-transform + S matrix only) ───────────────────────

rp_to_zp <- match(paste0(substr(rp_codes, 1, 2),
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), zp_codes)
rp_to_sp <- match(paste0(substr(rp_codes, 1, 1),
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
rp_to_pp <- match(paste0("Total",
                         str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)

rp_per_zp <- as.integer(table(rp_to_zp)[order(as.integer(names(table(rp_to_zp))))])
rp_per_sp <- as.integer(table(rp_to_sp)[order(as.integer(names(table(rp_to_sp))))])
rp_per_pp <- as.integer(table(rp_to_pp)[order(as.integer(names(table(rp_to_pp))))])

# ── Summing matrix S  [n_all x n_rp] ──────────────────────────────────────────

build_S <- function() {
  S <- matrix(0L, nrow = n_all, ncol = n_rp)
  S[1, ] <- 1L
  for (p in seq_len(n_pp))  S[1L + p,                    rp_to_pp == p] <- 1L
  for (s in seq_len(n_sp))  S[1L + n_pp + s,             rp_to_sp == s] <- 1L
  for (z in seq_len(n_zp))  S[1L + n_pp + n_sp + z,      rp_to_zp == z] <- 1L
  for (r in seq_len(n_rp))  S[1L + n_pp + n_sp + n_zp + r, r]           <- 1L
  S
}

S <- build_S()

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

task_id <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
origin  <- origins[task_id]
# origin <- 186

train_idx <- (origin - window_size + 1):origin
test_idx  <- (origin + 1):(origin + horizon)

# ── Scaling ────────────────────────────────────────────────────────────────────

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

  array[n_zp] int<lower=1> rp_per_zp;
  array[n_sp] int<lower=1> rp_per_sp;
  array[n_pp] int<lower=1> rp_per_pp;
  int<lower=1> n_rp_total;
}
parameters {
  // ── Region×purpose ────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_rp] alpha_rp;
  vector<lower=0, upper=1>[n_rp] gamma_raw_rp;
  vector<lower=0>[n_rp]          sigma_rp;
  vector[n_rp]                   l0_rp;
  matrix[n_rp, m - 1]            s0_rp_free;

  // ── Zone×purpose ──────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zp] alpha_zp;
  vector<lower=0, upper=1>[n_zp] gamma_raw_zp;
  vector<lower=0>[n_zp]          sigma_zp;
  vector[n_zp]                   l0_zp;
  matrix[n_zp, m - 1]            s0_zp_free;

  // ── State×purpose ─────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_sp] alpha_sp;
  vector<lower=0, upper=1>[n_sp] gamma_raw_sp;
  vector<lower=0>[n_sp]          sigma_sp;
  vector[n_sp]                   l0_sp;
  matrix[n_sp, m - 1]            s0_sp_free;

  // ── Purpose totals ────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_pp] alpha_pp;
  vector<lower=0, upper=1>[n_pp] gamma_raw_pp;
  vector<lower=0>[n_pp]          sigma_pp;
  vector[n_pp]                   l0_pp;
  matrix[n_pp, m - 1]            s0_pp_free;

  // ── Grand total ───────────────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;
  real                   l0_total;
  vector[m - 1]          s0_total_free;
}
transformed parameters {
  matrix[n_rp, m] s0_rp;
  matrix[n_zp, m] s0_zp;
  matrix[n_sp, m] s0_sp;
  matrix[n_pp, m] s0_pp;
  vector[m]       s0_total;

  for (n in 1:n_rp) { s0_rp[n, 1:(m-1)] = s0_rp_free[n]; s0_rp[n, m] = -sum(s0_rp_free[n]); }
  for (n in 1:n_zp) { s0_zp[n, 1:(m-1)] = s0_zp_free[n]; s0_zp[n, m] = -sum(s0_zp_free[n]); }
  for (n in 1:n_sp) { s0_sp[n, 1:(m-1)] = s0_sp_free[n]; s0_sp[n, m] = -sum(s0_sp_free[n]); }
  for (n in 1:n_pp) { s0_pp[n, 1:(m-1)] = s0_pp_free[n]; s0_pp[n, m] = -sum(s0_pp_free[n]); }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);

  vector[n_rp] gamma_rp;
  vector[n_zp] gamma_zp;
  vector[n_sp] gamma_sp;
  vector[n_pp] gamma_pp;
  real         gamma_total;

  for (n in 1:n_rp) gamma_rp[n] = gamma_raw_rp[n] * (1 - alpha_rp[n]);
  for (n in 1:n_zp) gamma_zp[n] = gamma_raw_zp[n] * (1 - alpha_zp[n]);
  for (n in 1:n_sp) gamma_sp[n] = gamma_raw_sp[n] * (1 - alpha_sp[n]);
  for (n in 1:n_pp) gamma_pp[n] = gamma_raw_pp[n] * (1 - alpha_pp[n]);
  gamma_total = gamma_raw_total * (1 - alpha_total);

  array[n_rp] vector[T_train + 1] l_rp;
  array[n_rp] vector[T_train + m] s_rp;
  for (n in 1:n_rp) {
    for (j in 1:m) s_rp[n][j] = s0_rp[n, j];
    l_rp[n][1] = l0_rp[n];
    for (i in 1:T_train) {
      real s_lag    = s_rp[n][i];
      l_rp[n][i+1] = alpha_rp[n] * (Y_rp[i, n] - s_lag) + (1 - alpha_rp[n]) * l_rp[n][i];
      s_rp[n][i+m] = gamma_rp[n] * (Y_rp[i, n] - l_rp[n][i]) + (1 - gamma_rp[n]) * s_lag;
    }
  }

  array[n_zp] vector[T_train + 1] l_zp;
  array[n_zp] vector[T_train + m] s_zp;
  for (n in 1:n_zp) {
    for (j in 1:m) s_zp[n][j] = s0_zp[n, j];
    l_zp[n][1] = l0_zp[n];
    for (i in 1:T_train) {
      real s_lag    = s_zp[n][i];
      l_zp[n][i+1] = alpha_zp[n] * (Y_zp[i, n] - s_lag) + (1 - alpha_zp[n]) * l_zp[n][i];
      s_zp[n][i+m] = gamma_zp[n] * (Y_zp[i, n] - l_zp[n][i]) + (1 - gamma_zp[n]) * s_lag;
    }
  }

  array[n_sp] vector[T_train + 1] l_sp;
  array[n_sp] vector[T_train + m] s_sp;
  for (n in 1:n_sp) {
    for (j in 1:m) s_sp[n][j] = s0_sp[n, j];
    l_sp[n][1] = l0_sp[n];
    for (i in 1:T_train) {
      real s_lag    = s_sp[n][i];
      l_sp[n][i+1] = alpha_sp[n] * (Y_sp[i, n] - s_lag) + (1 - alpha_sp[n]) * l_sp[n][i];
      s_sp[n][i+m] = gamma_sp[n] * (Y_sp[i, n] - l_sp[n][i]) + (1 - gamma_sp[n]) * s_lag;
    }
  }

  array[n_pp] vector[T_train + 1] l_pp;
  array[n_pp] vector[T_train + m] s_pp;
  for (n in 1:n_pp) {
    for (j in 1:m) s_pp[n][j] = s0_pp[n, j];
    l_pp[n][1] = l0_pp[n];
    for (i in 1:T_train) {
      real s_lag    = s_pp[n][i];
      l_pp[n][i+1] = alpha_pp[n] * (Y_pp[i, n] - s_lag) + (1 - alpha_pp[n]) * l_pp[n][i];
      s_pp[n][i+m] = gamma_pp[n] * (Y_pp[i, n] - l_pp[n][i]) + (1 - gamma_pp[n]) * s_lag;
    }
  }

  vector[T_train + 1] l_t;
  vector[T_train + m] s_t;
  for (j in 1:m) s_t[j] = s0_total[j];
  l_t[1] = l0_total;
  for (i in 1:T_train) {
    real s_lag = s_t[i];
    l_t[i+1]  = alpha_total * (Y_total[i] - s_lag) + (1 - alpha_total) * l_t[i];
    s_t[i+m]  = gamma_total * (Y_total[i] - l_t[i]) + (1 - gamma_total) * s_lag;
  }
}
model {
  // ── Independent priors — no shared hyperparameters ────────────────────────

  alpha_rp        ~ uniform(0, 1);
  gamma_raw_rp    ~ uniform(0, 1);
  alpha_zp        ~ uniform(0, 1);
  gamma_raw_zp    ~ uniform(0, 1);
  alpha_sp        ~ uniform(0, 1);
  gamma_raw_sp    ~ uniform(0, 1);
  alpha_pp        ~ uniform(0, 1);
  gamma_raw_pp    ~ uniform(0, 1);
  alpha_total     ~ uniform(0, 1);
  gamma_raw_total ~ uniform(0, 1);

  sigma_rp    ~ cauchy(0, 0.3);
  sigma_zp    ~ cauchy(0, 0.3);
  sigma_sp    ~ cauchy(0, 0.5);
  sigma_pp    ~ cauchy(0, 0.5);
  sigma_total ~ cauchy(0, 1);

  l0_rp ~ normal(0, 1);
  to_vector(s0_rp_free) ~ normal(0, 1);
  for (n in 1:n_zp) { l0_zp[n] ~ normal(0, rp_per_zp[n]); s0_zp_free[n] ~ normal(0, rp_per_zp[n]); }
  for (n in 1:n_sp) { l0_sp[n] ~ normal(0, rp_per_sp[n]); s0_sp_free[n] ~ normal(0, rp_per_sp[n]); }
  for (n in 1:n_pp) { l0_pp[n] ~ normal(0, rp_per_pp[n]); s0_pp_free[n] ~ normal(0, rp_per_pp[n]); }
  l0_total      ~ normal(0, n_rp_total);
  s0_total_free ~ normal(0, n_rp_total);

  // ── Likelihoods (equal weight, no lambda) ────────────────────────────────
  for (n in 1:n_rp)
    for (i in 1:T_train)
      target += normal_lpdf(Y_rp[i, n] | l_rp[n][i] + s_rp[n][i], sigma_rp[n]);
  for (n in 1:n_zp)
    for (i in 1:T_train)
      target += normal_lpdf(Y_zp[i, n] | l_zp[n][i] + s_zp[n][i], sigma_zp[n]);
  for (n in 1:n_sp)
    for (i in 1:T_train)
      target += normal_lpdf(Y_sp[i, n] | l_sp[n][i] + s_sp[n][i], sigma_sp[n]);
  for (n in 1:n_pp)
    for (i in 1:T_train)
      target += normal_lpdf(Y_pp[i, n] | l_pp[n][i] + s_pp[n][i], sigma_pp[n]);
  for (i in 1:T_train)
    target += normal_lpdf(Y_total[i] | l_t[i] + s_t[i], sigma_total);
}
generated quantities {
  matrix[T_train, n_rp] fitted_rp;
  matrix[T_train, n_zp] fitted_zp;
  matrix[T_train, n_sp] fitted_sp;
  matrix[T_train, n_pp] fitted_pp;
  vector[T_train]        fitted_total;

  for (n in 1:n_rp) for (i in 1:T_train) fitted_rp[i, n]  = l_rp[n][i] + s_rp[n][i];
  for (n in 1:n_zp) for (i in 1:T_train) fitted_zp[i, n]  = l_zp[n][i] + s_zp[n][i];
  for (n in 1:n_sp) for (i in 1:T_train) fitted_sp[i, n]  = l_sp[n][i] + s_sp[n][i];
  for (n in 1:n_pp) for (i in 1:T_train) fitted_pp[i, n]  = l_pp[n][i] + s_pp[n][i];
  for (i in 1:T_train)                   fitted_total[i]   = l_t[i]     + s_t[i];

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
      lf[h+1] = alpha_rp[n] * (Y_fc_rp[h, n] - sf[h]) + (1 - alpha_rp[n]) * lf[h];
    }
    Y_fc_rp[H, n] = normal_rng(lf[H] + sf[H], sigma_rp[n]);
  }
  for (n in 1:n_zp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_zp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_zp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_zp[h, n] = normal_rng(lf[h] + sf[h], sigma_zp[n]);
      lf[h+1] = alpha_zp[n] * (Y_fc_zp[h, n] - sf[h]) + (1 - alpha_zp[n]) * lf[h];
    }
    Y_fc_zp[H, n] = normal_rng(lf[H] + sf[H], sigma_zp[n]);
  }
  for (n in 1:n_sp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_sp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_sp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_sp[h, n] = normal_rng(lf[h] + sf[h], sigma_sp[n]);
      lf[h+1] = alpha_sp[n] * (Y_fc_sp[h, n] - sf[h]) + (1 - alpha_sp[n]) * lf[h];
    }
    Y_fc_sp[H, n] = normal_rng(lf[H] + sf[H], sigma_sp[n]);
  }
  for (n in 1:n_pp) {
    vector[H] lf; vector[H] sf;
    lf[1] = l_pp[n][T_train + 1];
    for (h in 1:H) sf[h] = s_pp[n][T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_pp[h, n] = normal_rng(lf[h] + sf[h], sigma_pp[n]);
      lf[h+1] = alpha_pp[n] * (Y_fc_pp[h, n] - sf[h]) + (1 - alpha_pp[n]) * lf[h];
    }
    Y_fc_pp[H, n] = normal_rng(lf[H] + sf[H], sigma_pp[n]);
  }
  {
    vector[H] lf; vector[H] sf;
    lf[1] = l_t[T_train + 1];
    for (h in 1:H) sf[h] = s_t[T_train + h - m];
    for (h in 1:(H-1)) {
      Y_fc_total[h] = normal_rng(lf[h] + sf[h], sigma_total);
      lf[h+1] = alpha_total * (Y_fc_total[h] - sf[h]) + (1 - alpha_total) * lf[h];
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
    rp_per_zp  = rp_per_zp,
    rp_per_sp  = rp_per_sp,
    rp_per_pp  = rp_per_pp,
    n_rp_total = n_rp
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
compute_energy_score <- function(draws_mat, actuals, n_subsample = NULL) {
  n_draws <- nrow(draws_mat)
  if (!is.null(n_subsample) && n_subsample < n_draws) {
    idx1 <- sample.int(n_draws, size = n_subsample, replace = FALSE)
    idx2 <- sample.int(n_draws, size = n_subsample, replace = FALSE)
    x    <- draws_mat[idx1, , drop = FALSE]
    xs   <- draws_mat[idx2, , drop = FALSE]
  } else {
    x    <- draws_mat
    xs   <- draws_mat[sample.int(n_draws), , drop = FALSE] # Crucial Fix: shuffle rows
  }
  diffs <- sweep(x, 2, as.numeric(actuals), "-")
  t1    <- mean(sqrt(rowSums(diffs^2)))
  t2    <- mean(sqrt(rowSums((x - xs)^2)))
  return(t1 - 0.5 * t2)
}

JSD_metric <- function(p_draws, q_draws, n_grid = 1000) {
  lo    <- min(c(p_draws, q_draws)) - 3 * sd(c(p_draws, q_draws))
  hi    <- max(c(p_draws, q_draws)) + 3 * sd(c(p_draws, q_draws))
  p_kde <- pmax(density(p_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  q_kde <- pmax(density(q_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  m_kde <- 0.5 * (p_kde + q_kde)
  dx    <- (hi - lo) / (n_grid - 1)
  jsd   <- 0.5 * (sum(p_kde * log(p_kde / m_kde)) + sum(q_kde * log(q_kde / m_kde))) * dx
  sqrt(max(0, min(log(2), jsd)))
}

# ── Single-pass extraction ─────────────────────────────────────────────────────

all_pars <- rstan::extract(stanFit, pars = c(
  "Y_fc_rp", "Y_fc_zp", "Y_fc_sp", "Y_fc_pp", "Y_fc_total",
  "fitted_rp", "fitted_zp", "fitted_sp", "fitted_pp", "fitted_total"
))

fc_rp_draws    <- all_pars$Y_fc_rp
fc_zp_draws    <- all_pars$Y_fc_zp
fc_sp_draws    <- all_pars$Y_fc_sp
fc_pp_draws    <- all_pars$Y_fc_pp
fc_total_draws <- all_pars$Y_fc_total

n_draws <- nrow(fc_rp_draws)

# ── Back-transform to original scale ──────────────────────────────────────────

fc_rp_orig    <- fc_rp_draws[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean
fc_total_orig <- fc_total_draws[, 1:horizon, drop = FALSE] * Y_sd + Y_mean * n_rp

fc_zp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_zp))
for (z in seq_len(n_zp))
  fc_zp_orig[, , z] <- fc_zp_draws[, 1:horizon, z] * Y_sd + Y_mean * rp_per_zp[z]

fc_sp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_sp))
for (s in seq_len(n_sp))
  fc_sp_orig[, , s] <- fc_sp_draws[, 1:horizon, s] * Y_sd + Y_mean * rp_per_sp[s]

fc_pp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_pp))
for (p in seq_len(n_pp))
  fc_pp_orig[, , p] <- fc_pp_draws[, 1:horizon, p] * Y_sd + Y_mean * rp_per_pp[p]

rm(fc_rp_draws, fc_zp_draws, fc_sp_draws, fc_pp_draws, fc_total_draws)

# ── Build posterior-mean residuals for MinT W matrix ──────────────────────────

pm_fitted_rp    <- apply(all_pars$fitted_rp,    2:3, mean)
pm_fitted_zp    <- apply(all_pars$fitted_zp,    2:3, mean)
pm_fitted_sp    <- apply(all_pars$fitted_sp,    2:3, mean)
pm_fitted_pp    <- apply(all_pars$fitted_pp,    2:3, mean)
pm_fitted_total <- apply(all_pars$fitted_total, 2,   mean)

rm(all_pars)

fit_rp_orig    <- pm_fitted_rp    * Y_sd + Y_mean
fit_total_orig <- pm_fitted_total * Y_sd + Y_mean * n_rp

fit_zp_orig <- pm_fitted_zp
for (z in seq_len(n_zp))
  fit_zp_orig[, z] <- pm_fitted_zp[, z] * Y_sd + Y_mean * rp_per_zp[z]

fit_sp_orig <- pm_fitted_sp
for (s in seq_len(n_sp))
  fit_sp_orig[, s] <- pm_fitted_sp[, s] * Y_sd + Y_mean * rp_per_sp[s]

fit_pp_orig <- pm_fitted_pp
for (p in seq_len(n_pp))
  fit_pp_orig[, p] <- pm_fitted_pp[, p] * Y_sd + Y_mean * rp_per_pp[p]

E <- cbind(
  Y_raw_total[train_idx] - fit_total_orig,
  Y_raw_pp[train_idx, ]  - fit_pp_orig,
  Y_raw_sp[train_idx, ]  - fit_sp_orig,
  Y_raw_zp[train_idx, ]  - fit_zp_orig,
  Y_raw_rp[train_idx, ]  - fit_rp_orig
)

rm(pm_fitted_rp, pm_fitted_zp, pm_fitted_sp, pm_fitted_pp, pm_fitted_total,
   fit_rp_orig, fit_zp_orig, fit_sp_orig, fit_pp_orig, fit_total_orig)

# ── Ledoit-Wolf shrinkage covariance (MinT-shr) ────────────────────────────────

lw_shrink <- function(E) {
  T   <- nrow(E)
  Sc  <- cov(E)
  d   <- diag(Sc)
  E_c <- sweep(E, 2, colMeans(E), "-")
  xs      <- E_c^2
  xss     <- crossprod(E_c)
  var_sij <- (T / (T - 1)^3) * (crossprod(xs) - (1 / T) * xss^2)
  denom   <- sum(xss[upper.tri(xss)]^2) / (T - 1)^2
  rho_hat <- if (denom < .Machine$double.eps) 0 else
    min(1, sum(var_sij[upper.tri(var_sij)]) / denom)
  (1 - rho_hat) * Sc + rho_hat * diag(d)
}

W <- lw_shrink(E)
rm(E)

# ── Precompute MinT projection matrix P = S(S'W^{-1}S)^{-1}S'W^{-1} ──────────

S_mat <- S
storage.mode(S_mat) <- "double"

WinvS <- tryCatch(
  solve(W, S_mat),
  error = function(e)
    solve(W + diag(1e-6 * mean(diag(W)), nrow(W)), S_mat)
)

StWinvS_chol <- chol(crossprod(S_mat, WinvS))

# P_mint: [n_all x n_all]
P_mint <- S_mat %*% chol2inv(StWinvS_chol) %*% t(WinvS)

rm(W, WinvS, StWinvS_chol)

# ── Define Index Positions Within Matrix Structure ───────────────────────────
# Consistent with: total | pp (n_pp) | sp (n_sp) | zp (n_zp) | rp (n_rp)

idx_total <- 1L
idx_pp    <- 1L + seq_len(n_pp)
idx_sp    <- 1L + n_pp + seq_len(n_sp)
idx_zp    <- 1L + n_pp + n_sp + seq_len(n_zp)
idx_rp    <- 1L + n_pp + n_sp + n_zp + seq_len(n_rp)

# ── Apply MinT to all draws for each horizon ───────────────────────────────────

fc_rp_mint    <- array(NA_real_, dim = c(n_draws, horizon, n_rp))
fc_zp_mint    <- array(NA_real_, dim = c(n_draws, horizon, n_zp))
fc_sp_mint    <- array(NA_real_, dim = c(n_draws, horizon, n_sp))
fc_pp_mint    <- array(NA_real_, dim = c(n_draws, horizon, n_pp))
fc_total_mint <- matrix(NA_real_, nrow = n_draws, ncol = horizon)

for (h in seq_len(horizon)) {
  # Base draw matrix for this horizon: [n_draws x n_all]
  base_h <- cbind(
    fc_total_orig[, h],
    fc_pp_orig[, h, ],
    fc_sp_orig[, h, ],
    fc_zp_orig[, h, ],
    fc_rp_orig[, h, ]
  )
  
  # Vectorised reconciliation
  rec_h <- base_h %*% t(P_mint)   # [n_draws x n_all]
  
  fc_total_mint[, h]  <- rec_h[, idx_total]
  fc_pp_mint[, h, ]   <- rec_h[, idx_pp,   drop = FALSE]
  fc_sp_mint[, h, ]   <- rec_h[, idx_sp,   drop = FALSE]
  fc_zp_mint[, h, ]   <- rec_h[, idx_zp,   drop = FALSE]
  fc_rp_mint[, h, ]   <- rec_h[, idx_rp,   drop = FALSE]
}

rm(P_mint)

saveRDS(
  list(
    fc_rp_mint    = fc_rp_mint,
    fc_zp_mint    = fc_zp_mint,
    fc_sp_mint    = fc_sp_mint,
    fc_pp_mint    = fc_pp_mint,
    fc_total_mint = fc_total_mint
  ),
  file.path(dir_mint_samples, sprintf("mint_timestep_%03d.rds", origin))
)

# ── Actuals ────────────────────────────────────────────────────────────────────

actuals_rp    <- Y_raw_rp[test_idx, ]
actuals_zp    <- Y_raw_zp[test_idx, ]
actuals_sp    <- Y_raw_sp[test_idx, ]
actuals_pp    <- Y_raw_pp[test_idx, ]
actuals_total <- Y_raw_total[test_idx]

# ── Accuracy helper (shared by base and MinT) ──────────────────────────────────

n_es_subsample <- 1000L
set.seed(origin)

make_accuracy <- function(fc_rp_arr, fc_zp_arr, fc_sp_arr, fc_pp_arr,
                          fc_total_mat, label_origin) {
  
  make_rows <- function(draws_arr, actuals_mat, codes, level) {
    lapply(seq_len(ncol(actuals_mat)), function(i) {
      lapply(seq_len(horizon), function(h) {
        draws  <- draws_arr[, h, i]
        actual <- actuals_mat[h, i]
        data.frame(
          origin     = label_origin,
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
    })
  }
  
  total_rows <- lapply(seq_len(horizon), function(h) {
    draws  <- fc_total_mat[, h]
    actual <- actuals_total[h]
    data.frame(
      origin     = label_origin,
      horizon    = h,
      series     = "Total",
      level      = "Total",
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
  
  univariate_df <- bind_rows(c(
    unlist(make_rows(fc_rp_arr, actuals_rp, rp_codes, "RegionPurpose"), recursive = FALSE),
    unlist(make_rows(fc_zp_arr, actuals_zp, zp_codes, "ZonePurpose"),   recursive = FALSE),
    unlist(make_rows(fc_sp_arr, actuals_sp, sp_codes, "StatePurpose"),  recursive = FALSE),
    unlist(make_rows(fc_pp_arr, actuals_pp, pp_codes, "Purpose"),       recursive = FALSE),
    total_rows
  ))
  
  multivariate_df <- bind_rows(lapply(seq_len(horizon), function(h) {
    rp_mat  <- fc_rp_arr[, h, ]
    zp_mat  <- fc_zp_arr[, h, ]
    sp_mat  <- fc_sp_arr[, h, ]
    pp_mat  <- fc_pp_arr[, h, ]
    draws_t <- fc_total_mat[, h]
    
    sum_rp <- rowSums(rp_mat)
    sum_zp <- rowSums(zp_mat)
    sum_sp <- rowSums(sp_mat)
    sum_pp <- rowSums(pp_mat)
    
    # CRITICAL FIX: Stacking order changed from (draws_t, pp_mat...) to match 
    # the exact total | pp | sp | zp | rp sequence used during the MinT calibration step.
    joint_mat     <- cbind(draws_t, pp_mat, sp_mat, zp_mat, rp_mat)
    actuals_joint <- c(actuals_total[h], actuals_pp[h, ], actuals_sp[h, ], actuals_zp[h, ], actuals_rp[h, ])
    
    data.frame(
      origin       = label_origin,
      horizon      = h,
      es           = compute_energy_score(joint_mat, actuals_joint,
                                          n_subsample = n_es_subsample),
      mai_rp_zp    = mean(sum_rp) - mean(sum_zp),
      mai_zp_sp    = mean(sum_zp) - mean(sum_sp),
      mai_sp_pp    = mean(sum_sp) - mean(sum_pp),
      mai_pp_total = mean(sum_pp) - mean(draws_t),
      mai_rp_total = mean(sum_rp) - mean(draws_t),
      jsd_rp_zp    = JSD_metric(sum_rp, sum_zp),
      jsd_zp_sp    = JSD_metric(sum_zp, sum_sp),
      jsd_sp_pp    = JSD_metric(sum_sp, sum_pp),
      jsd_pp_total = JSD_metric(sum_pp, draws_t),
      jsd_rp_total = JSD_metric(sum_rp, draws_t)
    )
  }))
  
  list(univariate = univariate_df, multivariate = multivariate_df)
}

# ── Compute and save base accuracy ─────────────────────────────────────────────

base_acc <- make_accuracy(
  fc_rp_orig, fc_zp_orig, fc_sp_orig, fc_pp_orig, fc_total_orig,
  label_origin = origin
)

saveRDS(
  base_acc,
  file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin))
)

# ── Compute and save MinT accuracy ─────────────────────────────────────────────

mint_acc <- make_accuracy(
  fc_rp_mint, fc_zp_mint, fc_sp_mint, fc_pp_mint, fc_total_mint,
  label_origin = origin
)

saveRDS(
  mint_acc,
  file.path(dir_mint_accuracy, sprintf("accuracy_timestep_%03d.rds", origin))
)