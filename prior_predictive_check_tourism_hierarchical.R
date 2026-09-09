# =============================================================================
# Prior predictive check — tourism ETS(A,N,A) hierarchy
# =============================================================================
# Reuses the same series-identification and standardization logic as the
# fitting script, but only to (a) pick representative nodes per level and
# (b) put the real data on the same standardized scale the priors are
# specified on. No S matrix / MinT machinery needed here.

library(dplyr)
library(stringr)
library(rstan)
library(ggplot2)

load("tourism.Rdata")

# =============================================================================
# 1. SERIES IDENTIFICATION (identical to the fitting script)
# =============================================================================

all_names <- names(tourism)

rp_codes <- all_names[grepl("^[A-G][A-Z][A-Z](Hol|Vis|Bus|Oth)$", all_names)]
zp_codes <- all_names[grepl("^[A-G][A-Z](Hol|Vis|Bus|Oth)$",      all_names)]
sp_codes <- all_names[grepl("^[A-G](Hol|Vis|Bus|Oth)$",           all_names)]
pp_codes <- all_names[grepl("^Total(Hol|Vis|Bus|Oth)$",           all_names)]

n_rp <- length(rp_codes); n_zp <- length(zp_codes)
n_sp <- length(sp_codes); n_pp <- length(pp_codes)

rp_to_zp <- match(paste0(substr(rp_codes, 1, 2), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), zp_codes)
rp_to_sp <- match(paste0(substr(rp_codes, 1, 1), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
rp_to_pp <- match(paste0("Total",                str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)

rp_per_zp <- as.integer(table(rp_to_zp)[order(as.integer(names(table(rp_to_zp))))])
rp_per_sp <- as.integer(table(rp_to_sp)[order(as.integer(names(table(rp_to_sp))))])
rp_per_pp <- as.integer(table(rp_to_pp)[order(as.integer(names(table(rp_to_pp))))])

message(sprintf("rp_per_zp range: %d - %d", min(rp_per_zp), max(rp_per_zp)))
message(sprintf("rp_per_sp range: %d - %d", min(rp_per_sp), max(rp_per_sp)))
message(sprintf("rp_per_pp range: %d - %d (constant means one pp node is enough)", min(rp_per_pp), max(rp_per_pp)))

# =============================================================================
# 2. PICK REPRESENTATIVE NODES
# =============================================================================
# rp: exactly exchangeable across nodes (identical priors) -> one is enough.
# pp: check the message above — if rp_per_pp is constant, one is enough here too.
# zp / sp: l0/s0 prior scale depends on rp_per_zp / rp_per_sp, which vary
# node-to-node, so we take the smallest- and largest-aggregation nodes to
# span the range rather than picking one arbitrarily.
set.seed(2025)

pick1 <- function(x) if (length(x) == 1) x else sample(x, 1)

purpose <- str_extract(sp_codes, "(Hol|Vis|Bus|Oth)$")
sp_to_pp <- match(paste0("Total", purpose), pp_codes)   # purpose-only match — no letter prefix
zp_to_sp <- match(paste0(substr(zp_codes, 1, 1), str_extract(zp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)

pp_chk_idx <- pick1(1:n_pp)

sp_children <- which(sp_to_pp == pp_chk_idx)
sp_chk_idx  <- pick1(sp_children)

zp_children <- which(zp_to_sp == sp_chk_idx)
zp_chk_idx  <- pick1(zp_children)

rp_children <- which(rp_to_zp == zp_chk_idx)
rp_chk_idx  <- if (length(rp_children) >= 2) sample(rp_children, 2) else rp_children

sp_to_pp_chk <- 1L
zp_to_sp_chk <- 1L
rp_to_zp_chk <- rep(1L, length(rp_chk_idx))

message(sprintf("Chain: pp=%s -> sp=%s -> zp=%s -> rp={%s}",
                pp_codes[pp_chk_idx], sp_codes[sp_chk_idx], zp_codes[zp_chk_idx],
                paste(rp_codes[rp_chk_idx], collapse = ", ")))

# =============================================================================
# 3. STANDARDIZE REAL DATA (same convention as the fitting script)
# =============================================================================

Y_raw_rp    <- tourism %>% select(all_of(rp_codes)) %>% as.matrix()
Y_raw_zp    <- tourism %>% select(all_of(zp_codes)) %>% as.matrix()
Y_raw_sp    <- tourism %>% select(all_of(sp_codes)) %>% as.matrix()
Y_raw_pp    <- tourism %>% select(all_of(pp_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

window_size <- 96
m           <- 12L
train_idx   <- 1:window_size    # representative window; swap in a real `origin` if you want a specific one

H <- 96L    # steps ahead to simulate for the prior predictive check

Y_mean <- mean(Y_raw_rp[train_idx, ])
Y_sd   <- sd(Y_raw_rp[train_idx, ])

Y_rp <- (Y_raw_rp[train_idx, ] - Y_mean) / Y_sd
Y_zp <- sweep(sweep(Y_raw_zp[train_idx, ], 2, Y_mean * rp_per_zp, "-"), 2, Y_sd, "/")
Y_sp <- sweep(sweep(Y_raw_sp[train_idx, ], 2, Y_mean * rp_per_sp, "-"), 2, Y_sd, "/")
Y_pp <- sweep(sweep(Y_raw_pp[train_idx, ], 2, Y_mean * rp_per_pp, "-"), 2, Y_sd, "/")
Y_total_std <- (Y_raw_total[train_idx] - Y_mean * n_rp) / Y_sd

# =============================================================================
# 4. SAMPLE THE PRIOR
# =============================================================================
# No likelihood, and only 1 + 2 + 2 + 1 + 1 = 7 representative series instead
# of n_rp + n_zp + n_sp + n_pp + 1 (likely several hundred) — this is the
# main lever that keeps this cheap relative to the fitting run. Compilation
# (one-off, ~1-2 min) dominates; sampling itself should be well under a
# minute for a few thousand draws given how little each replication costs
# (a handful of length-T state-space recursions, no posterior geometry to
# fight since there's no likelihood). Run the smoke test below first to get
# an actual per-iteration timing before committing to the full iter count.
model_string <- "

data {
  int<lower=1> T;
  int<lower=1> m;

  int<lower=1> n_pp_chk;
  int<lower=1> n_sp_chk;
  int<lower=1> n_zp_chk;
  int<lower=1> n_rp_chk;

  // mapping *within the representative subtree* (not the full tree) —
  // e.g. sp_to_pp_chk[1] = 1 means the one sp node's parent is pp_chk[1]
  array[n_sp_chk] int<lower=1, upper=n_pp_chk> sp_to_pp_chk;
  array[n_zp_chk] int<lower=1, upper=n_sp_chk> zp_to_sp_chk;
  array[n_rp_chk] int<lower=1, upper=n_zp_chk> rp_to_zp_chk;

  // l0/s0 prior scale still depends on how many rp roll up into each node
  array[n_pp_chk] int<lower=1> rp_per_pp_chk;
  array[n_sp_chk] int<lower=1> rp_per_sp_chk;
  array[n_zp_chk] int<lower=1> rp_per_zp_chk;

  int<lower=1> n_rp_total;   // total rp count in the FULL hierarchy
}

parameters {
  // ── Global hyper-means ──────────────────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_gamma_raw;

  // ── Level-specific hyper-means (each shrinks toward its parent's) ───────
  vector<lower=0, upper=1>[n_pp_chk] mu_alpha_pp;
  vector<lower=0, upper=1>[n_pp_chk] mu_gamma_raw_pp;

  vector<lower=0, upper=1>[n_sp_chk] mu_alpha_sp;
  vector<lower=0, upper=1>[n_sp_chk] mu_gamma_raw_sp;

  vector<lower=0, upper=1>[n_zp_chk] mu_alpha_zp;
  vector<lower=0, upper=1>[n_zp_chk] mu_gamma_raw_zp;

  // ── Region×purpose (bottom) ──────────────────────────────────────────────
  vector<lower=0, upper=1>[n_rp_chk] alpha_rp;
  vector<lower=0, upper=1>[n_rp_chk] gamma_raw_rp;
  vector<lower=0>[n_rp_chk]          sigma_rp;
  vector[n_rp_chk]                   l0_rp;
  matrix[n_rp_chk, m - 1]            s0_rp_free;

  // ── Zone×purpose ──────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_zp_chk] alpha_zp;
  vector<lower=0, upper=1>[n_zp_chk] gamma_raw_zp;
  vector<lower=0>[n_zp_chk]          sigma_zp;
  vector[n_zp_chk]                   l0_zp;
  matrix[n_zp_chk, m - 1]            s0_zp_free;

  // ── State×purpose ─────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_sp_chk] alpha_sp;
  vector<lower=0, upper=1>[n_sp_chk] gamma_raw_sp;
  vector<lower=0>[n_sp_chk]          sigma_sp;
  vector[n_sp_chk]                   l0_sp;
  matrix[n_sp_chk, m - 1]            s0_sp_free;

  // ── Purpose totals ────────────────────────────────────────────────────────
  vector<lower=0, upper=1>[n_pp_chk] alpha_pp;
  vector<lower=0, upper=1>[n_pp_chk] gamma_raw_pp;
  vector<lower=0>[n_pp_chk]          sigma_pp;
  vector[n_pp_chk]                   l0_pp;
  matrix[n_pp_chk, m - 1]            s0_pp_free;

  // ── Grand total ───────────────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_total;
  real<lower=0, upper=1> gamma_raw_total;
  real<lower=0>          sigma_total;
  real                   l0_total;
  vector[m - 1]          s0_total_free;
}

transformed parameters {
  matrix[n_rp_chk, m] s0_rp;
  matrix[n_zp_chk, m] s0_zp;
  matrix[n_sp_chk, m] s0_sp;
  matrix[n_pp_chk, m] s0_pp;
  vector[m]           s0_total;

  for (n in 1:n_rp_chk) { s0_rp[n, 1:(m-1)] = s0_rp_free[n]; s0_rp[n, m] = -sum(s0_rp_free[n]); }
  for (n in 1:n_zp_chk) { s0_zp[n, 1:(m-1)] = s0_zp_free[n]; s0_zp[n, m] = -sum(s0_zp_free[n]); }
  for (n in 1:n_sp_chk) { s0_sp[n, 1:(m-1)] = s0_sp_free[n]; s0_sp[n, m] = -sum(s0_sp_free[n]); }
  for (n in 1:n_pp_chk) { s0_pp[n, 1:(m-1)] = s0_pp_free[n]; s0_pp[n, m] = -sum(s0_pp_free[n]); }
  s0_total[1:(m-1)] = s0_total_free;
  s0_total[m]        = -sum(s0_total_free);

  vector[n_rp_chk] gamma_rp;
  vector[n_zp_chk] gamma_zp;
  vector[n_sp_chk] gamma_sp;
  vector[n_pp_chk] gamma_pp;
  real             gamma_total;

  for (n in 1:n_rp_chk) gamma_rp[n] = gamma_raw_rp[n] * (1 - alpha_rp[n]);
  for (n in 1:n_zp_chk) gamma_zp[n] = gamma_raw_zp[n] * (1 - alpha_zp[n]);
  for (n in 1:n_sp_chk) gamma_sp[n] = gamma_raw_sp[n] * (1 - alpha_sp[n]);
  for (n in 1:n_pp_chk) gamma_pp[n] = gamma_raw_pp[n] * (1 - alpha_pp[n]);
  gamma_total = gamma_raw_total * (1 - alpha_total);
}

model {
  // ── Same hyperprior chain as the fitting model, restricted to the
  //    representative subtree via the *_chk mapping arrays ────────────────
  mu_alpha     ~ uniform(0, 1);
  mu_gamma_raw ~ uniform(0, 1);

  mu_alpha_pp     ~ normal(mu_alpha,     0.1);
  mu_gamma_raw_pp ~ normal(mu_gamma_raw, 0.1);

  for (n in 1:n_sp_chk) {
    int p = sp_to_pp_chk[n];
    mu_alpha_sp[n]     ~ normal(mu_alpha_pp[p],     0.1);
    mu_gamma_raw_sp[n] ~ normal(mu_gamma_raw_pp[p], 0.1);
  }
  for (n in 1:n_zp_chk) {
    int s = zp_to_sp_chk[n];
    mu_alpha_zp[n]     ~ normal(mu_alpha_sp[s],     0.1);
    mu_gamma_raw_zp[n] ~ normal(mu_gamma_raw_sp[s], 0.1);
  }

  alpha_total     ~ normal(mu_alpha,     0.1);
  gamma_raw_total ~ normal(mu_gamma_raw, 0.1);
  sigma_total     ~ cauchy(0, 1);

  for (n in 1:n_pp_chk) {
    alpha_pp[n]     ~ normal(mu_alpha_pp[n],     0.1);
    gamma_raw_pp[n] ~ normal(mu_gamma_raw_pp[n], 0.1);
  }
  sigma_pp ~ cauchy(0, 0.5);

  for (n in 1:n_sp_chk) {
    int p = sp_to_pp_chk[n];
    alpha_sp[n]     ~ normal(mu_alpha_pp[p],     0.1);
    gamma_raw_sp[n] ~ normal(mu_gamma_raw_pp[p], 0.1);
  }
  sigma_sp ~ cauchy(0, 0.5);

  for (n in 1:n_zp_chk) {
    int s = zp_to_sp_chk[n];
    alpha_zp[n]     ~ normal(mu_alpha_sp[s],     0.1);
    gamma_raw_zp[n] ~ normal(mu_gamma_raw_sp[s], 0.1);
  }
  sigma_zp ~ cauchy(0, 0.3);

  for (n in 1:n_rp_chk) {
    int z = rp_to_zp_chk[n];
    alpha_rp[n]     ~ normal(mu_alpha_zp[z],     0.1);
    gamma_raw_rp[n] ~ normal(mu_gamma_raw_zp[z], 0.1);
  }
  sigma_rp ~ cauchy(0, 0.3);

  // ── l0/s0 priors, identical in form to the fitting model ─────────────────
  l0_rp ~ normal(0, 1);
  to_vector(s0_rp_free) ~ normal(0, 1);
  for (n in 1:n_zp_chk) { l0_zp[n] ~ normal(0, rp_per_zp_chk[n]); s0_zp_free[n] ~ normal(0, rp_per_zp_chk[n]); }
  for (n in 1:n_sp_chk) { l0_sp[n] ~ normal(0, rp_per_sp_chk[n]); s0_sp_free[n] ~ normal(0, rp_per_sp_chk[n]); }
  for (n in 1:n_pp_chk) { l0_pp[n] ~ normal(0, rp_per_pp_chk[n]); s0_pp_free[n] ~ normal(0, rp_per_pp_chk[n]); }
  l0_total      ~ normal(0, n_rp_total);
  s0_total_free ~ normal(0, n_rp_total);
}

generated quantities {
  matrix[T, n_rp_chk] Y_rp_sim;
  matrix[T, n_zp_chk] Y_zp_sim;
  matrix[T, n_sp_chk] Y_sp_sim;
  matrix[T, n_pp_chk] Y_pp_sim;
  vector[T]           Y_total_sim;

  // Same forward recursion as the independent-model PPC — the hierarchy
  // only affects how alpha/gamma/sigma/l0/s0 are drawn, not how a series
  // evolves once its own parameters are fixed.

  for (n in 1:n_rp_chk) {
    vector[T + 1] l_sim;
    vector[T + m] s_sim;
    for (j in 1:m) s_sim[j] = s0_rp[n, j];
    l_sim[1] = l0_rp[n];
    for (i in 1:T) {
      real e_i = normal_rng(0, sigma_rp[n]);
      Y_rp_sim[i, n] = l_sim[i] + s_sim[i] + e_i;
      l_sim[i+1]     = l_sim[i] + alpha_rp[n] * e_i;
      s_sim[i+m]     = s_sim[i] + gamma_rp[n] * e_i;
    }
  }

  for (n in 1:n_zp_chk) {
    vector[T + 1] l_sim;
    vector[T + m] s_sim;
    for (j in 1:m) s_sim[j] = s0_zp[n, j];
    l_sim[1] = l0_zp[n];
    for (i in 1:T) {
      real e_i = normal_rng(0, sigma_zp[n]);
      Y_zp_sim[i, n] = l_sim[i] + s_sim[i] + e_i;
      l_sim[i+1]     = l_sim[i] + alpha_zp[n] * e_i;
      s_sim[i+m]     = s_sim[i] + gamma_zp[n] * e_i;
    }
  }

  for (n in 1:n_sp_chk) {
    vector[T + 1] l_sim;
    vector[T + m] s_sim;
    for (j in 1:m) s_sim[j] = s0_sp[n, j];
    l_sim[1] = l0_sp[n];
    for (i in 1:T) {
      real e_i = normal_rng(0, sigma_sp[n]);
      Y_sp_sim[i, n] = l_sim[i] + s_sim[i] + e_i;
      l_sim[i+1]     = l_sim[i] + alpha_sp[n] * e_i;
      s_sim[i+m]     = s_sim[i] + gamma_sp[n] * e_i;
    }
  }

  for (n in 1:n_pp_chk) {
    vector[T + 1] l_sim;
    vector[T + m] s_sim;
    for (j in 1:m) s_sim[j] = s0_pp[n, j];
    l_sim[1] = l0_pp[n];
    for (i in 1:T) {
      real e_i = normal_rng(0, sigma_pp[n]);
      Y_pp_sim[i, n] = l_sim[i] + s_sim[i] + e_i;
      l_sim[i+1]     = l_sim[i] + alpha_pp[n] * e_i;
      s_sim[i+m]     = s_sim[i] + gamma_pp[n] * e_i;
    }
  }

  {
    vector[T + 1] l_sim;
    vector[T + m] s_sim;
    for (j in 1:m) s_sim[j] = s0_total[j];
    l_sim[1] = l0_total;
    for (i in 1:T) {
      real e_i = normal_rng(0, sigma_total);
      Y_total_sim[i] = l_sim[i] + s_sim[i] + e_i;
      l_sim[i+1]     = l_sim[i] + alpha_total * e_i;
      s_sim[i+m]     = s_sim[i] + gamma_total * e_i;
    }
  }
}
"

ppc_model <- stan_model(model_code = model_string)

ppc_data <- list(
  T = H,
  m = m,
  n_pp_chk = 1L, n_sp_chk = 1L, n_zp_chk = 1L, n_rp_chk = length(rp_chk_idx),
  sp_to_pp_chk = as.array(sp_to_pp_chk),
  zp_to_sp_chk = as.array(zp_to_sp_chk),
  rp_to_zp_chk = as.array(rp_to_zp_chk),
  rp_per_pp_chk = as.array(rp_per_pp[pp_chk_idx]),
  rp_per_sp_chk = as.array(rp_per_sp[sp_chk_idx]),
  rp_per_zp_chk = as.array(rp_per_zp[zp_chk_idx]),
  n_rp_total = n_rp
)
# as.array() matters here specifically because n_pp_chk = 1 (and n_rp_chk = 1):
# rstan drops the array dimension on length-1 numeric vectors when building
# the data list, so Stan sees a bare scalar where an array[1] is declared.
# zp_chk_idx/sp_chk_idx are length 2 here so wouldn't hit this, but as.array()
# is harmless either way and keeps this safe if you change n_zp_chk/n_sp_chk
# to 1 later.

# Smoke test — comment out once you've confirmed timing
smoke_start <- Sys.time()
smoke_fit <- sampling(ppc_model, data = ppc_data, chains = 1, iter = 300, warmup = 100, refresh = 0)
message(sprintf("Smoke test (300 iter, 1 chain): %.1f seconds", as.numeric(Sys.time() - smoke_start, units = "secs")))

ppc_fit <- sampling(
  object  = ppc_model,
  data    = ppc_data,
  seed    = 2025,
  chains  = 5,
  iter    = 2500,
  warmup  = 500,
  cores   = 5,
  refresh = 0
)

# =============================================================================
# 5. EXTRACT REPLICATIONS
# =============================================================================

Y_rp_sim    <- rstan::extract(ppc_fit, "Y_rp_sim")$Y_rp_sim        # [draws, T, 1]
Y_zp_sim    <- rstan::extract(ppc_fit, "Y_zp_sim")$Y_zp_sim        # [draws, T, 2]
Y_sp_sim    <- rstan::extract(ppc_fit, "Y_sp_sim")$Y_sp_sim        # [draws, T, 2]
Y_pp_sim    <- rstan::extract(ppc_fit, "Y_pp_sim")$Y_pp_sim        # [draws, T, 1]
Y_total_sim <- rstan::extract(ppc_fit, "Y_total_sim")$Y_total_sim  # [draws, T]

n_draws <- dim(Y_rp_sim)[1]
message(sprintf("Extracted %d replications", n_draws))

# =============================================================================
# 6. FAN CHART OVER THE H-STEP HORIZON
# =============================================================================
# Capping accumulation at H=12 steps instead of 96 keeps the bands in a
# readable range, so a direct ribbon/median-vs-actual plot is legible again
# here in a way it wasn't over the full training length.

qs <- c(0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995)

ribbon_from_draws <- function(sim_mat, label) {
  apply(sim_mat, 2, quantile, probs = qs) %>%
    t() %>%
    as.data.frame() %>%
    setNames(c("q005", "q025", "q05", "q25", "q50", "q75", "q95", "q975", "q995")) %>%
    mutate(t = row_number(), series = label)
}

# series_levels <- c(
#   "Total",
#   paste0("pp: ", pp_codes[pp_chk_idx[1]]),
#   paste0("sp: ", sp_codes[sp_chk_idx[1]], " (min rp)"),
#   paste0("sp: ", sp_codes[sp_chk_idx[2]], " (max rp)"),
#   paste0("zp: ", zp_codes[zp_chk_idx[1]], " (min rp)"),
#   paste0("zp: ", zp_codes[zp_chk_idx[2]], " (max rp)"),
#   paste0("rp: ", rp_codes[rp_chk_idx[1]])
# )
series_levels <- c(
  "Total",
  paste0("pp: ", pp_codes[pp_chk_idx[1]]),
  paste0("sp: ", sp_codes[sp_chk_idx[1]]),
  paste0("zp: ", zp_codes[zp_chk_idx[1]]),
  paste0("rp: ", rp_codes[rp_chk_idx[1]])
)

# ribbons <- bind_rows(
#   ribbon_from_draws(Y_total_sim,      series_levels[1]),
#   ribbon_from_draws(Y_pp_sim[, , 1],  series_levels[2]),
#   ribbon_from_draws(Y_sp_sim[, , 1],  series_levels[3]),
#   ribbon_from_draws(Y_sp_sim[, , 2],  series_levels[4]),
#   ribbon_from_draws(Y_zp_sim[, , 1],  series_levels[5]),
#   ribbon_from_draws(Y_zp_sim[, , 2],  series_levels[6]),
#   ribbon_from_draws(Y_rp_sim[, , 1],  series_levels[7])
# ) %>%
#   mutate(series = factor(series, levels = series_levels))
ribbons <- bind_rows(
  ribbon_from_draws(Y_total_sim,      series_levels[1]),
  ribbon_from_draws(Y_pp_sim[, , 1],  series_levels[2]),
  ribbon_from_draws(Y_sp_sim[, , 1],  series_levels[3]),
  ribbon_from_draws(Y_zp_sim[, , 1],  series_levels[4]),
  ribbon_from_draws(Y_rp_sim[, , 1],  series_levels[5])
) %>%
  mutate(series = factor(series, levels = series_levels))

# Overlay the first H periods of the real (standardized) data — same starting
# point (t=1) the prior on l0/s0 is centred around.
actual <- bind_rows(
  data.frame(t = seq_len(H), y = Y_total_std[seq_len(H)],              series = series_levels[1]),
  data.frame(t = seq_len(H), y = Y_pp[seq_len(H), pp_chk_idx[1]],       series = series_levels[2]),
  data.frame(t = seq_len(H), y = Y_sp[seq_len(H), sp_chk_idx[1]],       series = series_levels[3]),
  data.frame(t = seq_len(H), y = Y_zp[seq_len(H), zp_chk_idx[1]],       series = series_levels[4]),
  data.frame(t = seq_len(H), y = Y_rp[seq_len(H), rp_chk_idx[1]],       series = series_levels[5])
) %>%
  mutate(series = factor(series, levels = series_levels))

# Define manual fill scale palette mapping interval labels to colors/alphas
interval_colors <- c(
  "99% Interval" = alpha("steelblue", 0.15),
  "95% Interval" = alpha("steelblue", 0.25),
  "90% Interval" = alpha("steelblue", 0.40),
  "50% Interval" = alpha("steelblue", 0.60)
)

fan_plot <- ggplot(ribbons, aes(x = t)) +
  # Map fill aesthetic to strings so ggplot builds the discrete legend
  geom_ribbon(aes(ymin = q005, ymax = q995, fill = "99% Interval")) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = "95% Interval")) +
  geom_ribbon(aes(ymin = q05,  ymax = q95,  fill = "90% Interval")) +
  geom_ribbon(aes(ymin = q25,  ymax = q75,  fill = "50% Interval")) +
  # Line layers with legend key for the actual path
  geom_line(aes(y = q50, colour = "Prior Median"), linewidth = 0.6) +
  geom_line(data = actual, aes(x = t, y = y, colour = "Observed Data"), linewidth = 0.6) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  scale_fill_manual(
    name = "Predictive Intervals",
    values = interval_colors,
    breaks = c("99% Interval", "95% Interval", "90% Interval", "50% Interval")
  ) +
  scale_colour_manual(
    name = "Series",
    values = c("Prior Median" = "steelblue4", "Observed Data" = "black")
  ) +
  labs(
    title    = "Prior predictive check: Replications vs. observed data",
    subtitle = sprintf("%d prior predictive replications, length %d", n_draws, H),
    x = "t", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.margin = margin(t = -5)
  )

ggsave("prior_predictive_tourism_h12_fanchart.png", fan_plot, width = 10, height = 10, dpi = 150)

message("Saved: prior_predictive_tourism_h12_fanchart.png")
