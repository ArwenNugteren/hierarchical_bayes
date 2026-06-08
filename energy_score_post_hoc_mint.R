library(dplyr)
library(stringr)
library(rstan)

# ── 1. Setup Directories & Load Data ──────────────────────────────────────────
load("tourism.Rdata")
job_name <- "tourism_job_21996540"

dir_stanfits    <- file.path(job_name, "stanfits")
dir_output_base <- file.path(job_name, "post_hoc_es_base")
dir_output_mint <- file.path(job_name, "post_hoc_es_mint")

for (d in c(dir_output_base, dir_output_mint)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

all_names <- names(tourism)

rp_codes <- all_names[grepl("^[A-G][A-Z][A-Z](Hol|Vis|Bus|Oth)$", all_names)]
zp_codes <- all_names[grepl("^[A-G][A-Z](Hol|Vis|Bus|Oth)$",      all_names)]
sp_codes <- all_names[grepl("^[A-G](Hol|Vis|Bus|Oth)$",           all_names)]
pp_codes <- all_names[grepl("^Total(Hol|Vis|Bus|Oth)$",           all_names)]

n_rp  <- length(rp_codes)
n_zp  <- length(zp_codes)
n_sp  <- length(sp_codes)
n_pp  <- length(pp_codes)
n_all <- 1L + n_pp + n_sp + n_zp + n_rp

# Structural hierarchy mappings
rp_to_zp <- match(paste0(substr(rp_codes, 1, 2), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), zp_codes)
rp_to_sp <- match(paste0(substr(rp_codes, 1, 1), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
rp_to_pp <- match(paste0("Total", str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)

rp_per_zp <- as.integer(table(rp_to_zp)[order(as.integer(names(table(rp_to_zp))))])
rp_per_sp <- as.integer(table(rp_to_sp)[order(as.integer(names(table(rp_to_sp))))])
rp_per_pp <- as.integer(table(rp_to_pp)[order(as.integer(names(table(rp_to_pp))))])

# Raw data
Y_raw_rp    <- tourism %>% select(all_of(rp_codes)) %>% as.matrix()
Y_raw_zp    <- tourism %>% select(all_of(zp_codes)) %>% as.matrix()
Y_raw_sp    <- tourism %>% select(all_of(sp_codes)) %>% as.matrix()
Y_raw_pp    <- tourism %>% select(all_of(pp_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

# Index positions in joint vector (Aligned with build_S structure)
idx_total <- 1L
idx_pp    <- 1L + seq_len(n_pp)
idx_sp    <- 1L + n_pp + seq_len(n_sp)
idx_zp    <- 1L + n_pp + n_sp + seq_len(n_zp)
idx_rp    <- 1L + n_pp + n_sp + n_zp + seq_len(n_rp)

# Rolling window configuration
window_size    <- 96L
horizon        <- 12L
test_set       <- 121L
origins        <- 96:(96 + test_set - 1)
n_es_subsample <- 1000L   

# ── 2. Energy Score Function (No Scaling, Split-Sample) ───────────────────────
compute_energy_score <- function(draws_mat, actuals, n_subsample = NULL) {
  n_draws <- nrow(draws_mat)
  if (!is.null(n_subsample) && n_subsample < n_draws) {
    idx1 <- sample.int(n_draws, size = n_subsample, replace = FALSE)
    idx2 <- sample.int(n_draws, size = n_subsample, replace = FALSE)
    x    <- draws_mat[idx1, , drop = FALSE]
    xs   <- draws_mat[idx2, , drop = FALSE]
  } else {
    x    <- draws_mat
    xs   <- draws_mat[sample.int(n_draws), , drop = FALSE] 
  }
  
  diffs <- sweep(x, 2, as.numeric(actuals), "-")
  t1    <- mean(sqrt(rowSums(diffs^2)))
  t2    <- mean(sqrt(rowSums((x - xs)^2)))
  
  return(t1 - 0.5 * t2)
}

# ── 3. Ledoit-Wolf Shrinkage Covariance ───────────────────────────────────────
lw_shrink <- function(E_mat) {
  TT      <- nrow(E_mat)
  Sc      <- cov(E_mat)
  d       <- diag(Sc)
  E_c     <- sweep(E_mat, 2, colMeans(E_mat), "-")
  xs      <- E_c^2
  xss     <- crossprod(E_c)
  var_sij <- (TT / (TT - 1)^3) * (crossprod(xs) - (1 / TT) * xss^2)
  denom   <- sum(xss[upper.tri(xss)]^2) / (TT - 1)^2
  rho_hat <- if (denom < .Machine$double.eps) 0 else min(1, sum(var_sij[upper.tri(var_sij)]) / denom)
  (1 - rho_hat) * Sc + rho_hat * diag(d)
}

# ── 4. S matrix (Explicit structural layout matching build_S) ─────────────────
build_S <- function() {
  S <- matrix(0L, nrow = n_all, ncol = n_rp)
  S[idx_total, ] <- 1L
  for (p in seq_len(n_pp)) S[idx_pp[p], rp_to_pp == p] <- 1L
  for (s in seq_len(n_sp)) S[idx_sp[s], rp_to_sp == s] <- 1L
  for (z in seq_len(n_zp)) S[idx_zp[z], rp_to_zp == z] <- 1L
  for (r in seq_len(n_rp)) S[idx_rp[r], r]             <- 1L
  S
}

S_mat <- build_S()
storage.mode(S_mat) <- "double"

# ── 5. Main Loop ──────────────────────────────────────────────────────────────
for (origin in origins) {
  
  stanfit_path <- file.path(dir_stanfits, sprintf("stanfit_timestep_%03d.rds", origin))
  
  if (!file.exists(stanfit_path)) {
    warning(sprintf("Stanfit missing for origin %d, skipping.", origin))
    next
  }
  
  cat(sprintf("Processing origin %d...\n", origin))
  stanFit   <- readRDS(stanfit_path)
  train_idx <- (origin - window_size + 1):origin
  test_idx  <- (origin + 1):(origin + horizon)
  
  Y_mean <- mean(Y_raw_rp[train_idx, ])
  Y_sd   <- sd(Y_raw_rp[train_idx, ])
  
  # Single-pass parsing of all required extractions
  all_pars <- rstan::extract(stanFit, pars = c(
    "Y_fc_rp", "Y_fc_zp", "Y_fc_sp", "Y_fc_pp", "Y_fc_total",
    "fitted_rp", "fitted_zp", "fitted_sp", "fitted_pp", "fitted_total"
  ))
  
  n_draws <- nrow(all_pars$Y_fc_rp)
  
  # ── In-sample residuals for MinT W matrix ──────────────────────────────────
  fit_rp_orig    <- apply(all_pars$fitted_rp, 2:3, mean) * Y_sd + Y_mean
  fit_total_orig <- apply(all_pars$fitted_total, 2, mean) * Y_sd + Y_mean * n_rp
  
  fit_zp_orig <- apply(all_pars$fitted_zp, 2:3, mean)
  for (z in seq_len(n_zp)) fit_zp_orig[, z] <- fit_zp_orig[, z] * Y_sd + Y_mean * rp_per_zp[z]
  
  fit_sp_orig <- apply(all_pars$fitted_sp, 2:3, mean)
  for (s in seq_len(n_sp)) fit_sp_orig[, s] <- fit_sp_orig[, s] * Y_sd + Y_mean * rp_per_sp[s]
  
  fit_pp_orig <- apply(all_pars$fitted_pp, 2:3, mean)
  for (p in seq_len(n_pp)) fit_pp_orig[, p] <- fit_pp_orig[, p] * Y_sd + Y_mean * rp_per_pp[p]
  
  # Structural column binding: total | pp | sp | zp | rp
  E <- cbind(
    Y_raw_total[train_idx] - fit_total_orig,
    Y_raw_pp[train_idx, ]  - fit_pp_orig,
    Y_raw_sp[train_idx, ]  - fit_sp_orig,
    Y_raw_zp[train_idx, ]  - fit_zp_orig,
    Y_raw_rp[train_idx, ]  - fit_rp_orig
  )
  
  # ── MinT Projection Matrix Matrix Operators ─────────────────────────────────
  W_mat <- lw_shrink(E)
  
  WinvS <- tryCatch(
    solve(W_mat, S_mat),
    error = function(e) solve(W_mat + diag(1e-6 * mean(diag(W_mat)), n_all), S_mat)
  )
  StWinvS_chol <- chol(crossprod(S_mat, WinvS))
  
  # Global Projection matrix P (n_all x n_all) computed once per origin
  P_mint <- S_mat %*% chol2inv(StWinvS_chol) %*% t(WinvS)
  
  # ── Safe Unscaling of Forecast Draws ────────────────────────────────────────
  fc_rp_orig    <- all_pars$Y_fc_rp[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean
  fc_total_orig <- all_pars$Y_fc_total[, 1:horizon, drop = FALSE] * Y_sd + Y_mean * n_rp
  
  fc_zp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_zp))
  for (z in seq_len(n_zp)) {
    fc_zp_orig[, , z] <- all_pars$Y_fc_zp[, 1:horizon, z] * Y_sd + Y_mean * rp_per_zp[z]
  }
  
  fc_sp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_sp))
  for (s in seq_len(n_sp)) {
    fc_sp_orig[, , s] <- all_pars$Y_fc_sp[, 1:horizon, s] * Y_sd + Y_mean * rp_per_sp[s]
  }
  
  fc_pp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_pp))
  for (p in seq_len(n_pp)) {
    fc_pp_orig[, , p] <- all_pars$Y_fc_pp[, 1:horizon, p] * Y_sd + Y_mean * rp_per_pp[p]
  }
  
  # ── Horizon Evaluation Loop ─────────────────────────────────────────────────
  set.seed(origin)
  base_es_list <- vector("list", horizon)
  mint_es_list <- vector("list", horizon)
  
  for (h in seq_len(horizon)) {
    
    # Structural topology array mapping (total | pp | sp | zp | rp)
    joint_base <- cbind(
      fc_total_orig[, h],
      fc_pp_orig[, h, ],
      fc_sp_orig[, h, ],
      fc_zp_orig[, h, ],
      fc_rp_orig[, h, ]
    )
    
    # Vectorized MinT transformation over all draws at once
    joint_mint <- joint_base %*% t(P_mint)
    
    actuals_joint <- c(
      Y_raw_total[test_idx[h]],
      Y_raw_pp[test_idx[h], ],
      Y_raw_sp[test_idx[h], ],
      Y_raw_zp[test_idx[h], ],
      Y_raw_rp[test_idx[h], ]
    )
    
    base_es_list[[h]] <- data.frame(
      origin       = origin,
      horizon      = h,
      energy_score = compute_energy_score(joint_base, actuals_joint, n_subsample = n_es_subsample)
    )
    mint_es_list[[h]] <- data.frame(
      origin       = origin,
      horizon      = h,
      energy_score = compute_energy_score(joint_mint, actuals_joint, n_subsample = n_es_subsample)
    )
  }
  
  saveRDS(bind_rows(base_es_list), file.path(dir_output_base, sprintf("es_base_origin_%03d.rds", origin)))
  saveRDS(bind_rows(mint_es_list), file.path(dir_output_mint, sprintf("es_mint_origin_%03d.rds", origin)))
  
  # Explicit memory clear
  rm(stanFit, all_pars, fc_rp_orig, fc_zp_orig, fc_sp_orig, fc_pp_orig, fc_total_orig, E, W_mat, P_mint)
  gc()
}

cat("Post-hoc ES evaluation complete.\n")