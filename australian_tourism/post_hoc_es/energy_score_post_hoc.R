library(dplyr)
library(stringr)
library(rstan)
library(scoringRules)

# ── 1. Setup Directories & Load Data ──────────────────────────────────────────
load("tourism.Rdata")
job_name <- "tourism_job_22307354"  # soft coherency

dir_stanfits <- file.path(job_name, "stanfits")
dir_accuracy <- file.path(job_name, "es_post_hoc")

if (!dir.exists(dir_accuracy)) dir.create(dir_accuracy, recursive = TRUE)

all_names <- names(tourism)

rp_codes <- all_names[grepl("^[A-G][A-Z][A-Z](Hol|Vis|Bus|Oth)$", all_names)]
zp_codes <- all_names[grepl("^[A-G][A-Z](Hol|Vis|Bus|Oth)$",      all_names)]
sp_codes <- all_names[grepl("^[A-G](Hol|Vis|Bus|Oth)$",           all_names)]
pp_codes <- all_names[grepl("^Total(Hol|Vis|Bus|Oth)$",           all_names)]

n_rp  <- length(rp_codes)
n_zp  <- length(zp_codes)
n_sp  <- length(sp_codes)
n_pp  <- length(pp_codes)

# Hierarchical sub-counts for scaling
rp_to_zp <- match(paste0(substr(rp_codes, 1, 2), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), zp_codes)
rp_to_sp <- match(paste0(substr(rp_codes, 1, 1), str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), sp_codes)
rp_to_pp <- match(paste0("Total", str_extract(rp_codes, "(Hol|Vis|Bus|Oth)$")), pp_codes)

rp_per_zp <- as.integer(table(rp_to_zp)[order(as.integer(names(table(rp_to_zp))))])
rp_per_sp <- as.integer(table(rp_to_sp)[order(as.integer(names(table(rp_to_sp))))])
rp_per_pp <- as.integer(table(rp_to_pp)[order(as.integer(names(table(rp_to_pp))))])

# Raw data matrices
Y_raw_rp    <- tourism %>% select(all_of(rp_codes)) %>% as.matrix()
Y_raw_zp    <- tourism %>% select(all_of(zp_codes)) %>% as.matrix()
Y_raw_sp    <- tourism %>% select(all_of(sp_codes)) %>% as.matrix()
Y_raw_pp    <- tourism %>% select(all_of(pp_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

# Configurations
window_size    <- 96L
horizon        <- 12L
test_set       <- 121L
origins        <- 96:(96 + test_set - 1)
n_es_subsample <- 1000L   

# ── 2. Helpers (Energy Score & JSD Metric) ────────────────────────────────────
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

# ── 3. Main Loop over Forecast Origins ────────────────────────────────────────
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
  
  # Single-pass parsing (No fitted parameters requested)
  all_pars <- rstan::extract(stanFit, pars = c("Y_fc_rp", "Y_fc_zp", "Y_fc_sp", "Y_fc_pp", "Y_fc_total"))
  n_draws  <- nrow(all_pars$Y_fc_rp)
  
  # ── Safe Unscaling of Forecast Draws ────────────────────────────────────────
  fc_total_orig <- all_pars$Y_fc_total[, 1:horizon, drop = FALSE] * Y_sd + Y_mean * n_rp
  fc_rp_orig    <- all_pars$Y_fc_rp[, 1:horizon, , drop = FALSE] * Y_sd + Y_mean
  
  fc_zp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_zp))
  for (z in seq_len(n_zp)) {
    fc_zp_orig[, , z] <- matrix(all_pars$Y_fc_zp[, 1:horizon, z], nrow = n_draws, ncol = horizon) * Y_sd + Y_mean * rp_per_zp[z]
  }
  
  fc_sp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_sp))
  for (s in seq_len(n_sp)) {
    fc_sp_orig[, , s] <- matrix(all_pars$Y_fc_sp[, 1:horizon, s], nrow = n_draws, ncol = horizon) * Y_sd + Y_mean * rp_per_sp[s]
  }
  
  fc_pp_orig <- array(NA_real_, dim = c(n_draws, horizon, n_pp))
  for (p in seq_len(n_pp)) {
    fc_pp_orig[, , p] <- matrix(all_pars$Y_fc_pp[, 1:horizon, p], nrow = n_draws, ncol = horizon) * Y_sd + Y_mean * rp_per_pp[p]
  }
  
  rm(all_pars); gc()
  
  # ── Slice Target Actuals ───────────────────────────────────────────────────
  actuals_rp    <- Y_raw_rp[test_idx, , drop = FALSE]
  actuals_zp    <- Y_raw_zp[test_idx, , drop = FALSE]
  actuals_sp    <- Y_raw_sp[test_idx, , drop = FALSE]
  actuals_pp    <- Y_raw_pp[test_idx, , drop = FALSE]
  actuals_total <- Y_raw_total[test_idx]
  
  # ── 4. Core Metrics Calculations ───────────────────────────────────────────
  
  # Helper for univariate mapping
  make_rows <- function(draws_arr, actuals_mat, codes, level) {
    lapply(seq_len(ncol(actuals_mat)), function(i) {
      lapply(seq_len(horizon), function(h) {
        draws <- draws_arr[, h, i]
        actual <- actuals_mat[h, i]
        data.frame(
          origin = origin, horizon = h, series = codes[i], level = level, actual = actual,
          fc_mean = mean(draws), ae = abs(mean(draws) - actual), se = (mean(draws) - actual)^2,
          crps = scoringRules::crps_sample(actual, draws), logS = scoringRules::logs_sample(actual, draws),
          covered_90 = as.integer(actual >= quantile(draws, 0.05) & actual <= quantile(draws, 0.95)),
          covered_95 = as.integer(actual >= quantile(draws, 0.025) & actual <= quantile(draws, 0.975)),
          pit = mean(draws < actual)
        )
      })
    })
  }
  
  total_rows <- lapply(seq_len(horizon), function(h) {
    draws <- fc_total_orig[, h]; actual <- actuals_total[h]
    data.frame(
      origin = origin, horizon = h, series = "Total", level = "Total", actual = actual,
      fc_mean = mean(draws), ae = abs(mean(draws) - actual), se = (mean(draws) - actual)^2,
      crps = scoringRules::crps_sample(actual, draws), logS = scoringRules::logs_sample(actual, draws),
      covered_90 = as.integer(actual >= quantile(draws, 0.05) & actual <= quantile(draws, 0.95)),
      covered_95 = as.integer(actual >= quantile(draws, 0.025) & actual <= quantile(draws, 0.975)),
      pit = mean(draws < actual)
    )
  })
  
  univariate_df <- bind_rows(c(
    unlist(make_rows(fc_rp_orig, actuals_rp, rp_codes, "RegionPurpose"), recursive = FALSE),
    unlist(make_rows(fc_zp_orig, actuals_zp, zp_codes, "ZonePurpose"),   recursive = FALSE),
    unlist(make_rows(fc_sp_orig, actuals_sp, sp_codes, "StatePurpose"),  recursive = FALSE),
    unlist(make_rows(fc_pp_orig, actuals_pp, pp_codes, "Purpose"),       recursive = FALSE),
    total_rows
  ))
  
  # Original seed reset placement (outside horizon loop to match historical drift)
  set.seed(origin)
  
  multivariate_df <- bind_rows(lapply(seq_len(horizon), function(h) {
    # Drop the middle horizon dimension explicitly to create clean 2D matrices [n_draws, n_series]
    rp_mat  <- as.matrix(fc_rp_orig[, h, ])
    zp_mat  <- as.matrix(fc_zp_orig[, h, ])
    sp_mat  <- as.matrix(fc_sp_orig[, h, ])
    pp_mat  <- as.matrix(fc_pp_orig[, h, ])
    draws_t <- as.matrix(fc_total_orig[, h]) # Ensures total is a [n_draws, 1] matrix
    
    sum_rp <- rowSums(rp_mat); sum_zp <- rowSums(zp_mat)
    sum_sp <- rowSums(sp_mat); sum_pp <- rowSums(pp_mat)
    
    # Structural topology binding: total | pp | sp | zp | rp
    joint_mat     <- cbind(draws_t, pp_mat, sp_mat, zp_mat, rp_mat)
    actuals_joint <- c(actuals_total[h], actuals_pp[h, ], actuals_sp[h, ], actuals_zp[h, ], actuals_rp[h, ])
    
    data.frame(
      origin = origin, horizon = h,
      es = compute_energy_score(joint_mat, actuals_joint, n_subsample = n_es_subsample),
      mai_rp_zp = mean(sum_rp) - mean(sum_zp), mai_zp_sp = mean(sum_zp) - mean(sum_sp),
      mai_sp_pp = mean(sum_sp) - mean(sum_pp), mai_pp_total = mean(sum_pp) - mean(draws_t),
      mai_rp_total = mean(sum_rp) - mean(draws_t),
      jsd_rp_zp = JSD_metric(sum_rp, sum_zp), jsd_zp_sp = JSD_metric(sum_zp, sum_sp),
      jsd_sp_pp = JSD_metric(sum_sp, sum_pp), jsd_pp_total = JSD_metric(sum_pp, draws_t),
      jsd_rp_total = JSD_metric(sum_rp, draws_t)
    )
  }))
  
  # Save clean data list
  base_acc <- list(univariate = univariate_df, multivariate = multivariate_df)
  saveRDS(base_acc, file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin)))
  
  rm(stanFit, fc_rp_orig, fc_zp_orig, fc_sp_orig, fc_pp_orig, fc_total_orig)
  gc()
}