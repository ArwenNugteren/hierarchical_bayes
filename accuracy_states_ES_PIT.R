library(dplyr)
library(rstan)
library(ggplot2)
library(patchwork)
library(kldest)

# ── User-set paths and parameters ──────────────────────────────────────────────

dir_stanfits    <- "tourism_job_21345282/stanfits" 
dir_accuracy    <- "tourism_job_21345282/accuracy"    
dir_plots       <- "tourism_job_21345282/accuracy/plots" 

n_es_subsample <- 1000L                       # draws used for ES pairwise term

# ── Create output directories if needed ───────────────────────────────────────
dir.create(dir_plots,    recursive = TRUE, showWarnings = FALSE)

# ── Data and constants (must match the estimation script) ──────────────────────

load("tourism.Rdata")

state_codes <- names(tourism)[grepl("^[A-G]$", names(tourism))]
n_state     <- length(state_codes)
Y_raw_state <- tourism %>% select(all_of(state_codes)) %>% as.matrix()
Y_raw_total <- tourism %>% pull(Total)

m           <- 12L
window_size <- 96L
horizon     <- 12L
test_set    <- 121L
origins     <- 96:(96 + test_set - 1)

all_series  <- c(state_codes, "Total")        # 8-vector order used throughout
n_series    <- length(all_series)             # 8

# ── Helpers ────────────────────────────────────────────────────────────────────

# energy score calculation
energy_score <- function(draws_mat, actuals, sub_idx = NULL) {
  # 1. subset data (avoids several million calculations)
  if (!is.null(sub_idx)) {
    draws_mat <- draws_mat[sub_idx, , drop = FALSE]
  }
  
  # 2. Standardise columns by their cross-draw SD to maintain consistent scaling
  col_sds    <- apply(draws_mat, 2, sd)
  draws_sc   <- sweep(draws_mat, 2, col_sds, "/")
  
  # Ensure actuals is a flat vector matching the columns
  actuals_sc <- as.numeric(actuals) / col_sds
  
  # 3. Term 1: Distance from draws to actuals
  diffs1 <- sweep(draws_sc, 2, actuals_sc, "-")
  t1     <- mean(sqrt(rowSums(diffs1^2)))
  
  # 4. Term 2: Pairwise distance between draws (Memory efficient)
  # dist() calculates the exact Euclidean distances between all rows
  pairwise_dists <- dist(draws_sc, method = "euclidean")
  t2             <- mean(pairwise_dists)
  
  # Return standard Energy Score formulation (t1 - 0.5 * t2)
  return(t1 - 0.5 * t2)
}

# PIT jitter (breaks ties for KS test using small jitter)
jitter_pit <- function(x, amount = 1e-4) {
  pmax(0, pmin(1, x + runif(length(x), -amount, amount)))
}

# ── Per-origin accuracy loop ───────────────────────────────────────────────────

stanfit_files <- list.files(dir_stanfits,
                            pattern    = "^stanfit_timestep_[0-9]+\\.rds$",
                            full.names = TRUE)

for (f in stanfit_files) {
  
  # derive origin from filename
  origin <- as.integer(regmatches(f, regexpr("[0-9]+(?=\\.rds$)", f, perl = TRUE)))
  
  train_idx <- (origin - window_size + 1):origin
  test_idx  <- (origin + 1):(origin + horizon)
  
  # back-transform constants for this window
  Y_mean_w <- mean(Y_raw_state[train_idx, ])
  Y_sd_w   <- sd(Y_raw_state[train_idx, ])
  
  actuals_state <- Y_raw_state[test_idx, ]      # horizon x n_state
  actuals_total <- Y_raw_total[test_idx]        # horizon
  
  stanFit <- readRDS(f)
  
  # raw draws: (n_draws x horizon x n_state) and (n_draws x horizon)
  fc_state_draws <- rstan::extract(stanFit, pars = "Y_fc_state")$Y_fc_state
  fc_total_draws <- rstan::extract(stanFit, pars = "Y_fc_total")$Y_fc_total
  
  n_draws <- nrow(fc_state_draws)
  
  # back-transform to original scale
  fc_state_orig <- fc_state_draws[, 1:horizon, , drop = FALSE] * Y_sd_w + Y_mean_w
  fc_total_orig <- fc_total_draws[, 1:horizon,  drop = FALSE]  * Y_sd_w + Y_mean_w * n_state
  
  # fixed subsample index for ES pairwise term
  set.seed(origin)
  sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)
  
  horizon_results <- lapply(1:horizon, function(h) {
    
    # ── Univariate summaries per series ───────────────────────────────────────
    
    series_rows <- lapply(seq_len(n_state), function(s) {
      draws  <- fc_state_orig[, h, s]
      actual <- actuals_state[h, s]
      data.frame(
        origin  = origin,
        horizon = h,
        series  = state_codes[s],
        actual  = actual,
        fc_mean = mean(draws),
        ae      = abs(mean(draws) - actual),
        se      = (mean(draws) - actual)^2,
        pit     = mean(draws < actual)
      )
    })
    
    draws_t  <- fc_total_orig[, h]
    actual_t <- actuals_total[h]
    
    total_row <- data.frame(
      origin  = origin,
      horizon = h,
      series  = "Total",
      actual  = actual_t,
      fc_mean = mean(draws_t),
      ae      = abs(mean(draws_t) - actual_t),
      se      = (mean(draws_t) - actual_t)^2,
      pit     = mean(draws_t < actual_t)
    )
    
    univariate_df <- bind_rows(c(series_rows, list(total_row)))
    
    # ── Energy Score ──────────────────────────────────────────────────────────
    
    state_mat     <- fc_state_orig[, h, ]           # n_draws x n_state
    joint_mat     <- cbind(state_mat, draws_t)      # n_draws x n_series
    actuals_joint <- c(actuals_state[h, ], actual_t)
    
    es_val <- energy_score(joint_mat, actuals_joint, sub_idx)
    
    # ── Incoherence ───────────────────────────────────────────────────────────
    # mean forecast incoherence: E[total] minus E[sum of states]
    sum_states_draws <- rowSums(state_mat)           # n_draws
    incoherence_mean <- mean(draws_t) - mean(sum_states_draws)
    
    # ── KL divergence: KL(total || sum of states) ─────────────────────────────
    # p = total forecast distribution
    # q = sum-of-states forecast distribution (approximating distribution)
    # direction: KL(total || sum_states)
    kl_val <- tryCatch({
      kldest::kld_est_nn(
        matrix(draws_t,          ncol = 1),   # p: total draws
        matrix(sum_states_draws, ncol = 1)    # q: sum-of-states draws
      )
    }, error = function(e) NA_real_)
    
    list(
      univariate   = univariate_df,
      multivariate = data.frame(
        origin           = origin,
        horizon          = h,
        es               = es_val,
        incoherence_mean = incoherence_mean,
        kl_total_vs_sum  = kl_val
      )
    )
  })
  
  origin_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
  origin_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))
  
  saveRDS(
    list(univariate = origin_univariate, multivariate = origin_multivariate),
    file.path(dir_accuracy, sprintf("accuracy_timestep_%03d.rds", origin))
  )
  
  message(sprintf("Origin %d done.", origin))
}

# ── Aggregate ──────────────────────────────────────────────────────────────────

acc_files <- list.files(dir_accuracy,
                        pattern    = "^accuracy_timestep_[0-9]+\\.rds$",
                        full.names = TRUE)

all_uni   <- bind_rows(lapply(acc_files, function(f) readRDS(f)$univariate))
all_multi <- bind_rows(lapply(acc_files, function(f) readRDS(f)$multivariate))

# ── summary_by_horizon ────────────────────────────────────────────────────────
# MAE/RMSE averaged over all series and origins at each horizon
# ES, MAI, KL are joint/multivariate measures summarised per horizon

summary_by_horizon <- all_uni %>%
  group_by(horizon) %>%
  summarise(MAE = mean(ae), RMSE = sqrt(mean(se)), .groups = "drop") %>%
  left_join(
    all_multi %>%
      group_by(horizon) %>%
      summarise(
        ES              = mean(es),
        MAI             = mean(abs(incoherence_mean)),
        KL_total_vs_sum = mean(kl_total_vs_sum, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "horizon"
  )

# ── summary_by_series_horizon ─────────────────────────────────────────────────

summary_by_series_horizon <- all_uni %>%
  group_by(series, horizon) %>%
  summarise(MAE = mean(ae), RMSE = sqrt(mean(se)), .groups = "drop")

# ── summary_overall: rows = States / Total / Overall ─────────────────────────
# ES, MAI, KL are joint measures so only meaningful at the Overall level

summary_overall <- bind_rows(
  
  # States: MAE/RMSE for state series only; joint metrics NA
  all_uni %>%
    filter(series != "Total") %>%
    summarise(
      level           = "States",
      MAE             = mean(ae),
      RMSE            = sqrt(mean(se)),
      ES              = NA_real_,
      MAI             = NA_real_,
      KL_total_vs_sum = NA_real_
    ),
  
  # Total: MAE/RMSE for total series only; joint metrics NA
  all_uni %>%
    filter(series == "Total") %>%
    summarise(
      level           = "Total",
      MAE             = mean(ae),
      RMSE            = sqrt(mean(se)),
      ES              = NA_real_,
      MAI             = NA_real_,
      KL_total_vs_sum = NA_real_
    ),
  
  # Overall: MAE/RMSE over everything; ES, MAI, KL over all origins x horizons
  bind_cols(
    all_uni %>%
      summarise(
        level = "Overall",
        MAE   = mean(ae),
        RMSE  = sqrt(mean(se))
      ),
    all_multi %>%
      summarise(
        ES              = mean(es),
        MAI             = mean(abs(incoherence_mean)),
        KL_total_vs_sum = mean(kl_total_vs_sum, na.rm = TRUE)
      )
  )
)

# ── PIT data partitions ────────────────────────────────────────────────────────

# pooled overall
pit_states <- all_uni %>% filter(series != "Total") %>% pull(pit)
pit_total  <- all_uni %>% filter(series == "Total") %>% pull(pit)
pit_all    <- all_uni %>% pull(pit)

# per-horizon, split by level
pit_h_states <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, series != "Total") %>% pull(pit))
pit_h_total  <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, series == "Total") %>% pull(pit))
pit_h_all    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h) %>% pull(pit))

# ── KS tests against Uniform(0,1) ─────────────────────────────────────────────

ks_states <- ks.test(jitter_pit(pit_states), "punif")
ks_total  <- ks.test(jitter_pit(pit_total),  "punif")
ks_all    <- ks.test(jitter_pit(pit_all),    "punif")

ks_h_states <- lapply(pit_h_states, function(p) ks.test(jitter_pit(p), "punif"))
ks_h_total  <- lapply(pit_h_total,  function(p) ks.test(jitter_pit(p), "punif"))
ks_h_all    <- lapply(pit_h_all,    function(p) ks.test(jitter_pit(p), "punif"))

format_ks <- function(ks) sprintf("D=%.3f, p=%.3f", ks$statistic, ks$p.value)

# ── PIT plot helper ────────────────────────────────────────────────────────────

pit_histogram <- function(pit_vals, title, ks_obj, n_bins = 10) {
  df            <- data.frame(pit = pit_vals)
  expected_line <- length(pit_vals) / n_bins
  ggplot(df, aes(x = pit)) +
    geom_histogram(breaks   = seq(0, 1, length.out = n_bins + 1),
                   fill     = "#4C72B0", colour = "white", alpha = 0.85) +
    geom_hline(yintercept = expected_line, linetype = "dashed",
               colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = 0.98, y = Inf,
             label = format_ks(ks_obj),
             hjust = 1, vjust = 1.4, size = 2.8, colour = "grey30") +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = title, x = "PIT", y = "Count") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"))
}

# ── Build and save PDF ─────────────────────────────────────────────────────────

pdf(file.path(dir_plots, "pit_histograms.pdf"), width = 11, height = 4)

# pages 1-12: one page per horizon, 3 panels (States | Total | Overall)
for (h in 1:horizon) {
  p_s <- pit_histogram(pit_h_states[[h]],
                       sprintf("h = %d  |  States",  h), ks_h_states[[h]])
  p_t <- pit_histogram(pit_h_total[[h]],
                       sprintf("h = %d  |  Total",   h), ks_h_total[[h]])
  p_a <- pit_histogram(pit_h_all[[h]],
                       sprintf("h = %d  |  Overall", h), ks_h_all[[h]])
  print(
    (p_s | p_t | p_a) +
      plot_annotation(
        title = sprintf("PIT Histograms - Horizon h = %d", h),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )
}

# page 13: pooled overall (States | Total | Overall)
p_states_all <- pit_histogram(pit_states, "States (all horizons & origins)", ks_states)
p_total_all  <- pit_histogram(pit_total,  "Total (all horizons & origins)",  ks_total)
p_all_all    <- pit_histogram(pit_all,    "Overall (all series, horizons & origins)", ks_all)

print(
  (p_states_all | p_total_all | p_all_all) +
    plot_annotation(
      title = "PIT Histograms - Pooled Overall",
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
)

dev.off()

message("PIT plots saved to ", file.path(dir_plots, "pit_histograms.pdf"))

# ── Save aggregated results ────────────────────────────────────────────────────

save(
  all_uni, all_multi,
  summary_by_series_horizon, summary_by_horizon, summary_overall,
  file = file.path(dir_accuracy, "accuracy_es_pit.RData")
)
