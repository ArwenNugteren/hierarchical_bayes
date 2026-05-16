# ── Aggregate results ──────────────────────────────────────────────────────────
library(dplyr)
files <- list.files(file.path("accuracy_states_independent_coherency"), 
                    pattern    = "^accuracy_timestep_[0-9]+\\.rds$", 
                    full.names = TRUE)
all_results_df <- bind_rows(lapply(files, readRDS))

# ── Coherence error (per timestep x horizon) ──────────────────────────────────
coherence_df <- all_results_df %>%
  group_by(origin, horizon) %>%
  summarise(
    coherence_error = fc_mean[series == "Total"] - 
      sum(fc_mean[series != "Total"]),
    .groups = "drop"
  )

# ── By series and horizon (coherence joined by horizon) ───────────────────────
coherence_by_horizon <- coherence_df %>%
  group_by(horizon) %>%
  summarise(
    MAE_coherence  = mean(abs(coherence_error)),
    RMSE_coherence = sqrt(mean(coherence_error^2)),
    .groups = "drop"
  )

summary_by_series_horizon <- all_results_df %>%
  group_by(series, horizon) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    Coverage = mean(covered),
    n        = n(),
    .groups  = "drop"
  ) %>%
  left_join(coherence_by_horizon, by = "horizon")
# NOTE: coherence is the same value for all series within a horizon —
# it reflects the aggregate mismatch, not a per-series contribution

# ── By horizon only ───────────────────────────────────────────────────────────
summary_by_horizon <- all_results_df %>%
  group_by(horizon) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    Coverage = mean(covered),
    n        = n(),
    .groups  = "drop"
  ) %>%
  left_join(coherence_by_horizon, by = "horizon")

# ── By series only + Overall ──────────────────────────────────────────────────
coherence_overall <- coherence_df %>%
  summarise(
    MAE_coherence  = mean(abs(coherence_error)),
    RMSE_coherence = sqrt(mean(coherence_error^2))
  )

summary_overall <- all_results_df %>%
  group_by(series) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    Coverage = mean(covered),
    .groups  = "drop"
  ) %>%
  bind_rows(
    all_results_df %>%
      summarise(
        series   = "Overall",
        MAE      = mean(ae),
        RMSE     = sqrt(mean(se)),
        CRPS     = mean(crps),
        Coverage = mean(covered)
      )
  ) %>%
  # Coherence columns are NA for individual series, populated only for Overall
  mutate(
    MAE_coherence  = if_else(series == "Overall", coherence_overall$MAE_coherence,  NA_real_),
    RMSE_coherence = if_else(series == "Overall", coherence_overall$RMSE_coherence, NA_real_)
  )

save(summary_by_series_horizon, summary_by_horizon, summary_overall, all_results_df,
     file = "accuracy_state_base.RData")

print(summary_overall)
# print(summary_by_horizon)
# print(summary_by_series_horizon)
