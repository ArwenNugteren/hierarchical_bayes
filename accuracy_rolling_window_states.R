# ── Aggregate results ──────────────────────────────────────────────────────────
library(dplyr)
files <- list.files(file.path("accuracy"), 
                    pattern    = "^accuracy_timestep_[0-9]+\\.rds$", 
                    full.names = TRUE)
all_results_df <- bind_rows(lapply(files, readRDS))

# ── By series and horizon ──────────────────────────────────────────────────────
summary_by_series_horizon <- all_results_df %>%
  group_by(series, horizon) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    Coverage = mean(covered),
    n        = n(),
    .groups  = "drop"
  )

# ── By horizon only (pooled across all series) ────────────────────────────────
summary_by_horizon <- all_results_df %>%
  group_by(horizon) %>%
  summarise(
    MAE      = mean(ae),
    RMSE     = sqrt(mean(se)),
    CRPS     = mean(crps),
    Coverage = mean(covered),
    n        = n(),
    .groups  = "drop"
  )

# ── By series only (pooled across all horizons) ───────────────────────────────
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
  )

save(summary_by_series_horizon, summary_by_horizon, summary_overall, all_results_df,
     file = "state_hierarchy_full_t.RData")

print(summary_overall)
print(summary_by_horizon)
print(summary_by_series_horizon)
