# ── Aggregate results ──────────────────────────────────────────────────────────
library(dplyr)

files <- list.files(file.path("accuracy"), 
                    pattern    = "^accuracy_timestep_[0-9]+\\.rds$", 
                    full.names = TRUE)
all_results_df <- bind_rows(lapply(files, readRDS))

# ── Classify series by level ───────────────────────────────────────────────────

all_results_df <- all_results_df %>%
  mutate(level = case_when(
    series == "Total"                ~ "Total",
    nchar(series) == 1               ~ "State",
    nchar(series) == 2               ~ "Zone",
    TRUE                             ~ "Other"
  ))

# ── Summary functions ──────────────────────────────────────────────────────────

summarise_metrics <- function(df) {
  df %>%
    summarise(
      MAE      = mean(ae),
      RMSE     = sqrt(mean(se)),
      CRPS     = mean(crps),
      Coverage = mean(covered),
      .groups  = "drop"
    )
}

# ── By series and horizon ──────────────────────────────────────────────────────

summary_by_series_horizon <- all_results_df %>%
  group_by(level, series, horizon) %>%
  summarise_metrics()

# ── By horizon only — pooled across all series ────────────────────────────────

summary_by_horizon_all <- all_results_df %>%
  group_by(horizon) %>%
  summarise_metrics() %>%
  mutate(level = "All", series = "Overall")

# ── By horizon — pooled within each level ─────────────────────────────────────

summary_by_horizon_level <- all_results_df %>%
  group_by(level, horizon) %>%
  summarise_metrics() %>%
  mutate(series = "Overall (level)")

# ── By series only — pooled across all horizons ───────────────────────────────

summary_by_series <- all_results_df %>%
  group_by(level, series) %>%
  summarise_metrics()

# ── Overall — pooled across everything ────────────────────────────────────────

summary_overall_all <- all_results_df %>%
  summarise_metrics() %>%
  mutate(level = "All", series = "Overall")

# ── Overall — pooled within each level ────────────────────────────────────────

summary_overall_level <- all_results_df %>%
  group_by(level) %>%
  summarise_metrics() %>%
  mutate(series = "Overall (level)")

# ── Combined overall table ─────────────────────────────────────────────────────

summary_overall <- bind_rows(summary_by_series, summary_overall_level, summary_overall_all)

# ── Combined by-horizon table ──────────────────────────────────────────────────

summary_by_horizon <- bind_rows(summary_by_horizon_level, summary_by_horizon_all)

# ── Save ───────────────────────────────────────────────────────────────────────

save(summary_by_series_horizon, summary_by_horizon, summary_by_series,
     summary_overall, all_results_df,
     file = "zone_hierarchy_accuracy.RData")

print(summary_overall, n=39)
print(summary_by_horizon)
print(summary_by_series_horizon)