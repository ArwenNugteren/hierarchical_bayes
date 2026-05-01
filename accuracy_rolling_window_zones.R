# ── Aggregate results ──────────────────────────────────────────────────────────
library(dplyr)

files <- list.files(file.path("accuracy_zones_top_lambda2"), 
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

# ── Coherence (Total vs sum of bottom-level series) ───────────────────────────
# "Bottom" here = Zone level (nchar == 2); adjust if your bottom level differs

compute_coherence <- function(df) {
  df %>%
    group_by(origin, horizon) %>%
    summarise(
      coherence_error = fc_mean[series == "Total"] -
        sum(fc_mean[level == "Zone"]),
      .groups = "drop"
    )
}

coherence_by_horizon <- function(coh_df) {
  coh_df %>%
    group_by(horizon) %>%
    summarise(
      MAE_coherence  = mean(abs(coherence_error)),
      RMSE_coherence = sqrt(mean(coherence_error^2)),
      .groups = "drop"
    )
}

coherence_overall <- function(coh_df) {
  coh_df %>%
    summarise(
      MAE_coherence  = mean(abs(coherence_error)),
      RMSE_coherence = sqrt(mean(coherence_error^2))
    )
}

coh_df      <- compute_coherence(all_results_df)
coh_by_h    <- coherence_by_horizon(coh_df)
coh_overall <- coherence_overall(coh_df)

# ── Join coherence into summaries ──────────────────────────────────────────────

# By series and horizon: coherence joins by horizon, repeats across series/level
summary_by_series_horizon <- all_results_df %>%
  group_by(level, series, horizon) %>%
  summarise_metrics() %>%
  left_join(coh_by_h, by = "horizon")

# By horizon
summary_by_horizon_all <- all_results_df %>%
  group_by(horizon) %>%
  summarise_metrics() %>%
  mutate(level = "All", series = "Overall") %>%
  left_join(coh_by_h, by = "horizon")

summary_by_horizon_level <- all_results_df %>%
  group_by(level, horizon) %>%
  summarise_metrics() %>%
  mutate(series = "Overall (level)") %>%
  left_join(coh_by_h, by = "horizon")
# coherence is level-agnostic (Total vs bottom), same value repeated per level

summary_by_horizon <- bind_rows(summary_by_horizon_level, summary_by_horizon_all)

# Overall: coherence only meaningful at the aggregate "All" row
summary_by_series <- all_results_df %>%
  group_by(level, series) %>%
  summarise_metrics() %>%
  mutate(MAE_coherence  = NA_real_,
         RMSE_coherence = NA_real_)

summary_overall_level <- all_results_df %>%
  group_by(level) %>%
  summarise_metrics() %>%
  mutate(series         = "Overall (level)",
         MAE_coherence  = NA_real_,
         RMSE_coherence = NA_real_)

summary_overall_all <- all_results_df %>%
  summarise_metrics() %>%
  mutate(level          = "All",
         series         = "Overall",
         MAE_coherence  = coh_overall$MAE_coherence,
         RMSE_coherence = coh_overall$RMSE_coherence)

summary_overall <- bind_rows(summary_by_series, summary_overall_level, summary_overall_all)

#### multilevel coherence
# ── Multi-level coherence (Total → State → Zone) ──────────────────────────────
compute_coherence_multilevel <- function(df) {
  df %>%
    group_by(origin, horizon) %>%
    summarise(
      # Top vs State: Total should equal sum of all States
      coherence_error_total_state = fc_mean[series == "Total"] -
        sum(fc_mean[level == "State"]),
      # State vs Zone: sum of States should equal sum of all Zones
      coherence_error_state_zone  = sum(fc_mean[level == "State"]) -
        sum(fc_mean[level == "Zone"]),
      # Top vs Bottom: Total should equal sum of all Zones (derived, but useful)
      coherence_error_total_zone  = fc_mean[series == "Total"] -
        sum(fc_mean[level == "Zone"]),
      .groups = "drop"
    )
}

summarise_coherence_multilevel <- function(coh_df, group_var = NULL) {
  grp <- if (is.null(group_var)) coh_df else group_by(coh_df, across(all_of(group_var)))
  grp %>%
    summarise(
      MAE_coh_total_state  = mean(abs(coherence_error_total_state)),
      RMSE_coh_total_state = sqrt(mean(coherence_error_total_state^2)),
      MAE_coh_state_zone   = mean(abs(coherence_error_state_zone)),
      RMSE_coh_state_zone  = sqrt(mean(coherence_error_state_zone^2)),
      MAE_coh_total_zone   = mean(abs(coherence_error_total_zone)),
      RMSE_coh_total_zone  = sqrt(mean(coherence_error_total_zone^2)),
      .groups = "drop"
    )
}

coh_ml_df      <- compute_coherence_multilevel(all_results_df)
coh_ml_by_h    <- summarise_coherence_multilevel(coh_ml_df, group_var = "horizon")
coh_ml_overall <- summarise_coherence_multilevel(coh_ml_df)

# ── Save ───────────────────────────────────────────────────────────────────────
save(summary_by_series_horizon, summary_by_horizon, summary_by_series,
     summary_overall, all_results_df,
     file = "zone_hierarchy_accuracy.RData")

print(summary_overall, n=39)
print(coh_ml_overall)
