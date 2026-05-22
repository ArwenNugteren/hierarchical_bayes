library(dplyr)
library(ggplot2)
library(patchwork)

# ── User-set paths ─────────────────────────────────────────────────────────────

dir_accuracy <- file.path("accuracy_regions_hierarchy")
dir_plots    <- file.path(dir_accuracy, "plots")

dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

horizon <- 12L

# ── Load all per-origin accuracy files ────────────────────────────────────────

acc_files <- list.files(dir_accuracy,
                        pattern    = "^accuracy_timestep_[0-9]+\\.rds$",
                        full.names = TRUE)

all_uni   <- bind_rows(lapply(acc_files, function(f) readRDS(f)$univariate))
all_multi <- bind_rows(lapply(acc_files, function(f) readRDS(f)$multivariate))

# ── Summary by horizon ────────────────────────────────────────────────────────

summary_by_horizon <- all_uni %>%
  group_by(horizon) %>%
  summarise(
    MAE       = mean(ae),
    RMSE      = sqrt(mean(se)),
    CRPS      = mean(crps),
    LogS      = mean(logS[is.finite(logS)]),
    Cov90     = mean(covered_90),
    Cov95     = mean(covered_95),
    .groups = "drop"
  ) %>%
  left_join(
    all_multi %>%
      group_by(horizon) %>%
      summarise(
        ES               = mean(es),
        MAI_region_zone  = mean(abs(mai_region_zone)),
        MAI_zone_state   = mean(abs(mai_zone_state)),
        MAI_state_total  = mean(abs(mai_state_total)),
        MAI_region_total = mean(abs(mai_region_total)),
        JSD_region_zone  = mean(jsd_region_zone,  na.rm = TRUE),
        JSD_zone_state   = mean(jsd_zone_state,   na.rm = TRUE),
        JSD_state_total  = mean(jsd_state_total,  na.rm = TRUE),
        JSD_region_total = mean(jsd_region_total, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "horizon"
  )

# ── Summary by level and horizon ──────────────────────────────────────────────

summary_by_level_horizon <- all_uni %>%
  group_by(level, horizon) %>%
  summarise(
    MAE   = mean(ae),
    RMSE  = sqrt(mean(se)),
    CRPS  = mean(crps),
    LogS  = mean(logS[is.finite(logS)]),
    Cov90 = mean(covered_90),
    Cov95 = mean(covered_95),
    .groups = "drop"
  )

# ── Summary by series and horizon ─────────────────────────────────────────────

summary_by_series_horizon <- all_uni %>%
  group_by(level, series, horizon) %>%
  summarise(
    MAE   = mean(ae),
    RMSE  = sqrt(mean(se)),
    CRPS  = mean(crps),
    LogS  = mean(logS[is.finite(logS)]),
    Cov90 = mean(covered_90),
    Cov95 = mean(covered_95),
    .groups = "drop"
  )

# ── Summary overall ───────────────────────────────────────────────────────────
# Rows: Region / Zone / State / Total / Overall
# Joint metrics (ES, MAI, JSD) only populated at Overall level

multi_overall <- all_multi %>%
  summarise(
    ES               = mean(es),
    MAI_region_zone  = mean(abs(mai_region_zone)),
    MAI_zone_state   = mean(abs(mai_zone_state)),
    MAI_state_total  = mean(abs(mai_state_total)),
    MAI_region_total = mean(abs(mai_region_total)),
    JSD_region_zone  = mean(jsd_region_zone,  na.rm = TRUE),
    JSD_zone_state   = mean(jsd_zone_state,   na.rm = TRUE),
    JSD_state_total  = mean(jsd_state_total,  na.rm = TRUE),
    JSD_region_total = mean(jsd_region_total, na.rm = TRUE)
  )

make_level_row <- function(lev) {
  all_uni %>%
    filter(level == lev) %>%
    summarise(
      level            = lev,
      MAE              = mean(ae),
      RMSE             = sqrt(mean(se)),
      CRPS             = mean(crps),
      LogS             = mean(logS[is.finite(logS)]),
      Cov90            = mean(covered_90),
      Cov95            = mean(covered_95),
      ES               = NA_real_,
      MAI_region_zone  = NA_real_,
      MAI_zone_state   = NA_real_,
      MAI_state_total  = NA_real_,
      MAI_region_total = NA_real_,
      JSD_region_zone  = NA_real_,
      JSD_zone_state   = NA_real_,
      JSD_state_total  = NA_real_,
      JSD_region_total = NA_real_
    )
}

summary_overall <- bind_rows(
  make_level_row("Region"),
  make_level_row("Zone"),
  make_level_row("State"),
  make_level_row("Total"),
  bind_cols(
    all_uni %>%
      summarise(
        level = "Overall",
        MAE   = mean(ae),
        RMSE  = sqrt(mean(se)),
        CRPS  = mean(crps),
        LogS  = mean(logS[is.finite(logS)]),
        Cov90 = mean(covered_90),
        Cov95 = mean(covered_95)
      ),
    multi_overall
  )
)

# ── PIT helpers ───────────────────────────────────────────────────────────────

jitter_pit <- function(x, amount = 1e-4) {
  pmax(0, pmin(1, x + runif(length(x), -amount, amount)))
}

format_ks <- function(ks) sprintf("D=%.3f, p=%.3f", ks$statistic, ks$p.value)

pit_histogram <- function(pit_vals, title, ks_obj, n_bins = 10) {
  df            <- data.frame(pit = pit_vals)
  expected_line <- length(pit_vals) / n_bins
  ggplot(df, aes(x = pit)) +
    geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                   fill = "#4C72B0", colour = "white", alpha = 0.85) +
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

# ── PIT partitions ────────────────────────────────────────────────────────────

pit_region  <- all_uni %>% filter(level == "Region") %>% pull(pit)
pit_zone    <- all_uni %>% filter(level == "Zone")   %>% pull(pit)
pit_state   <- all_uni %>% filter(level == "State")  %>% pull(pit)
pit_total   <- all_uni %>% filter(level == "Total")  %>% pull(pit)
pit_all     <- all_uni %>% pull(pit)

pit_h_region <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "Region") %>% pull(pit))
pit_h_zone   <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "Zone")   %>% pull(pit))
pit_h_state  <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "State")  %>% pull(pit))
pit_h_total  <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "Total")  %>% pull(pit))
pit_h_all    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h) %>% pull(pit))

# ── KS tests ──────────────────────────────────────────────────────────────────

ks_region <- ks.test(jitter_pit(pit_region), "punif")
ks_zone   <- ks.test(jitter_pit(pit_zone),   "punif")
ks_state  <- ks.test(jitter_pit(pit_state),  "punif")
ks_total  <- ks.test(jitter_pit(pit_total),  "punif")
ks_all    <- ks.test(jitter_pit(pit_all),    "punif")

ks_h_region <- lapply(pit_h_region, function(p) ks.test(jitter_pit(p), "punif"))
ks_h_zone   <- lapply(pit_h_zone,   function(p) ks.test(jitter_pit(p), "punif"))
ks_h_state  <- lapply(pit_h_state,  function(p) ks.test(jitter_pit(p), "punif"))
ks_h_total  <- lapply(pit_h_total,  function(p) ks.test(jitter_pit(p), "punif"))
ks_h_all    <- lapply(pit_h_all,    function(p) ks.test(jitter_pit(p), "punif"))

# ── Build and save PDF ─────────────────────────────────────────────────────────

pdf(file.path(dir_plots, "pit_histograms.pdf"), width = 18, height = 4)

# pages 1-12: one page per horizon, 5 panels (Region | Zone | State | Total | Overall)
for (h in 1:horizon) {
  p_r <- pit_histogram(pit_h_region[[h]],
                       sprintf("h = %d  |  Region",  h), ks_h_region[[h]])
  p_z <- pit_histogram(pit_h_zone[[h]],
                       sprintf("h = %d  |  Zone",    h), ks_h_zone[[h]])
  p_s <- pit_histogram(pit_h_state[[h]],
                       sprintf("h = %d  |  State",   h), ks_h_state[[h]])
  p_t <- pit_histogram(pit_h_total[[h]],
                       sprintf("h = %d  |  Total",   h), ks_h_total[[h]])
  p_a <- pit_histogram(pit_h_all[[h]],
                       sprintf("h = %d  |  Overall", h), ks_h_all[[h]])
  print(
    (p_r | p_z | p_s | p_t | p_a) +
      plot_annotation(
        title = sprintf("PIT Histograms - Horizon h = %d", h),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )
}

# final page: pooled overall, 5 panels
p_region_all <- pit_histogram(pit_region, "Region (all horizons & origins)",          ks_region)
p_zone_all   <- pit_histogram(pit_zone,   "Zone   (all horizons & origins)",          ks_zone)
p_state_all  <- pit_histogram(pit_state,  "State  (all horizons & origins)",          ks_state)
p_total_all  <- pit_histogram(pit_total,  "Total  (all horizons & origins)",          ks_total)
p_all_all    <- pit_histogram(pit_all,    "Overall (all levels, horizons & origins)", ks_all)

print(
  (p_region_all | p_zone_all | p_state_all | p_total_all | p_all_all) +
    plot_annotation(
      title = "PIT Histograms - Pooled Overall",
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
)

dev.off()

message("PIT plots saved to ", file.path(dir_plots, "pit_histograms.pdf"))

# ── Save ──────────────────────────────────────────────────────────────────────

save(
  all_uni, all_multi,
  summary_overall,
  summary_by_horizon,
  summary_by_level_horizon,
  summary_by_series_horizon,
  file = file.path(dir_accuracy, "accuracy_zones.RData")
)

# ── Print ──────────────────────────────────────────────────────────────────────
reduced_summary <- summary_overall %>%
  select(-MAI_region_zone, -MAI_zone_state, -MAI_state_total,
         -JSD_region_zone, -JSD_zone_state, -JSD_state_total) %>%
  mutate(
    across(c(MAE, RMSE, CRPS, LogS, ES, MAI_region_total), round, 1),  # large-ish values
    across(c(Cov90, Cov95, JSD_region_total),               round, 3)  # proportions / small decimals
  )

print(reduced_summary)

write.csv(reduced_summary, file = file.path(dir_accuracy, "reduced_summary.csv"))
