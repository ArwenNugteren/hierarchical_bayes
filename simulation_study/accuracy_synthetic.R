library(dplyr)
library(ggplot2)
library(patchwork)

# ── User-set paths ─────────────────────────────────────────────────────────────

dir_accuracy <- file.path("accuracy_wellspecified_joint_coherency0.8")
dir_plots    <- file.path(dir_accuracy, "plots")

dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)

horizon <- 5L

# ── Load all per-roll accuracy files ──────────────────────────────────────────

acc_files <- list.files(dir_accuracy,
                        pattern    = "^accuracy_roll_[0-9]+\\.rds$",
                        full.names = TRUE)

all_uni   <- bind_rows(lapply(acc_files, function(f) readRDS(f)$univariate))
all_multi <- bind_rows(lapply(acc_files, function(f) readRDS(f)$multivariate))

# ── Summary by horizon ────────────────────────────────────────────────────────

summary_by_horizon <- all_uni %>%
  group_by(horizon) %>%
  summarise(
    MAE   = mean(ae),
    RMSE  = sqrt(mean(se)),
    CRPS  = mean(crps),
    LogS  = mean(logS[is.finite(logS)]),
    Cov90 = mean(covered_90),
    Cov95 = mean(covered_95),
    .groups = "drop"
  ) %>%
  left_join(
    all_multi %>%
      group_by(horizon) %>%
      summarise(
        ES          = mean(es),
        MAI_bot_mid = mean(abs(mai_bot_mid)),
        MAI_mid_top = mean(abs(mai_mid_top)),
        MAI_bot_top = mean(abs(mai_bot_top)),
        JSD_bot_mid = mean(jsd_bot_mid, na.rm = TRUE),
        JSD_mid_top = mean(jsd_mid_top, na.rm = TRUE),
        JSD_bot_top = mean(jsd_bot_top, na.rm = TRUE),
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

multi_overall <- all_multi %>%
  summarise(
    ES          = mean(es),
    MAI_bot_mid = mean(abs(mai_bot_mid)),
    MAI_mid_top = mean(abs(mai_mid_top)),
    MAI_bot_top = mean(abs(mai_bot_top)),
    JSD_bot_mid = mean(jsd_bot_mid, na.rm = TRUE),
    JSD_mid_top = mean(jsd_mid_top, na.rm = TRUE),
    JSD_bot_top = mean(jsd_bot_top, na.rm = TRUE)
  )

make_level_row <- function(lev) {
  all_uni %>%
    filter(level == lev) %>%
    summarise(
      level       = lev,
      MAE         = mean(ae),
      RMSE        = sqrt(mean(se)),
      CRPS        = mean(crps),
      LogS        = mean(logS[is.finite(logS)]),
      Cov90       = mean(covered_90),
      Cov95       = mean(covered_95),
      ES          = NA_real_,
      MAI_bot_mid = NA_real_,
      MAI_mid_top = NA_real_,
      MAI_bot_top = NA_real_,
      JSD_bot_mid = NA_real_,
      JSD_mid_top = NA_real_,
      JSD_bot_top = NA_real_
    )
}

summary_overall <- bind_rows(
  make_level_row("bottom"),
  make_level_row("mid"),
  make_level_row("top"),
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

pit_bottom <- all_uni %>% filter(level == "bottom") %>% pull(pit)
pit_mid    <- all_uni %>% filter(level == "mid")    %>% pull(pit)
pit_top    <- all_uni %>% filter(level == "top")    %>% pull(pit)
pit_all    <- all_uni %>% pull(pit)

pit_h_bottom <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "bottom") %>% pull(pit))
pit_h_mid    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "mid")    %>% pull(pit))
pit_h_top    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "top")    %>% pull(pit))
pit_h_all    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h) %>% pull(pit))

# ── KS tests ──────────────────────────────────────────────────────────────────

ks_bottom <- ks.test(jitter_pit(pit_bottom), "punif")
ks_mid    <- ks.test(jitter_pit(pit_mid),    "punif")
ks_top    <- ks.test(jitter_pit(pit_top),    "punif")
ks_all    <- ks.test(jitter_pit(pit_all),    "punif")

ks_h_bottom <- lapply(pit_h_bottom, function(p) ks.test(jitter_pit(p), "punif"))
ks_h_mid    <- lapply(pit_h_mid,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_top    <- lapply(pit_h_top,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_all    <- lapply(pit_h_all,    function(p) ks.test(jitter_pit(p), "punif"))

# ── Build and save PDF ─────────────────────────────────────────────────────────

pdf(file.path(dir_plots, "pit_histograms.pdf"), width = 14, height = 4)

for (h in 1:horizon) {
  p_bot <- pit_histogram(pit_h_bottom[[h]],
                         sprintf("h=%d | Bottom (AA/AB/BA/BB)", h), ks_h_bottom[[h]])
  p_mid <- pit_histogram(pit_h_mid[[h]],
                         sprintf("h=%d | Mid (A/B)",            h), ks_h_mid[[h]])
  p_top <- pit_histogram(pit_h_top[[h]],
                         sprintf("h=%d | Top (Total)",          h), ks_h_top[[h]])
  p_all <- pit_histogram(pit_h_all[[h]],
                         sprintf("h=%d | Overall",              h), ks_h_all[[h]])
  print(
    (p_bot | p_mid | p_top | p_all) +
      plot_annotation(
        title = sprintf("PIT Histograms - Horizon h = %d", h),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )
}

p_bot_all <- pit_histogram(pit_bottom, "Bottom (all horizons & rolls)", ks_bottom)
p_mid_all <- pit_histogram(pit_mid,    "Mid    (all horizons & rolls)", ks_mid)
p_top_all <- pit_histogram(pit_top,    "Top    (all horizons & rolls)", ks_top)
p_all_all <- pit_histogram(pit_all,    "Overall (all levels, horizons & rolls)", ks_all)

print(
  (p_bot_all | p_mid_all | p_top_all | p_all_all) +
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
  file = file.path(dir_accuracy, "accuracy_synthetic.RData")
)

# ── Print ──────────────────────────────────────────────────────────────────────

reduced_summary <- summary_overall %>%
  select(-MAI_bot_mid, -MAI_mid_top, -JSD_bot_mid, -JSD_mid_top) %>%
  mutate(
    across(c(MAE, RMSE, CRPS, LogS, ES, MAI_bot_top), round, 3),
    across(c(Cov90, Cov95, JSD_bot_top),               round, 3)
  )
# reduced_summary <- summary_overall %>%
#   select(MAE, RMSE, CRPS, LogS, ES, MAI_bot_top, Cov90, Cov95, JSD_bot_top)

print(reduced_summary)

write.csv(reduced_summary,
          file = file.path(dir_accuracy, "reduced_summary.csv"),
          row.names = FALSE)