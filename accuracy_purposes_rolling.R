library(dplyr)
library(ggplot2)
library(patchwork)

# ── User-set paths ─────────────────────────────────────────────────────────────

dir_accuracy <- file.path("accuracy_purposes_independentv2")
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
        ES               = mean(es),
        MAI_rp_zp        = mean(abs(mai_rp_zp)),
        MAI_zp_sp        = mean(abs(mai_zp_sp)),
        MAI_sp_pp        = mean(abs(mai_sp_pp)),
        MAI_pp_total     = mean(abs(mai_pp_total)),
        MAI_rp_total     = mean(abs(mai_rp_total)),
        JSD_rp_zp        = mean(jsd_rp_zp,    na.rm = TRUE),
        JSD_zp_sp        = mean(jsd_zp_sp,    na.rm = TRUE),
        JSD_sp_pp        = mean(jsd_sp_pp,    na.rm = TRUE),
        JSD_pp_total     = mean(jsd_pp_total, na.rm = TRUE),
        JSD_rp_total     = mean(jsd_rp_total, na.rm = TRUE),
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
    ES               = mean(es),
    MAI_rp_zp        = mean(abs(mai_rp_zp)),
    MAI_zp_sp        = mean(abs(mai_zp_sp)),
    MAI_sp_pp        = mean(abs(mai_sp_pp)),
    MAI_pp_total     = mean(abs(mai_pp_total)),
    MAI_rp_total     = mean(abs(mai_rp_total)),
    JSD_rp_zp        = mean(jsd_rp_zp,    na.rm = TRUE),
    JSD_zp_sp        = mean(jsd_zp_sp,    na.rm = TRUE),
    JSD_sp_pp        = mean(jsd_sp_pp,    na.rm = TRUE),
    JSD_pp_total     = mean(jsd_pp_total, na.rm = TRUE),
    JSD_rp_total     = mean(jsd_rp_total, na.rm = TRUE)
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
      MAI_rp_zp        = NA_real_,
      MAI_zp_sp        = NA_real_,
      MAI_sp_pp        = NA_real_,
      MAI_pp_total     = NA_real_,
      MAI_rp_total     = NA_real_,
      JSD_rp_zp        = NA_real_,
      JSD_zp_sp        = NA_real_,
      JSD_sp_pp        = NA_real_,
      JSD_pp_total     = NA_real_,
      JSD_rp_total     = NA_real_
    )
}

summary_overall <- bind_rows(
  make_level_row("RegionPurpose"),
  make_level_row("ZonePurpose"),
  make_level_row("StatePurpose"),
  make_level_row("Purpose"),
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

pit_rp      <- all_uni %>% filter(level == "RegionPurpose") %>% pull(pit)
pit_zp      <- all_uni %>% filter(level == "ZonePurpose")   %>% pull(pit)
pit_sp      <- all_uni %>% filter(level == "StatePurpose")  %>% pull(pit)
pit_pp      <- all_uni %>% filter(level == "Purpose")       %>% pull(pit)
pit_total   <- all_uni %>% filter(level == "Total")         %>% pull(pit)
pit_all     <- all_uni %>% pull(pit)

pit_h_rp    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "RegionPurpose") %>% pull(pit))
pit_h_zp    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "ZonePurpose")   %>% pull(pit))
pit_h_sp    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "StatePurpose")  %>% pull(pit))
pit_h_pp    <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "Purpose")       %>% pull(pit))
pit_h_total <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h, level == "Total")         %>% pull(pit))
pit_h_all   <- lapply(1:horizon, function(h)
  all_uni %>% filter(horizon == h) %>% pull(pit))

# ── KS tests ──────────────────────────────────────────────────────────────────

ks_rp    <- ks.test(jitter_pit(pit_rp),    "punif")
ks_zp    <- ks.test(jitter_pit(pit_zp),    "punif")
ks_sp    <- ks.test(jitter_pit(pit_sp),    "punif")
ks_pp    <- ks.test(jitter_pit(pit_pp),    "punif")
ks_total <- ks.test(jitter_pit(pit_total), "punif")
ks_all   <- ks.test(jitter_pit(pit_all),   "punif")

ks_h_rp    <- lapply(pit_h_rp,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_zp    <- lapply(pit_h_zp,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_sp    <- lapply(pit_h_sp,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_pp    <- lapply(pit_h_pp,    function(p) ks.test(jitter_pit(p), "punif"))
ks_h_total <- lapply(pit_h_total, function(p) ks.test(jitter_pit(p), "punif"))
ks_h_all   <- lapply(pit_h_all,   function(p) ks.test(jitter_pit(p), "punif"))

# ── Build and save PDF ─────────────────────────────────────────────────────────

pdf(file.path(dir_plots, "pit_histograms.pdf"), width = 21, height = 4)

# pages 1-12: one page per horizon, 6 panels
for (h in 1:horizon) {
  p_rp <- pit_histogram(pit_h_rp[[h]],
                        sprintf("h=%d | Region×Purpose", h), ks_h_rp[[h]])
  p_zp <- pit_histogram(pit_h_zp[[h]],
                        sprintf("h=%d | Zone×Purpose",   h), ks_h_zp[[h]])
  p_sp <- pit_histogram(pit_h_sp[[h]],
                        sprintf("h=%d | State×Purpose",  h), ks_h_sp[[h]])
  p_pp <- pit_histogram(pit_h_pp[[h]],
                        sprintf("h=%d | Purpose",        h), ks_h_pp[[h]])
  p_t  <- pit_histogram(pit_h_total[[h]],
                        sprintf("h=%d | Total",          h), ks_h_total[[h]])
  p_a  <- pit_histogram(pit_h_all[[h]],
                        sprintf("h=%d | Overall",        h), ks_h_all[[h]])
  print(
    (p_rp | p_zp | p_sp | p_pp | p_t | p_a) +
      plot_annotation(
        title = sprintf("PIT Histograms - Horizon h = %d", h),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )
}

# final page: pooled overall, 6 panels
p_rp_all <- pit_histogram(pit_rp,    "Region×Purpose (all horizons & origins)", ks_rp)
p_zp_all <- pit_histogram(pit_zp,    "Zone×Purpose   (all horizons & origins)", ks_zp)
p_sp_all <- pit_histogram(pit_sp,    "State×Purpose  (all horizons & origins)", ks_sp)
p_pp_all <- pit_histogram(pit_pp,    "Purpose        (all horizons & origins)", ks_pp)
p_t_all  <- pit_histogram(pit_total, "Total          (all horizons & origins)", ks_total)
p_a_all  <- pit_histogram(pit_all,   "Overall (all levels, horizons & origins)", ks_all)

print(
  (p_rp_all | p_zp_all | p_sp_all | p_pp_all | p_t_all | p_a_all) +
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
  file = file.path(dir_accuracy, "accuracy_purpose.RData")
)

# ── Print ──────────────────────────────────────────────────────────────────────

reduced_summary <- summary_overall %>%
  select(-MAI_rp_zp, -MAI_zp_sp, -MAI_sp_pp, -MAI_pp_total,
         -JSD_rp_zp, -JSD_zp_sp, -JSD_sp_pp, -JSD_pp_total) %>%
  mutate(
    across(c(MAE, RMSE, CRPS, LogS, ES, MAI_rp_total), round, 1),
    across(c(Cov90, Cov95, JSD_rp_total),               round, 3)
  )

print(reduced_summary)

write.csv(reduced_summary, file = file.path(dir_accuracy, "reduced_summary.csv"))