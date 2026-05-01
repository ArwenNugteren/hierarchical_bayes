# ══════════════════════════════════════════════════════════════════════════════
# MCB Test: Friedman + Nemenyi post-hoc comparison of forecasting methods
# Hierarchy: Total → States → Zones
# Metrics: MAE, RMSE, CRPS, Coverage (distance from 90%)
# Granularities:
#   (1) Overall      — one observation per origin (pooled over horizons & series)
#   (2) By series    — one observation per origin × series
#   (3) Full         — one observation per origin × horizon × series
#   Each granularity is also run separately for three hierarchy levels:
#     (a) All series   (b) Total + States only   (c) Zones only
# Output: one PDF per metric (3 CD diagrams per granularity), + RData + CSV
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(stringr)
library(tsutils)     # nemenyi()

# ── 0. Configuration ───────────────────────────────────────────────────────────
bayesian_runs <- list(
  "Bayes default" = "accuracy_zones",
  "Bayes total_lambda=2" = "accuracy_zones_top_lambda2",
  "Bayes state_lambda=2" = "accuracy_zones_state_lambda2"
)

# Path to fable RData file for the zone hierarchy.
# Must contain: all_base_df and all_mint_df with columns:
#   origin, horizon, series, actual, forecast, abs_error, sq_error, crps, covered
fable_rdata_path <- "fable_zone_hierarchy_121_origins.RData"

# Output directory
output_dir <- "mcb_results_zones"
dir.create(output_dir, showWarnings = FALSE)


# ── 1. Load Bayesian results ───────────────────────────────────────────────────

load_bayesian <- function(folder, method_name) {
  files <- list.files(folder,
                      pattern    = "^accuracy_timestep_[0-9]+\\.rds$",
                      full.names = TRUE)
  if (length(files) == 0)
    stop(sprintf("No accuracy_timestep_*.rds files found in '%s'", folder))
  
  bind_rows(lapply(files, readRDS)) %>%
    mutate(method = method_name)
}

bayes_list <- mapply(load_bayesian,
                     folder      = bayesian_runs,
                     method_name = names(bayesian_runs),
                     SIMPLIFY    = FALSE)

bayes_df <- bind_rows(bayes_list)
# Expected columns: origin, horizon, series, actual, fc_mean,
#                   ae, se, crps, covered, method


# ── 2. Load fable results ──────────────────────────────────────────────────────

load(fable_rdata_path)   # loads all_base_df and all_mint_df

prep_fable <- function(df, method_name) {
  df %>%
    transmute(
      origin  = origin,
      horizon = horizon,
      series  = series,
      actual  = actual,
      fc_mean = forecast,
      ae      = abs_error,
      se      = sq_error,
      crps    = crps,
      covered = covered,
      method  = method_name
    )
}

fable_base_df <- prep_fable(all_base_df, "ETS (base)")
fable_mint_df <- prep_fable(all_mint_df, "ETS (MinT)")


# ── 3. Combine all methods ─────────────────────────────────────────────────────

all_df <- bind_rows(bayes_df, fable_base_df, fable_mint_df) %>%
  mutate(
    cov_dist = abs(covered - 0.90),
    # Classify each series by hierarchy level
    level = case_when(
      series == "Total"                          ~ "Total/State",
      str_detect(series, "^[A-G]$")             ~ "Total/State",
      str_detect(series, "^[A-G][A-Z]$")        ~ "Zone",
      TRUE                                       ~ "Other"
    )
  )

methods <- unique(all_df$method)
cat("Methods found:\n")
print(methods)

cat("\nSeries levels found:\n")
print(table(all_df$level))


# ── 4. Helper: build wide matrix for nemenyi() ────────────────────────────────

build_rank_matrix <- function(long_df, metric_col) {
  long_df %>%
    select(uid, method, value = all_of(metric_col)) %>%
    pivot_wider(names_from = method, values_from = value) %>%
    arrange(uid) %>%
    select(-uid) %>%
    as.matrix()
}


# ── 5. Build aggregated data frames at each granularity ───────────────────────
# Each function accepts a (possibly pre-filtered) all_df subset.

make_overall_df <- function(df) {
  df %>%
    group_by(origin, method) %>%
    summarise(
      MAE     = mean(ae,       na.rm = TRUE),
      RMSE    = sqrt(mean(se,  na.rm = TRUE)),
      CRPS    = mean(crps,     na.rm = TRUE),
      CovDist = mean(cov_dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(uid = as.character(origin))
}

make_series_df <- function(df) {
  df %>%
    group_by(origin, series, method) %>%
    summarise(
      MAE     = mean(ae,       na.rm = TRUE),
      RMSE    = sqrt(mean(se,  na.rm = TRUE)),
      CRPS    = mean(crps,     na.rm = TRUE),
      CovDist = mean(cov_dist, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(uid = paste(origin, series, sep = "_"))
}

make_full_df <- function(df) {
  df %>%
    mutate(
      MAE     = ae,
      RMSE    = sqrt(se),
      CRPS    = crps,
      CovDist = cov_dist,
      uid     = paste(origin, horizon, series, sep = "_")
    )
}


# ── 6. Define series subsets ───────────────────────────────────────────────────
# "all"          — every series
# "total_state"  — Total + state-level only  (level == "Total/State")
# "zone"         — zone-level only           (level == "Zone")

series_subsets <- list(
  list(label = "All series",         filter_expr = quote(level %in% c("Total/State", "Zone"))),
  list(label = "Total & states only",filter_expr = quote(level == "Total/State")),
  list(label = "Zones only",         filter_expr = quote(level == "Zone"))
)


# ── 7. Run MCB (Friedman + Nemenyi) ───────────────────────────────────────────

run_mcb <- function(long_df, metric_col, granularity_label) {
  
  mat <- build_rank_matrix(long_df, metric_col)
  
  complete_rows <- complete.cases(mat)
  n_dropped <- nrow(mat) - sum(complete_rows)
  if (n_dropped > 0) {
    warning(sprintf(
      "[%s / %s] Dropping %d incomplete rows (missing method values).",
      granularity_label, metric_col, n_dropped
    ))
    mat <- mat[complete_rows, , drop = FALSE]
  }
  
  if (nrow(mat) < 3) {
    warning(sprintf("[%s / %s] Too few observations (%d) — skipping.",
                    granularity_label, metric_col, nrow(mat)))
    return(NULL)
  }
  
  result <- tryCatch(
    tsutils::nemenyi(mat, conf.level = 0.95, plottype = "none"),
    error = function(e) {
      warning(sprintf("[%s / %s] nemenyi() failed: %s",
                      granularity_label, metric_col, conditionMessage(e)))
      NULL
    }
  )
  
  result
}


# ── 8. Main loop: metric → series subset → granularity ────────────────────────

metrics       <- c("MAE", "RMSE", "CRPS", "CovDist")
metric_labels <- c(
  MAE     = "MAE",
  RMSE    = "RMSE",
  CRPS    = "CRPS",
  CovDist = "Coverage distance from 90%"
)

granularity_makers <- list(
  list(label = "Overall (pooled horizons & series)", fn = make_overall_df),
  list(label = "By series (pooled horizons)",        fn = make_series_df),
  list(label = "Full (origin × horizon × series)",   fn = make_full_df)
)

# Storage
all_mcb_results <- list()
summary_tables  <- list()

for (metric in metrics) {
  
  # One PDF per metric; rows = series subsets, cols = granularities
  n_rows <- length(series_subsets)
  n_cols <- length(granularity_makers)
  
  pdf_path <- file.path(output_dir, sprintf("CD_%s.pdf", metric))
  pdf(pdf_path, width = 5 * n_cols, height = 5 * n_rows)
  par(mfrow = c(n_rows, n_cols))
  
  for (ss in series_subsets) {
    
    # Filter all_df to this series subset
    subset_df <- all_df %>% filter(eval(ss$filter_expr))
    
    for (gran in granularity_makers) {
      
      # Build aggregated df for this granularity
      agg_df <- gran$fn(subset_df)
      
      key <- paste(metric, ss$label, gran$label, sep = " | ")
      cat(sprintf("\n── MCB: %s | %s | %s ──\n", metric, ss$label, gran$label))
      
      res <- run_mcb(
        long_df           = agg_df,
        metric_col        = metric,
        granularity_label = paste(ss$label, gran$label, sep = " / ")
      )
      
      all_mcb_results[[key]] <- res
      
      plot_title <- sprintf("%s\n%s\n%s",
                            metric_labels[metric], ss$label, gran$label)
      
      if (!is.null(res)) {
        
        # Re-build clean matrix for plotting
        mat_plot <- build_rank_matrix(agg_df, metric)
        mat_plot <- mat_plot[complete.cases(mat_plot), , drop = FALSE]
        
        tsutils::nemenyi(
          mat_plot,
          conf.level = 0.95,
          plottype   = "vmcb",
          main       = plot_title
        )
        
        # Ranks table
        mean_ranks <- res$means
        tbl <- data.frame(
          metric       = metric,
          series_scope = ss$label,
          granularity  = gran$label,
          method       = names(mean_ranks),
          mean_rank    = as.numeric(mean_ranks),
          stringsAsFactors = FALSE
        ) %>%
          arrange(mean_rank)
        
        # Pairwise p-values
        pvals <- res$p.value
        if (!is.null(pvals)) {
          pval_long <- as.data.frame(as.table(pvals)) %>%
            setNames(c("method_1", "method_2", "p_value")) %>%
            mutate(
              metric       = metric,
              series_scope = ss$label,
              granularity  = gran$label
            ) %>%
            filter(!is.na(p_value))
        } else {
          pval_long <- NULL
        }
        
        summary_tables[[key]] <- list(ranks = tbl, pvalues = pval_long)
        
        cat("  Average ranks:\n");  print(tbl)
        if (!is.null(pval_long)) { cat("  Pairwise p-values:\n"); print(pval_long) }
        
      } else {
        plot.new()
        title(main = paste0(plot_title, "\n[FAILED]"))
      }
    }
  }
  
  dev.off()
  cat(sprintf("\nSaved: %s\n", pdf_path))
}


# ── 9. Save CSV summaries ──────────────────────────────────────────────────────

all_ranks <- bind_rows(lapply(summary_tables, function(x) x$ranks))
all_pvals <- bind_rows(lapply(summary_tables, function(x) x$pvalues))

write.csv(all_ranks, file.path(output_dir, "mcb_mean_ranks.csv"),    row.names = FALSE)
write.csv(all_pvals, file.path(output_dir, "mcb_pairwise_pvals.csv"), row.names = FALSE)

cat("\nSaved CSVs:\n")
cat(" ", file.path(output_dir, "mcb_mean_ranks.csv"), "\n")
cat(" ", file.path(output_dir, "mcb_pairwise_pvals.csv"), "\n")


# ── 10. Save full results as RData ─────────────────────────────────────────────

save(all_mcb_results, all_ranks, all_pvals, all_df,
     file = file.path(output_dir, "mcb_results_zones.RData"))

cat("\nSaved RData:", file.path(output_dir, "mcb_results_zones.RData"), "\n")
cat("\nDone.\n")