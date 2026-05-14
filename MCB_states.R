# ══════════════════════════════════════════════════════════════════════════════
# MCB Test: Friedman + Nemenyi post-hoc comparison of forecasting methods
# Metrics: MAE, RMSE, CRPS, Coverage (distance from 90%)
# Granularities:
#   (1) Overall      — one observation per origin (pooled over horizons & series)
#   (2) By series    — one observation per origin, per series
#   (3) Full         — one observation per origin × horizon × series
# Output: one PDF per metric (3 CD diagrams inside), + RData + CSV summaries
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(tsutils)     # nemenyi()
library(PMCMRplus)   # frdAllPairsNemenyiTest() — used internally for p-values

# ── 0. Configuration ───────────────────────────────────────────────────────────
# Each Bayesian run: named list entry mapping display name → folder path.
# The folder must contain files matching accuracy_timestep_NNN.rds as produced
# by the PBS array job.  Add or remove entries freely; the rest of the script
# adapts automatically.

bayesian_runs <- list(
  "Bayes unweighted" = "accuracy_states_normal",
  "Bayes coherency" = "accuracy_states_coherency_new"
  # "Bayes lambda=5" = "accuracy_states_lambda5"
)

# Path to the fable RData file (produced by the fable script).
# Must contain objects: all_base_df  and  all_mint_df
fable_rdata_path <- "fable_state_hierarchy_121_origins.RData"   # e.g. "fable_state_hierarchy_121_origins.RData"

# Output directory for PDFs, CSVs and RData
output_dir <- "mcb_results_coherence"
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
# Expected columns from Bayesian runs: origin, horizon, series, actual,
#   fc_mean, ae, se, crps, covered

# ── 2. Load fable results ──────────────────────────────────────────────────────

load(fable_rdata_path)   # loads all_base_df and all_mint_df

# Fable columns expected: origin, horizon, series, actual, forecast,
#   abs_error, sq_error, crps, covered

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
    cov_dist = abs(covered - 0.90)   # distance from nominal 90% coverage
  )

methods <- unique(all_df$method)
cat("Methods found:\n")
print(methods)


# ── 4. Helper: build wide matrix for nemenyi() ─────────────────────────────────
# nemenyi() expects a matrix with rows = datasets (e.g. origins) and
# columns = methods.  Ties in ranks are handled by the function.

build_rank_matrix <- function(long_df, metric_col) {
  long_df %>%
    select(uid, method, value = all_of(metric_col)) %>%
    pivot_wider(names_from = method, values_from = value) %>%
    arrange(uid) %>%
    select(-uid) %>%
    as.matrix()
}


# ── 5. Compute per-observation metrics at each granularity ────────────────────

# ── 5a. Overall: one row per origin (pooled over all horizons & series) ────────
overall_df <- all_df %>%
  group_by(origin, method) %>%
  summarise(
    MAE      = mean(ae,       na.rm = TRUE),
    RMSE     = sqrt(mean(se, na.rm = TRUE)),
    CRPS     = mean(crps,     na.rm = TRUE),
    CovDist  = mean(cov_dist, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(uid = as.character(origin))

# ── 5b. By series: one row per origin × series ────────────────────────────────
series_df <- all_df %>%
  group_by(origin, series, method) %>%
  summarise(
    MAE      = mean(ae,       na.rm = TRUE),
    RMSE     = sqrt(mean(se, na.rm = TRUE)),
    CRPS     = mean(crps,     na.rm = TRUE),
    CovDist  = mean(cov_dist, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(uid = paste(origin, series, sep = "_"))

# ── 5c. Full: one row per origin × horizon × series ───────────────────────────
full_df <- all_df %>%
  mutate(
    MAE     = ae,
    RMSE    = sqrt(se),
    CRPS    = crps,
    CovDist = cov_dist,
    uid     = paste(origin, horizon, series, sep = "_")
  )


# ── 6. Run MCB (Friedman + Nemenyi) for one granularity + metric ──────────────

run_mcb <- function(long_df, metric_col, granularity_label) {
  
  mat <- build_rank_matrix(long_df, metric_col)
  
  # Drop any rows with NA (can happen if a method is missing for some origins)
  complete_rows <- complete.cases(mat)
  if (sum(complete_rows) < nrow(mat)) {
    warning(sprintf(
      "[%s / %s] Dropping %d incomplete rows (missing method values).",
      granularity_label, metric_col,
      nrow(mat) - sum(complete_rows)
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


# ── 7. Collect all results and produce outputs ─────────────────────────────────

metrics        <- c("MAE", "RMSE", "CRPS", "CovDist")
metric_labels  <- c(MAE = "MAE", RMSE = "RMSE", CRPS = "CRPS",
                    CovDist = "Coverage distance from 90%")

granularities <- list(
  list(label = "Overall (pooled horizons & series)", df = overall_df),
  list(label = "By series (pooled horizons)",        df = series_df),
  list(label = "Full (origin × horizon × series)",   df = full_df)
)

# Storage for all results (for RData save)
all_mcb_results <- list()

# Summary table list (for CSV)
summary_tables <- list()

for (metric in metrics) {
  
  pdf_path <- file.path(output_dir, sprintf("CD_%s.pdf", metric))
  pdf(pdf_path, width = 10, height = 5)
  par(mfrow = c(1, 3))
  
  for (gran in granularities) {
    
    key <- paste(metric, gran$label, sep = " | ")
    cat(sprintf("\n── Running MCB: %s | %s ──\n", metric, gran$label))
    
    res <- run_mcb(long_df        = gran$df,
                   metric_col     = metric,
                   granularity_label = gran$label)
    
    all_mcb_results[[key]] <- res
    
    if (!is.null(res)) {
      
      # CD diagram
      tsutils::nemenyi(
        build_rank_matrix(gran$df, metric)[complete.cases(build_rank_matrix(gran$df, metric)), ],
        conf.level = 0.95,
        plottype   = "vmcb",
        main       = sprintf("%s\n%s", metric_labels[metric], gran$label)
      )
      
      # Summary table
      mean_ranks <- res$means   # named vector of average ranks
      tbl <- data.frame(
        metric       = metric,
        granularity  = gran$label,
        method       = names(mean_ranks),
        mean_rank    = as.numeric(mean_ranks),
        stringsAsFactors = FALSE
      ) %>%
        arrange(mean_rank)
      
      # Append p-value matrix as a long table
      pvals <- res$p.value
      if (!is.null(pvals)) {
        pval_long <- as.data.frame(as.table(pvals)) %>%
          setNames(c("method_1", "method_2", "p_value")) %>%
          mutate(
            metric      = metric,
            granularity = gran$label
          ) %>%
          filter(!is.na(p_value))
      } else {
        pval_long <- NULL
      }
      
      summary_tables[[key]] <- list(ranks = tbl, pvalues = pval_long)
      
      cat("  Average ranks:\n")
      print(tbl)
      if (!is.null(pval_long)) {
        cat("  Pairwise p-values:\n")
        print(pval_long)
      }
      
    } else {
      # Blank plot placeholder
      plot.new()
      title(main = sprintf("%s\n%s\n[FAILED]",
                           metric_labels[metric], gran$label))
    }
  }
  
  dev.off()
  cat(sprintf("\nSaved: %s\n", pdf_path))
}


# ── 8. Save summary tables as CSV ─────────────────────────────────────────────

all_ranks  <- bind_rows(lapply(summary_tables, function(x) x$ranks))
all_pvals  <- bind_rows(lapply(summary_tables, function(x) x$pvalues))

write.csv(all_ranks,  file.path(output_dir, "mcb_mean_ranks.csv"),   row.names = FALSE)
write.csv(all_pvals,  file.path(output_dir, "mcb_pairwise_pvals.csv"), row.names = FALSE)

cat("\nSaved CSVs:\n")
cat(" ", file.path(output_dir, "mcb_mean_ranks.csv"), "\n")
cat(" ", file.path(output_dir, "mcb_pairwise_pvals.csv"), "\n")


# ── 9. Save full results as RData ──────────────────────────────────────────────

save(all_mcb_results, all_ranks, all_pvals, all_df,
     file = file.path(output_dir, "mcb_results_coherence.RData"))

cat("\nSaved RData:", file.path(output_dir, "mcb_results_coherence.RData"), "\n")
cat("\nDone.\n")