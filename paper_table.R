library(dplyr)
library(tidyr)
library(stringr)

# =============================================================================
# 0. CONFIGURATION & PATHS
# =============================================================================

methods <- c(
  "Independent Bayes",
  "Independent Bayes MinT",
  "Hierarchical Bayes",
  "Hierarchical Bayes MinT",
  "Hierarchical Bayes softer coherency",
  "Hierarchical Bayes tighter coherency",
  "Hierarchical Bayes total-weighted",
  "Hierarchical Bayes purposes-weighted"
)

folders <- c(
  "accuracy_purposes_independentv2",
  "accuracy_purposes_independent_MinTv2",
  "accuracy_purposes_hierarchy_v2",
  "accuracy_purposes_hierarchy_MinT",
  "accuracy_purposes_hierarchy_coherency0.005v2",
  "accuracy_purposes_hierarchy_coherency0.002v2",
  "accuracy_purposes_hierarchy_coherency0.002_totalweight_nrp",
  "accuracy_purposes_hierarchy_coherency0.002_purposesweight_nregions"
)

stopifnot(length(methods) == length(folders))

dir_out <- "comparison_plots_updated"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Structural hierarchy ordering matching original raw factors
base_levels <- c("Overall", "Total", "Purpose", "StatePurpose", "ZonePurpose", "RegionPurpose")

# Clean presentation labels lookup dictionary
level_labels <- c(
  "Overall"       = "Overall System",
  "Total"         = "Total (Australia)",
  "Purpose"       = "Purposes",
  "StatePurpose"  = "States",
  "ZonePurpose"   = "Zones",
  "RegionPurpose" = "Regions"
)

# =============================================================================
# 1. LOAD RAW DATA FUNCTION (Preserves Series & Windows for True RMSE Calculation)
# =============================================================================

load_method <- function(folder, method_name) {
  acc_files <- list.files(folder,
                          pattern    = "^accuracy_timestep_[0-9]+\\.rds$",
                          full.names = TRUE)
  if (length(acc_files) == 0) {
    warning("No accuracy files found in: ", folder)
    return(NULL)
  }
  
  load_and_tag <- function(f) {
    ts_num <- as.integer(stringr::str_extract(basename(f), "[0-9]+"))
    res    <- readRDS(f)
    uni    <- res$univariate  %>% dplyr::mutate(timestep = ts_num)
    multi  <- res$multivariate %>% dplyr::mutate(timestep = ts_num)
    list(uni = uni, multi = multi)
  }
  
  loaded_list <- lapply(acc_files, load_and_tag)
  all_uni     <- dplyr::bind_rows(lapply(loaded_list, `[[`, "uni"))
  all_multi   <- dplyr::bind_rows(lapply(loaded_list, `[[`, "multi"))
  
  # Extract long series logs without averaging timesteps inside the file reader
  uni_clean <- all_uni %>%
    dplyr::mutate(LogS = ifelse(is.finite(logS), logS, NA_real_)) %>%
    dplyr::select(level, series, horizon, timestep, ae, se, crps, LogS, covered_90, covered_95)
  
  # Multivariate metrics are structural system-wide values; average over windows here safely
  multi_clean <- all_multi %>%
    dplyr::group_by(horizon) %>%
    dplyr::summarise(
      across(starts_with(c("mai_", "jsd_")), ~ mean(abs(.x), na.rm = TRUE)),
      es = mean(es, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(method = method_name)
  
  list(uni = uni_clean %>% dplyr::mutate(method = method_name), multi = multi_clean)
}

# =============================================================================
# 2. DATA INGESTION LOOPS
# =============================================================================

all_uni_list   <- list()
all_multi_list <- list()

for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  res <- load_method(folders[i], methods[i])
  if (!is.null(res)) {
    all_uni_list[[i]]   <- res$uni
    all_multi_list[[i]] <- res$multi
  }
}

all_uni_raw   <- dplyr::bind_rows(all_uni_list)
all_multi_raw <- dplyr::bind_rows(all_multi_list)

if (nrow(all_uni_raw) == 0) {
  stop("Critical Error: No data was loaded. Check your folder paths and RDS filenames.")
}

# =============================================================================
# 3. CRITICAL ENGINE STEP: TRUE SERIES-LEVEL RMSE OVER 121 ROLLING WINDOWS
# =============================================================================

# Step 3a: Compute accurate summary performance vectors for each individual series asset 
series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MAE   = mean(ae, na.rm = TRUE),
    Series_RMSE  = sqrt(mean(se, na.rm = TRUE)), # Correct: Root of mean squared error over windows
    Series_CRPS  = mean(crps, na.rm = TRUE),
    Series_LogS  = mean(LogS, na.rm = TRUE),
    Series_Cov90 = mean(covered_90, na.rm = TRUE) * 100,
    Series_Cov95 = mean(covered_95, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Step 3b: Take the simple average across all individual series within their tier tags
level_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, level, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple average over series
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3c: Build system-wide "Overall" level by averaging all nodes across the entire model structure
overall_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple average across all system series nodes
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

# Stitch elements to formulate final tracking structures
master_horizon_summary <- dplyr::bind_rows(level_horizon_summary, overall_horizon_summary)

# =============================================================================
# 4. COMPUTE MASTER PERFORMANCE SUMMARY TABLE (Averaged Across Horizons)
# =============================================================================

univariate_macro_summary <- master_horizon_summary %>%
  dplyr::filter(level %in% base_levels) %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE   = mean(MAE, na.rm = TRUE),
    RMSE  = mean(RMSE, na.rm = TRUE), # Simple average across horizons 1-12
    CRPS  = mean(CRPS, na.rm = TRUE),
    Cov90 = mean(Cov90, na.rm = TRUE),
    Cov95 = mean(Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# Fetch system macro divergences for the "Overall" summary matrix rows
multivariate_overall <- all_multi_raw %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    MAI = mean(abs(mai_rp_total), na.rm = TRUE),
    JSD = mean(jsd_rp_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

master_grid <- expand.grid(method = methods, level = base_levels, stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble()

comprehensive_summary_table <- master_grid %>%
  dplyr::left_join(univariate_macro_summary, by = c("method", "level")) %>%
  dplyr::left_join(multivariate_overall, by = c("method", "level")) %>%
  dplyr::mutate(
    Presentation_Level = factor(level_labels[level], levels = level_labels),
    method             = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Presentation_Level, method) %>%
  dplyr::select(Level = Presentation_Level, Method = method, MAE, RMSE, CRPS, Cov90, Cov95, MAI, JSD)

# =============================================================================
# 5. PRINT & EXPORT CORE TABLES
# =============================================================================

cat("\n--- MASTER PERFORMANCE SUMMARY BY HIERARCHICAL TIER ---\n")
comprehensive_summary_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

write.csv(comprehensive_summary_table, 
          file.path(dir_out, "comprehensive_hierarchical_performance_table.csv"), 
          row.names = FALSE)

# Generate Legacy Width Coverage Matrix Layout
coverage_95_by_level <- comprehensive_summary_table %>%
  dplyr::select(Level, Method, Cov95) %>%
  tidyr::pivot_wider(names_from = Level, values_from = Cov95)

cat("\n--- 95% EMPIRICAL COVERAGE BY TIER VERIFICATION (%) ---\n")
coverage_95_by_level %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  print(n = Inf)

write.csv(coverage_95_by_level, 
          file.path(dir_out, "coverage_95_by_level_summary.csv"), 
          row.names = FALSE)

# =============================================================================
# 6. RE-ENGINEER COHESION SEGMENTS FOR ALL TIER COMBINATIONS (Wide Format)
# =============================================================================

# =============================================================================
# REPAIRED COHESION SEGMENTS FOR ALL TIER COMBINATIONS
# =============================================================================

clean_pair_labels <- function(raw_col) {
  raw_col %>%
    stringr::str_remove("^(mai_|jsd_)") %>%
    # Use explicit anchors or replace all instances safely so names don't collapse
    stringr::str_replace_all("_", " -> ") %>% 
    stringr::str_to_upper()
}

long_divergences <- all_multi_raw %>%
  dplyr::select(-any_of(c("es", "horizon"))) %>%
  # Apply absolute values to all MAI components prior to summarizing
  dplyr::mutate(across(starts_with("mai_"), abs)) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::pivot_longer(cols = -method, names_to = "raw_metric", values_to = "score") %>%
  dplyr::mutate(
    Metric_Type     = ifelse(stringr::str_starts(raw_metric, "mai_"), "MAI", "JSD"),
    Structural_Pair = clean_pair_labels(raw_metric),
    method          = factor(method, levels = methods)
  )

# MAI Wide Table Execution
mai_wide_table <- long_divergences %>%
  dplyr::filter(Metric_Type == "MAI") %>%
  dplyr::select(Method = method, Structural_Pair, score) %>%
  tidyr::pivot_wider(names_from = Structural_Pair, values_from = score) %>%
  dplyr::arrange(Method)

# JSD Wide Table Execution
jsd_wide_table <- long_divergences %>%
  dplyr::filter(Metric_Type == "JSD") %>%
  dplyr::select(Method = method, Structural_Pair, score) %>%
  tidyr::pivot_wider(names_from = Structural_Pair, values_from = score) %>%
  dplyr::arrange(Method)

master_wide_export <- dplyr::bind_rows(
  dplyr::mutate(mai_wide_table, Metric = "MAI") %>% dplyr::select(Metric, Method, everything()),
  dplyr::mutate(jsd_wide_table, Metric = "JSD") %>% dplyr::select(Metric, Method, everything())
)

cat("\n--- MEAN ABSOLUTE INCOHERENCY (MAI) COHESION WIDE MATRIX ---\n")
mai_wide_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

cat("\n--- JENSEN-SHANNON DIVERGENCE (JSD) COHESION WIDE MATRIX ---\n")
jsd_wide_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

write.csv(master_wide_export, 
          file.path(dir_out, "structural_coherency_wide_table.csv"), 
          row.names = FALSE)

message("\nAll 121-window data paths synchronized and compiled successfully.")