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
  "accuracy_purposes_independent_v3",
  "accuracy_purposes_independent_MinT_v3",
  "accuracy_purposes_hierarchy_v2",
  "accuracy_purposes_hierarchy_MinT_v2",
  "accuracy_purposes_hierarchy_coherency0.005v3",
  "accuracy_purposes_hierarchy_coherency0.002v2",
  "accuracy_purposes_hierarchy_coherency0.002_totalweight_nrp",
  "accuracy_purposes_hierarchy_coherency0.002_purposesweight_nregions"
)

stopifnot(length(methods) == length(folders))

dir_out <- "comparison_plots_updated"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Base tiers (EXCLUDING "Overall" from the cross-sectional comprehensive list)
base_levels <- c("Total", "Purpose", "StatePurpose", "ZonePurpose", "RegionPurpose")

level_labels <- c(
  "Total"         = "Total (Australia)",
  "Purpose"       = "Purposes",
  "StatePurpose"  = "States",
  "ZonePurpose"   = "Zones",
  "RegionPurpose" = "Regions"
)

# =============================================================================
# 1. COMPUTE DYNAMIC IN-SAMPLE SEASONAL NAÏVE SCALES (m = 12, T = 96)
# =============================================================================

message("Generating dynamic in-sample seasonal baseline matrices from tourism.Rdata...")
load("tourism.Rdata") # Contains 'tourism' data frame

# Reconstruct exact tracking index sets matching your HPC simulation script
window_size <- 96
test_set    <- 121
origins     <- 96:(96 + test_set - 1) # Vector of 121 loop origins
m_period    <- 12L

# FIX: Select ONLY columns that are numeric before converting to a matrix
tourism_numeric <- tourism %>% dplyr::select(where(is.numeric))
tourism_mat     <- as.matrix(tourism_numeric)

# Process every timestep origin window sequentially to capture active scaling constants
in_sample_seasonal_scales <- tibble::tibble(timestep = origins) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_vector = list({
      # Extract the explicit historical training slice assigned to this rolling origin
      train_idx   <- (timestep - window_size + 1):timestep
      train_slice <- tourism_mat[train_idx, ]
      
      # Compute seasonal differences: Y_t - Y_{t-12}
      n_rows       <- nrow(train_slice)
      seas_diff_sq <- (train_slice[(m_period + 1):n_rows, , drop = FALSE] - 
                         train_slice[1:(n_rows - m_period), , drop = FALSE])^2
      
      # Compute mean squared error of the seasonal persistence model
      colMeans(seas_diff_sq, na.rm = TRUE)
    })
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest_wider(scale_vector) %>%
  tidyr::pivot_longer(
    cols = -timestep,
    names_to = "series",
    values_to = "in_sample_seasonal_naive_mse"
  )

# =============================================================================
# 2. LOAD RAW DATA FUNCTION (Merges Seasonal Scalers)
# =============================================================================

load_method <- function(folder, method_name, scale_df) {
  acc_files <- list.files(folder, pattern = "^accuracy_timestep_[0-9]+\\.rds$", full.names = TRUE)
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
  
  uni_clean <- all_uni %>%
    dplyr::mutate(LogS = ifelse(is.finite(logS), logS, NA_real_)) %>%
    dplyr::left_join(scale_df, by = c("timestep", "series")) %>%
    dplyr::select(level, series, horizon, timestep, ae, se, in_sample_seasonal_naive_mse, crps, LogS, covered_90, covered_95)
  
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
# 3. DATA INGESTION LOOPS
# =============================================================================

all_uni_list   <- list()
all_multi_list <- list()

for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  res <- load_method(folders[i], methods[i], in_sample_seasonal_scales)
  if (!is.null(res)) {
    all_uni_list[[i]]   <- res$uni
    all_multi_list[[i]] <- res$multi
  }
}

all_uni_raw   <- dplyr::bind_rows(all_uni_list)
all_multi_raw <- dplyr::bind_rows(all_multi_list)

if (nrow(all_uni_raw) == 0) {
  stop("Critical Error: No data was loaded. Check paths and filenames.")
}

# =============================================================================
# 4. COMPUTE TRUE SERIES-LEVEL OUT-OF-SAMPLE METRICS & RMSSE
# =============================================================================

# Step 4a: Individual Series Asset Summaries over all windows
series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MAE   = mean(ae, na.rm = TRUE),
    Series_MSE   = mean(se, na.rm = TRUE),
    Series_RMSE  = sqrt(Series_MSE),
    
    # TRUE SEASONAL RMSSE: Root of (Out-Of-Sample MSE / In-Sample Seasonal Naive MSE)
    Series_RMSSE = sqrt(Series_MSE / ifelse(in_sample_seasonal_naive_mse == 0, NA_real_, in_sample_seasonal_naive_mse)),
    
    Series_CRPS  = mean(crps, na.rm = TRUE),
    Series_LogS  = mean(LogS, na.rm = TRUE),
    Series_Cov90 = mean(covered_90, na.rm = TRUE) * 100,
    Series_Cov95 = mean(covered_95, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Step 4b: Macro Tier Aggregations (Excluding Overall)
level_horizon_summary <- series_horizon_metrics %>%
  dplyr::filter(level %in% base_levels) %>%
  dplyr::group_by(method, level, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE),
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# 5. NEW STEP: ISOLATE SYSTEM-WIDE OVERALL METRICS TABLE (WITH RMSSE & ES)
# =============================================================================

overall_system_summary <- series_horizon_metrics %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSSE = mean(Series_RMSSE, na.rm = TRUE), # Safe cross-level seasonal scaled metric
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# Extract and link multivariate global Energy Scores (ES)
multivariate_global_summary <- all_multi_raw %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(ES = mean(es, na.rm = TRUE), .groups = "drop")

overall_system_table <- overall_system_summary %>%
  dplyr::left_join(multivariate_global_summary, by = "method") %>%
  dplyr::mutate(Method = factor(method, levels = methods)) %>%
  dplyr::select(Method, MAE, RMSSE, CRPS, LogS, Cov90, Cov95, ES) %>%
  dplyr::arrange(Method)

# =============================================================================
# 6. COMPUTE TIERED PERFORMANCE SUMMARY TABLE (Excluding Overall Level)
# =============================================================================

univariate_macro_summary <- level_horizon_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE   = mean(MAE, na.rm = TRUE),
    RMSE  = mean(RMSE, na.rm = TRUE),
    CRPS  = mean(CRPS, na.rm = TRUE),
    Cov90 = mean(Cov90, na.rm = TRUE),
    Cov95 = mean(Cov95, na.rm = TRUE),
    .groups = "drop"
  )

master_grid <- expand.grid(method = methods, level = base_levels, stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble()

comprehensive_summary_table <- master_grid %>%
  dplyr::left_join(univariate_macro_summary, by = c("method", "level")) %>%
  dplyr::mutate(
    Presentation_Level = factor(level_labels[level], levels = level_labels),
    Method             = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Presentation_Level, Method) %>%
  dplyr::select(Level = Presentation_Level, Method, MAE, RMSE, CRPS, Cov90, Cov95)

# =============================================================================
# 7. PRINT & EXPORT FINAL DATA OUTPUT MATRICES
# =============================================================================

cat("\n--- OVERALL SYSTEM SCALE-SAFE PERFORMANCE MATRIX (HORIZONS 1-12 WITH RMSSE & ES) ---\n")
overall_system_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n--- TIERED UNIVARIATE SUMMARY MATRIX (EXCLUDING OVERALL LEVEL) ---\n")
comprehensive_summary_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

# Save Files to Export Path
write.csv(overall_system_table,        file.path(dir_out, "overall_system_scale_safe_table.csv"), row.names = FALSE)
write.csv(comprehensive_summary_table, file.path(dir_out, "comprehensive_hierarchical_performance_table.csv"), row.names = FALSE)

# =============================================================================
# 8. RE-ENGINEER COHESION SEGMENTS FOR ALL TIER COMBINATIONS (Wide Format)
# =============================================================================

clean_pair_labels <- function(raw_col) {
  raw_col %>%
    stringr::str_remove("^(mai_|jsd_)") %>%
    stringr::str_replace_all("_", " -> ") %>% 
    stringr::str_to_upper()
}

long_divergences <- all_multi_raw %>%
  dplyr::select(-any_of(c("es", "horizon"))) %>%
  dplyr::mutate(across(starts_with("mai_"), abs)) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::pivot_longer(cols = -method, names_to = "raw_metric", values_to = "score") %>%
  dplyr::mutate(
    Metric_Type     = ifelse(stringr::str_starts(raw_metric, "mai_"), "MAI", "JSD"),
    Structural_Pair = clean_pair_labels(raw_metric),
    Method          = factor(method, levels = methods)
  )

mai_wide_table <- long_divergences %>%
  dplyr::filter(Metric_Type == "MAI") %>%
  dplyr::select(Method, Structural_Pair, score) %>%
  tidyr::pivot_wider(names_from = Structural_Pair, values_from = score) %>%
  dplyr::arrange(Method)

jsd_wide_table <- long_divergences %>%
  dplyr::filter(Metric_Type == "JSD") %>%
  dplyr::select(Method, Structural_Pair, score) %>%
  tidyr::pivot_wider(names_from = Structural_Pair, values_from = score) %>%
  dplyr::arrange(Method)

cat("\n--- MEAN ABSOLUTE INCOHERENCY (MAI) COHESION WIDE MATRIX ---\n")
mai_wide_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

cat("\n--- JENSEN-SHANNON DIVERGENCE (JSD) COHESION WIDE MATRIX ---\n")
jsd_wide_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

write.csv(mai_wide_table, file.path(dir_out, "cohesion_wide_MAI.csv"), row.names = FALSE)
write.csv(jsd_wide_table, file.path(dir_out, "cohesion_wide_JSD.csv"), row.names = FALSE)

message("\nAll rolling window seasonal metrics calculated and exported successfully.")