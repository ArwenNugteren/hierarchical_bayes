library(dplyr)
library(tidyr)
library(stringr)

# =============================================================================
# 0. CONFIGURATION & PATHS
# =============================================================================

# methods <- c(
#   "Independent Bayes - 10",
#   # "Independent Bayes MinT - 10",
#   "Hierarchical Bayes - 10",
#   # "Hierarchical Bayes MinT - 10",
#   "Independent Bayes - 20",
#   # "Independent Bayes MinT - 20",
#   "Hierarchical Bayes - 20",
#   # "Hierarchical Bayes MinT - 20",
#   "Independent Bayes - 50",
#   # "Independent Bayes MinT - 50",
#   "Hierarchical Bayes - 50"
#   # "Hierarchical Bayes MinT - 50"
# )
# 
# folders <- c(
#   "accuracy_misspecified_independent10",
#   # "accuracy_misspecified_independent10_mint",
#   "accuracy_misspecified_hierarchy10",
#   # "accuracy_misspecified_hierarchy10_mint",
#   "accuracy_misspecified_independent20",
#   # "accuracy_misspecified_independent20_mint",
#   "accuracy_misspecified_hierarchy20",
#   # "accuracy_misspecified_hierarchy20_mint",
#   "accuracy_misspecified_independent_v2",
#   # "accuracy_misspecified_independent_mint_v2",
#   "accuracy_misspecified_hierarchical_v3"
#   # "accuracy_misspecified_hierarchical_mint_v2"
# )
methods <- c(
  "Independent Bayes",
  "Independent Bayes MinT",
  "Hierarchical Bayes",
  "Hierarchical Bayes MinT",
  "Hierarchical Bayes softer coherency",
  "Hierarchical Bayes tighter coherency",
  "Hierarchical Bayes total-weighted",
  "Hierarchical Bayes bottom-weighted"
)

folders <- c(
  "accuracy_misspecified_independent_v2",
  "accuracy_misspecified_independent_mint_v2",
  "accuracy_misspecified_hierarchical_v3",
  "accuracy_misspecified_hierarchical_mint_v3",
  "accuracy_misspecified_coherency1_v3",
  "accuracy_misspecified_coherencysqrt0.5_v3",
  "accuracy_misspecified_coherency0.707_topweighted_v3",
  "accuracy_misspecified_coherency0.707_bottomweighted_v2"
)


stopifnot(length(methods) == length(folders))

dir_out   <- "comparison_plots_missspecified"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

horizon   <- 5L   
target_hs <- c(1, 2, 3, 4, 5) 

# Structural constants matching the training configuration
TRAIN_SIZE <- 50L

# Structural tiers (Excluding "Overall" entirely from cross-sectional tables)
base_levels  <- c("top", "mid", "bottom")
level_labels <- c(
  "top"    = "Top Level",
  "mid"    = "Middle Level",
  "bottom" = "Bottom Level"
)

# =============================================================================
# 1. GENERATE TRUE IN-SAMPLE NAÏVE BASELINE SCALES FROM HISTORICAL WINDOWS
# =============================================================================
synthetic_data <- read.csv("synthetic_hierarchy.csv")
colnames(synthetic_data) <- c("Total", "A", "B", "AA", "AB", "BA", "BB")

# Convert to matrix for highly efficient slicing inside the rolling loop
data_mat_raw <- as.matrix(synthetic_data)

# Compute the mean squared in-sample first-differences for each training window
# Each roll window starts at index: roll, and ends at: roll + TRAIN_SIZE - 1
in_sample_scales <- tibble::tibble(roll = 1L:500L) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_vector = list({
      # Extract the exact 50-step historical matrix used to train this specific roll
      train_window <- data_mat_raw[roll:(roll + TRAIN_SIZE - 1L), ]
      
      # Calculate (Y_t - Y_{t-1})^2 across the training matrix
      # diff calculates the difference between each element and the next
      diff_sq <- diff(train_window)^2
      
      # Compute the mean squared error of this historical persistence baseline
      colMeans(diff_sq, na.rm = TRUE)
    })
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest_wider(scale_vector) %>%
  tidyr::pivot_longer(
    cols = -roll,
    names_to = "series",
    values_to = "in_sample_naive_mse"
  )

# =============================================================================
# 2. LOAD RAW DATA FUNCTION (Merges In-Sample Scaling Constants)
# =============================================================================

load_method <- function(folder, method_name, scale_df) {
  acc_files <- list.files(folder, pattern = "^accuracy_roll_[0-9]+\\.rds$", full.names = TRUE)
  if (length(acc_files) == 0) {
    warning("No accuracy files found in: ", folder)
    return(NULL)
  }
  
  load_and_tag <- function(f) {
    roll_num <- as.integer(stringr::str_extract(basename(f), "[0-9]+"))
    res      <- readRDS(f)
    uni      <- res$univariate   %>% dplyr::mutate(roll = roll_num)
    multi    <- res$multivariate %>% dplyr::mutate(roll = roll_num)
    list(uni = uni, multi = multi)
  }
  
  loaded_list <- lapply(acc_files, load_and_tag)
  all_uni     <- dplyr::bind_rows(lapply(loaded_list, `[[`, "uni"))
  all_multi   <- dplyr::bind_rows(lapply(loaded_list, `[[`, "multi"))
  
  # Inject the historical scales into the model's out-of-sample tracking frame
  uni_clean <- all_uni %>%
    dplyr::left_join(scale_df, by = c("roll", "series")) %>%
    dplyr::mutate(
      LogS = ifelse(is.finite(logS), logS, NA_real_)
    ) %>%
    # Explicitly retain model squared error (se) and historical base scale
    dplyr::select(level, series, horizon, roll, ae, se, in_sample_naive_mse, crps, LogS, covered_90, covered_95)
  
  multi_clean <- all_multi %>%
    dplyr::group_by(horizon) %>%
    dplyr::summarise(
      ES          = mean(es, na.rm = TRUE),
      MAI_bot_mid = mean(abs(mai_bot_mid), na.rm = TRUE),
      MAI_mid_top = mean(abs(mai_mid_top), na.rm = TRUE),
      MAI_bot_top = mean(abs(mai_bot_top), na.rm = TRUE),
      JSD_bot_mid = mean(jsd_bot_mid, na.rm = TRUE),
      JSD_mid_top = mean(jsd_mid_top, na.rm = TRUE),
      JSD_bot_top = mean(jsd_bot_top, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(method = method_name)
  
  list(uni = uni_clean %>% dplyr::mutate(method = method_name), multi = multi_clean)
}

# =============================================================================
# 3. INGEST DATA FROM ALL FOLDERS
# =============================================================================

all_uni_list   <- list()
all_multi_list <- list()

for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  # CHANGED: Passed in_sample_scales here instead of naive_baselines
  res <- load_method(folders[i], methods[i], in_sample_scales)
  if (!is.null(res)) {
    all_uni_list[[i]]   <- res$uni
    all_multi_list[[i]] <- res$multi
  }
}

all_uni_raw   <- dplyr::bind_rows(all_uni_list)
all_multi_raw <- dplyr::bind_rows(all_multi_list)

# =============================================================================
# 4. COMPUTE TRUE SERIES-LEVEL RMSE & RMSSE OVER ROLLING WINDOWS
# =============================================================================

# Step 4a: Calculate specific metrics per unique series at each horizon over all rolling windows
series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MAE   = mean(ae, na.rm = TRUE),
    Series_MSE   = mean(se, na.rm = TRUE),
    Series_RMSE  = sqrt(Series_MSE), 
    
    # TRUE RMSSE: Root of (Mean Out-of-Sample MSE / Mean In-Sample Naive MSE)
    Series_RMSSE = sqrt(Series_MSE / ifelse(in_sample_naive_mse == 0, NA_real_, in_sample_naive_mse)),
    
    Series_CRPS  = mean(crps, na.rm = TRUE),
    Series_LogS  = mean(LogS, na.rm = TRUE),
    Series_Cov90 = mean(covered_90, na.rm = TRUE) * 100,
    Series_Cov95 = mean(covered_95, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Step 4b: Average across all series within each specific structural level and horizon
level_horizon_summary <- series_horizon_metrics %>%
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

# Step 4c: Build the scale-safe global "Overall System Table" spanning Horizons 1-5
overall_system_summary <- series_horizon_metrics %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSSE = mean(Series_RMSSE, na.rm = TRUE), # Safe cross-level scale-invariant metric
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# Extract and attach the global multivariate Energy Score (ES) directly 
multivariate_global_summary <- all_multi_raw %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(ES = mean(ES, na.rm = TRUE), .groups = "drop")

overall_system_table <- overall_system_summary %>%
  dplyr::left_join(multivariate_global_summary, by = "method") %>%
  dplyr::mutate(Method = factor(method, levels = methods)) %>%
  dplyr::select(Method, MAE, RMSSE, CRPS, LogS, Cov90, Cov95, ES) %>%
  dplyr::arrange(Method)

# =============================================================================
# 5. COMPUTE COMPREHENSIVE TIERED PERFORMANCE TABLE (Excluding Overall)
# =============================================================================

univariate_macro_summary <- level_horizon_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE   = mean(MAE, na.rm = TRUE),
    RMSE  = mean(RMSE, na.rm = TRUE), 
    CRPS  = mean(CRPS, na.rm = TRUE),
    LogS  = mean(LogS, na.rm = TRUE),
    Cov90 = mean(Cov90, na.rm = TRUE),
    Cov95 = mean(Cov95, na.rm = TRUE),
    .groups = "drop"
  )

master_spine <- expand.grid(method = methods, level = base_levels, stringsAsFactors = FALSE)

comprehensive_summary_table <- master_spine %>%
  dplyr::as_tibble() %>%
  dplyr::left_join(univariate_macro_summary, by = c("method", "level")) %>%
  dplyr::mutate(
    Presentation_Level = factor(level_labels[level], levels = level_labels),
    Method             = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Presentation_Level, Method) %>%
  dplyr::select(Level = Presentation_Level, Method, MAE, RMSE, CRPS, LogS, Cov90, Cov95)

# =============================================================================
# 6. SEPARATE ACCURACY SPREADS (Excluding Overall)
# =============================================================================

window_aggregations <- level_horizon_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    RMSE_h1_5 = mean(RMSE, na.rm = TRUE),
    CRPS_h1_5 = mean(CRPS, na.rm = TRUE),
    .groups = "drop"
  )

pivot_simulation_metrics <- function(df, column_name) {
  df %>%
    dplyr::filter(horizon %in% target_hs) %>%
    dplyr::mutate(h_lbl = paste0("h", horizon)) %>%
    dplyr::select(method, level, h_lbl, value = dplyr::all_of(column_name)) %>%
    tidyr::pivot_wider(names_from = h_lbl, values_from = value)
}

rmse_sim_wide <- pivot_simulation_metrics(level_horizon_summary, "RMSE") %>%
  dplyr::left_join(window_aggregations %>% dplyr::select(method, level, `h1-5 Avg` = RMSE_h1_5), by = c("method", "level"))

crps_sim_wide <- pivot_simulation_metrics(level_horizon_summary, "CRPS") %>%
  dplyr::left_join(window_aggregations %>% dplyr::select(method, level, `h1-5 Avg` = CRPS_h1_5), by = c("method", "level"))

format_master_sim_layout <- function(df) {
  df %>%
    dplyr::mutate(
      Level  = factor(level_labels[level], levels = level_labels),
      Method = factor(method, levels = methods)
    ) %>%
    dplyr::arrange(Level, Method) %>%
    dplyr::select(Level, Method, h1, h2, h3, h5, `h1-5 Avg`)
}

rmse_separated_table <- format_master_sim_layout(rmse_sim_wide)
crps_separated_table <- format_master_sim_layout(crps_sim_wide)

# =============================================================================
# 7. MULTIVARIATE COHESION WIDE TABLES (MAI & JSD)
# =============================================================================

raw_cohesion_means <- all_multi_raw %>%
  tidyr::pivot_longer(cols = starts_with(c("MAI_", "JSD_")), names_to = "raw_metric", values_to = "score") %>%
  dplyr::group_by(method, raw_metric) %>%
  dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    Family = ifelse(stringr::str_starts(raw_metric, "MAI_"), "MAI", "JSD"),
    Pair   = stringr::str_remove(raw_metric, "^(MAI_|JSD_)") %>% stringr::str_to_upper(),
    Method = factor(method, levels = methods)
  )

mai_cohesion_wide <- raw_cohesion_means %>%
  dplyr::filter(Family == "MAI") %>%
  dplyr::select(Method, Pair, score) %>%
  tidyr::pivot_wider(names_from = Pair, values_from = score) %>%
  dplyr::arrange(Method)

jsd_cohesion_wide <- raw_cohesion_means %>%
  dplyr::filter(Family == "JSD") %>%
  dplyr::select(Method, Pair, score) %>%
  tidyr::pivot_wider(names_from = Pair, values_from = score) %>%
  dplyr::arrange(Method)

# =============================================================================
# 8. SYSTEM RENDER PRINT & FILE EXPORT
# =============================================================================

cat("\n--- OVERALL SCALE-SAFE SYSTEM MATRIX (HORIZONS 1-5 SUMMARY WITH RMSSE & ES) ---\n")
overall_system_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n--- COMPREHENSIVE SIMULATION MATRIX BY HIERARCHICAL TIER (EXCLUDING OVERALL) ---\n")
comprehensive_summary_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n--- DETAILED RMSE EXTENSION MATRIX ---\n")
rmse_separated_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n--- DETAILED CRPS EXTENSION MATRIX ---\n")
crps_separated_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n--- MEAN ABSOLUTE INCOHERENCY (MAI) COHESION MATRIX ---\n")
mai_cohesion_wide %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

cat("\n--- JENSEN-SHANNON DIVERGENCE (JSD) COHESION MATRIX ---\n")
jsd_cohesion_wide %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print(n = Inf)

# CSV File Export Sequence
write.csv(overall_system_table,        file.path(dir_out, "overall_system_table.csv"),     row.names = FALSE)
write.csv(comprehensive_summary_table, file.path(dir_out, "comprehensive_summary_hierarchical.csv"), row.names = FALSE)
write.csv(rmse_separated_table,        file.path(dir_out, "separated_horizon_table_RMSE.csv"),        row.names = FALSE)
write.csv(crps_separated_table,        file.path(dir_out, "separated_horizon_table_CRPS.csv"),        row.names = FALSE)
write.csv(mai_cohesion_wide,           file.path(dir_out, "cohesion_wide_MAI.csv"),                  row.names = FALSE)
write.csv(jsd_cohesion_wide,           file.path(dir_out, "cohesion_wide_JSD.csv"),                  row.names = FALSE)