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
  "Hierarchical Bayes total-weighted"
)

folders <- c(
  "accuracy_independent_wellspecified_joint",
  "accuracy_independent_wellspecified_joint_mint",
  "accuracy_wellspecified_joint",
  "accuracy_wellspecified_joint_mint",
  "accuracy_wellspecified_joint_coherency1",
  "accuracy_wellspecified_joint_coherency0.707",
  "accuracy_wellspecified_joint_coherency0.707_topweighted"
)

stopifnot(length(methods) == length(folders))

dir_out   <- "comparison_plots_wellspecified"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

horizon   <- 5L   
target_hs <- c(1, 2, 3, 5) 

base_levels  <- c("Overall", "top", "mid", "bottom")
level_labels <- c(
  "Overall" = "Overall System",
  "top"     = "Top Level",
  "mid"     = "Middle Level",
  "bottom"  = "Bottom Level"
)

# =============================================================================
# 1. LOAD RAW DATA FUNCTION (Preserves Series and Rolls for Proper RMSE)
# =============================================================================

load_method <- function(folder, method_name) {
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
  
  # Return the granular series data with roll intact so we can compute RMSE properly downstream
  uni_clean <- all_uni %>%
    dplyr::mutate(LogS = ifelse(is.finite(logS), logS, NA_real_)) %>%
    dplyr::select(level, series, horizon, roll, ae, se, crps, LogS, covered_90, covered_95)
  
  # Multivariate metrics are global properties per roll/horizon, so we can average them here safely
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
# 2. INGEST DATA FROM ALL FOLDERS
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

# =============================================================================
# 3. CRITICAL STEP: COMPUTE TRUE SERIES-LEVEL RMSE OVER ROLLING WINDOWS
# =============================================================================

# Step 3a: Calculate specific metrics per unique series at each horizon over all rolling windows
series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MAE   = mean(ae, na.rm = TRUE),
    Series_RMSE  = sqrt(mean(se, na.rm = TRUE)), # Correct: Root of mean squared error over rolls
    Series_CRPS  = mean(crps, na.rm = TRUE),
    Series_LogS  = mean(LogS, na.rm = TRUE),
    Series_Cov90 = mean(covered_90, na.rm = TRUE) * 100,
    Series_Cov95 = mean(covered_95, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Step 3b: Average across all series within each specific structural level and horizon
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

# Step 3c: Build the system-wide "Overall" level by averaging across all series across all levels
overall_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple average over all series in the system
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

# Combine into one master structural baseline frame
master_horizon_summary <- dplyr::bind_rows(level_horizon_summary, overall_horizon_summary)

# =============================================================================
# 4. COMPUTE LARGE HIERARCHICAL COMPREHENSIVE PERFORMANCE TABLE (Averaged Over H)
# =============================================================================

univariate_macro_summary <- master_horizon_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE   = mean(MAE, na.rm = TRUE),
    RMSE  = mean(RMSE, na.rm = TRUE), # Simple average over horizons
    CRPS  = mean(CRPS, na.rm = TRUE),
    LogS  = mean(LogS, na.rm = TRUE),
    Cov90 = mean(Cov90, na.rm = TRUE),
    Cov95 = mean(Cov95, na.rm = TRUE),
    .groups = "drop"
  )

multivariate_global_summary <- all_multi_raw %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(ES = mean(ES, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(level = "Overall")

master_spine <- expand.grid(method = methods, level = base_levels, stringsAsFactors = FALSE)

comprehensive_summary_table <- master_spine %>%
  dplyr::as_tibble() %>%
  dplyr::left_join(univariate_macro_summary, by = c("method", "level")) %>%
  dplyr::left_join(multivariate_global_summary, by = c("method", "level")) %>%
  dplyr::mutate(
    Presentation_Level = factor(level_labels[level], levels = level_labels),
    Method             = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Presentation_Level, Method) %>%
  dplyr::select(Level = Presentation_Level, Method, MAE, RMSE, CRPS, LogS, Cov90, Cov95, ES)

# =============================================================================
# 5. SEPARATE ACCURACY SPREADS (RMSE & CRPS SEPARATED BY TIER & HORIZON)
# =============================================================================

window_aggregations <- master_horizon_summary %>%
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

rmse_sim_wide <- pivot_simulation_metrics(master_horizon_summary, "RMSE") %>%
  dplyr::left_join(window_aggregations %>% dplyr::select(method, level, `h1-5 Avg` = RMSE_h1_5), by = c("method", "level"))

crps_sim_wide <- pivot_simulation_metrics(master_horizon_summary, "CRPS") %>%
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
# 6. MULTIVARIATE COHESION WIDE TABLES (MAI & JSD)
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
# 7. SYSTEM RENDER PRINT & FILE EXPORT
# =============================================================================

cat("\n--- COMPREHENSIVE SIMULATION MATRIX BY HIERARCHICAL LEVEL ---\n")
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
write.csv(comprehensive_summary_table, file.path(dir_out, "comprehensive_summary_hierarchical.csv"), row.names = FALSE)
write.csv(rmse_separated_table,        file.path(dir_out, "separated_horizon_table_RMSE.csv"),        row.names = FALSE)
write.csv(crps_separated_table,        file.path(dir_out, "separated_horizon_table_CRPS.csv"),        row.names = FALSE)
write.csv(mai_cohesion_wide,           file.path(dir_out, "cohesion_wide_MAI.csv"),                  row.names = FALSE)
write.csv(jsd_cohesion_wide,           file.path(dir_out, "cohesion_wide_JSD.csv"),                  row.names = FALSE)

