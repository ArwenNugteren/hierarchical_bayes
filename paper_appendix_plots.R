library(dplyr)
library(tidyr)
library(stringr)

# =============================================================================
# 0. CONFIGURATION & TARGET HORIZONS
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
  "accuracy_purposes_independent",
  "accuracy_purposes_independent_MinT",
  "accuracy_purposes_hierarchy_v2",
  "accuracy_purposes_hierarchy_MinT",
  "accuracy_purposes_hierarchy_coherency0.005",
  "accuracy_purposes_hierarchy_coherency0.002",
  "accuracy_purposes_hierarchy_coherency0.002_totalweight_nrp",
  "accuracy_purposes_hierarchy_coherency0.002_purposesweight_nregions"
)

dir_out <- "comparison_plots"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

base_levels <- c("Overall", "Total", "Purpose", "StatePurpose", "ZonePurpose", "RegionPurpose")

level_labels <- c(
  "Overall"       = "Overall System",
  "Total"         = "Total (Australia)",
  "Purpose"       = "Purposes",
  "StatePurpose"  = "States",
  "ZonePurpose"   = "Zones",
  "RegionPurpose" = "Regions"
)

# =============================================================================
# 1. DATA INGESTION ENGINE (Preserves Granular Windows For Correct Math Progression)
# =============================================================================

load_method <- function(folder, method_name) {
  acc_files <- list.files(folder, pattern = "^accuracy_timestep_[0-9]+\\.rds$", full.names = TRUE)
  if (length(acc_files) == 0) return(NULL)
  
  load_and_tag <- function(f) {
    ts_num <- as.integer(stringr::str_extract(basename(f), "[0-9]+"))
    res    <- readRDS(f)
    uni    <- res$univariate %>% dplyr::mutate(timestep = ts_num)
    list(uni = uni)
  }
  
  loaded_list <- lapply(acc_files, load_and_tag)
  all_uni     <- dplyr::bind_rows(lapply(loaded_list, `[[`, "uni"))
  
  # Structural Tier Definitions extracted long (Keep raw se and crps for rolling windows)
  uni_clean <- all_uni %>%
    dplyr::select(level, series, horizon, timestep, se, crps) %>%
    dplyr::mutate(method = method_name)
  
  return(uni_clean)
}

# Run Ingestion Loop
all_data_list <- list()
for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  all_data_list[[i]] <- load_method(folders[i], methods[i])
}
all_data_raw <- dplyr::bind_rows(all_data_list)

# =============================================================================
# 2. CRITICAL ENGINE STEP: TIME-SERIES ROLLING WINDOW AGGREGATION
# =============================================================================

# Step A: Compute accurate summary performance vectors for each individual series asset 
series_horizon_metrics <- all_data_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_RMSE  = sqrt(mean(se, na.rm = TRUE)), # Correct: Root of mean squared error over windows
    Series_CRPS  = mean(crps, na.rm = TRUE),
    .groups      = "drop"
  )

# Step B: Take the simple average across all individual series within their tier tags
level_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, level, horizon) %>%
  dplyr::summarise(
    Avg_RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple average over asset nodes
    Avg_CRPS  = mean(Series_CRPS, na.rm = TRUE),
    .groups   = "drop"
  )

# Step C: Build systemic "Overall" level by averaging all nodes across the entire model structure
overall_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    Avg_RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple average across all system asset nodes
    Avg_CRPS  = mean(Series_CRPS, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

# Stitch elements together to produce the final tracking base
horizon_base_summary <- dplyr::bind_rows(level_horizon_summary, overall_horizon_summary)

# =============================================================================
# 3. GRANULAR HORIZON SUMMARY PROCESSING ENGINE
# =============================================================================

# Step A: Compute the custom multi-horizon average blocks (1-6 and 1-12)
horizon_averages_summary <- horizon_base_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    RMSE_h1_6  = mean(Avg_RMSE[horizon %in% 1:6],  na.rm = TRUE),
    RMSE_h1_12 = mean(Avg_RMSE[horizon %in% 1:12], na.rm = TRUE),
    CRPS_h1_6  = mean(Avg_CRPS[horizon %in% 1:6],  na.rm = TRUE),
    CRPS_h1_12 = mean(Avg_CRPS[horizon %in% 1:12], na.rm = TRUE),
    .groups = "drop"
  )

# Step B: Filter base summaries to requested milestone horizons and pivot wide
pivot_milestones <- function(df, target_metric) {
  df %>%
    dplyr::filter(horizon %in% c(1, 2, 3, 6, 12)) %>%
    dplyr::mutate(horizon_lbl = paste0("h", horizon)) %>%
    dplyr::select(method, level, horizon_lbl, value = dplyr::all_of(target_metric)) %>%
    tidyr::pivot_wider(names_from = horizon_lbl, values_from = value)
}

rmse_milestones_wide <- pivot_milestones(horizon_base_summary, "Avg_RMSE")
crps_milestones_wide <- pivot_milestones(horizon_base_summary, "Avg_CRPS")

# =============================================================================
# 4. BUILD MASTER RMSE REPORTING FRAME
# =============================================================================

rmse_table_master <- rmse_milestones_wide %>%
  dplyr::left_join(
    horizon_averages_summary %>% dplyr::select(method, level, `h1-6 Avg` = RMSE_h1_6, `h1-12 Avg` = RMSE_h1_12),
    by = c("method", "level")
  ) %>%
  dplyr::mutate(
    Level  = factor(level_labels[level], levels = level_labels),
    Method = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Level, Method) %>%
  dplyr::select(Level, Method, h1, h2, h3, h4 = h6, h12, `h1-6 Avg`, `h1-12 Avg`) %>%
  dplyr::rename(`h6` = h4)

# =============================================================================
# 5. BUILD MASTER CRPS REPORTING FRAME
# =============================================================================

crps_table_master <- crps_milestones_wide %>%
  dplyr::left_join(
    horizon_averages_summary %>% dplyr::select(method, level, `h1-6 Avg` = CRPS_h1_6, `h1-12 Avg` = CRPS_h1_12),
    by = c("method", "level")
  ) %>%
  dplyr::mutate(
    Level  = factor(level_labels[level], levels = level_labels),
    Method = factor(method, levels = methods)
  ) %>%
  dplyr::arrange(Level, Method) %>%
  dplyr::select(Level, Method, h1, h2, h3, h4 = h6, h12, `h1-6 Avg`, `h1-12 Avg`) %>%
  dplyr::rename(`h6` = h4)

# =============================================================================
# 6. PRINT REPORT OUT AND EXPORT FILES
# =============================================================================

cat("\n--- APPENDIX HORIZON TABLE: RMSE ---\n")
rmse_table_master %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

write.csv(rmse_table_master, file.path(dir_out, "separate_summary_table_RMSE.csv"), row.names = FALSE)

cat("\n--- APPENDIX HORIZON TABLE: CRPS ---\n")
crps_table_master %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

write.csv(crps_table_master, file.path(dir_out, "separate_summary_table_CRPS.csv"), row.names = FALSE)