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

# edit with your folder names for the accuracy metrics
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

dir_out <- "paper_plots"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Structural hierarchy ordering matching original raw factors (EXCLUDING "Overall")
base_levels <- c("Total", "Purpose", "StatePurpose", "ZonePurpose", "RegionPurpose")

# Clean presentation labels lookup dictionary
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
load("data/tourism.Rdata") # Loads 'tourism' frame into environment

window_size <- 96
test_set    <- 121
origins     <- 96:(96 + test_set - 1)
m_period    <- 12L

# Isolate numeric metrics exclusively to prevent structural matrix conversion failures
tourism_numeric <- tourism %>% dplyr::select(where(is.numeric))
tourism_mat     <- as.matrix(tourism_numeric)

in_sample_seasonal_scales <- tibble::tibble(timestep = origins) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_vector = list({
      train_idx   <- (timestep - window_size + 1):timestep
      train_slice <- tourism_mat[train_idx, ]
      
      n_rows       <- nrow(train_slice)
      seas_diff_sq <- (train_slice[(m_period + 1):n_rows, , drop = FALSE] - 
                         train_slice[1:(n_rows - m_period), , drop = FALSE])^2
      
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
# 2. DATA INGESTION ENGINE (Preserves Granular Windows and Merges Scalers)
# =============================================================================

load_method <- function(folder, method_name, scale_df) {
  acc_files <- list.files(folder, pattern = "^accuracy_timestep_[0-9]+\\.rds$", full.names = TRUE)
  if (length(acc_files) == 0) return(NULL)
  
  load_and_tag <- function(f) {
    ts_num <- as.integer(stringr::str_extract(basename(f), "[0-9]+"))
    res    <- readRDS(f)
    uni    <- res$univariate   %>% dplyr::mutate(timestep = ts_num)
    multi  <- res$multivariate %>% dplyr::mutate(timestep = ts_num)
    list(uni = uni, multi = multi)
  }
  
  loaded_list <- lapply(acc_files, load_and_tag)
  all_uni     <- dplyr::bind_rows(lapply(loaded_list, `[[`, "uni"))
  all_multi   <- dplyr::bind_rows(lapply(loaded_list, `[[`, "multi"))
  
  uni_clean <- all_uni %>%
    dplyr::select(level, series, horizon, timestep, se, crps) %>%
    dplyr::left_join(scale_df, by = c("timestep", "series")) %>%
    dplyr::mutate(method = method_name)
  
  multi_clean <- all_multi %>%
    dplyr::select(horizon, timestep, es) %>%
    dplyr::mutate(method = method_name)
  
  return(list(uni = uni_clean, multi = multi_clean))
}

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

# =============================================================================
# 3. ROLLING WINDOW TIMESTEP AGGREGATION ENGINE (SERIES LEVEL FIRST)
# =============================================================================

series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MSE   = mean(se, na.rm = TRUE),
    Series_RMSE  = sqrt(Series_MSE),
    Series_RMSSE = sqrt(Series_MSE / ifelse(in_sample_seasonal_naive_mse == 0, NA_real_, in_sample_seasonal_naive_mse)),
    Series_CRPS  = mean(crps, na.rm = TRUE),
    .groups      = "drop"
  )

# =============================================================================
# 4. CONSTRUCT MULTI-HORIZON AVERAGES AND PIVOT UTILITIES
# =============================================================================

# Wide Pivot Function for Targets
pivot_milestones <- function(df, target_metric) {
  df %>%
    dplyr::filter(horizon %in% c(1, 2, 3, 6, 12)) %>%
    dplyr::mutate(horizon_lbl = paste0("h", horizon)) %>%
    dplyr::select(method, level, horizon_lbl, value = dplyr::all_of(target_metric)) %>%
    tidyr::pivot_wider(names_from = horizon_lbl, values_from = value)
}

# =============================================================================
# 5. GENERATE TIERED PERFORMANCE SUMMARIES (EXCLUDING OVERALL LEVEL)
# =============================================================================

level_horizon_summary <- series_horizon_metrics %>%
  dplyr::filter(level %in% base_levels) %>%
  dplyr::group_by(method, level, horizon) %>%
  dplyr::summarise(
    Avg_RMSE  = mean(Series_RMSE, na.rm = TRUE),
    Avg_CRPS  = mean(Series_CRPS, na.rm = TRUE),
    .groups   = "drop"
  )

horizon_averages_summary <- level_horizon_summary %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    RMSE_h1_6   = mean(Avg_RMSE[horizon %in% 1:6],  na.rm = TRUE),
    RMSE_h1_12  = mean(Avg_RMSE[horizon %in% 1:12], na.rm = TRUE),
    CRPS_h1_6   = mean(Avg_CRPS[horizon %in% 1:6],  na.rm = TRUE),
    CRPS_h1_12  = mean(Avg_CRPS[horizon %in% 1:12], na.rm = TRUE),
    .groups     = "drop"
  )

rmse_milestones_wide <- pivot_milestones(level_horizon_summary, "Avg_RMSE")
# FIXED: Changed := to <- 
crps_milestones_wide <- pivot_milestones(level_horizon_summary, "Avg_CRPS")

# Finalize Sub-Tier Structural Frames
rmse_table_master <- rmse_milestones_wide %>%
  dplyr::left_join(horizon_averages_summary %>% dplyr::select(method, level, `h1-6 Avg` = RMSE_h1_6, `h1-12 Avg` = RMSE_h1_12), by = c("method", "level")) %>%
  dplyr::mutate(Level = factor(level_labels[level], levels = level_labels), Method = factor(method, levels = methods)) %>%
  dplyr::arrange(Level, Method) %>%
  dplyr::select(Level, Method, h1, h2, h3, h6, h12, `h1-6 Avg`, `h1-12 Avg`)

crps_table_master <- crps_milestones_wide %>%
  dplyr::left_join(horizon_averages_summary %>% dplyr::select(method, level, `h1-6 Avg` = CRPS_h1_6, `h1-12 Avg` = CRPS_h1_12), by = c("method", "level")) %>%
  dplyr::mutate(Level = factor(level_labels[level], levels = level_labels), Method = factor(method, levels = methods)) %>%
  dplyr::arrange(Level, Method) %>%
  dplyr::select(Level, Method, h1, h2, h3, h6, h12, `h1-6 Avg`, `h1-12 Avg`)

# =============================================================================
# 6. GENERATE ISOLATED SYSTEM-WIDE OVERALL MATRIX (WITH RMSSE & GLOBAL ES)
# =============================================================================

overall_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    Avg_RMSSE = mean(Series_RMSSE, na.rm = TRUE),
    .groups   = "drop"
  )

multivariate_global_summary <- all_multi_raw %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(Avg_ES = mean(es, na.rm = TRUE), .groups = "drop")

overall_master_summary <- overall_horizon_summary %>%
  dplyr::left_join(multivariate_global_summary, by = c("method", "horizon")) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    RMSSE_h1_6   = mean(Avg_RMSSE[horizon %in% 1:6],  na.rm = TRUE),
    RMSSE_h1_12  = mean(Avg_RMSSE[horizon %in% 1:12], na.rm = TRUE),
    ES_h1_6      = mean(Avg_ES[horizon %in% 1:6],     na.rm = TRUE),
    ES_h1_12     = mean(Avg_ES[horizon %in% 1:12],    na.rm = TRUE),
    .groups      = "drop"
  )

overall_rmsse_wide <- overall_horizon_summary %>%
  dplyr::filter(horizon %in% c(1, 2, 3, 6, 12)) %>%
  dplyr::mutate(horizon_lbl = paste0("h", horizon)) %>%
  dplyr::select(method, horizon_lbl, value = Avg_RMSSE) %>%
  tidyr::pivot_wider(names_from = horizon_lbl, values_from = value)

overall_es_wide <- multivariate_global_summary %>%
  dplyr::filter(horizon %in% c(1, 2, 3, 6, 12)) %>%
  dplyr::mutate(horizon_lbl = paste0("h", horizon)) %>%
  dplyr::select(method, horizon_lbl, value = Avg_ES) %>%
  tidyr::pivot_wider(names_from = horizon_lbl, values_from = value)

overall_rmsse_table <- overall_rmsse_wide %>%
  dplyr::left_join(overall_master_summary %>% dplyr::select(method, `h1-6 Avg` = RMSSE_h1_6, `h1-12 Avg` = RMSSE_h1_12), by = "method") %>%
  dplyr::mutate(Method = factor(method, levels = methods)) %>%
  dplyr::arrange(Method) %>%
  dplyr::select(Method, h1, h2, h3, h6, h12, `h1-6 Avg`, `h1-12 Avg`)

overall_es_table <- overall_es_wide %>%
  dplyr::left_join(overall_master_summary %>% dplyr::select(method, `h1-6 Avg` = ES_h1_6, `h1-12 Avg` = ES_h1_12), by = "method") %>%
  dplyr::mutate(Method = factor(method, levels = methods)) %>%
  dplyr::arrange(Method) %>%
  dplyr::select(Method, h1, h2, h3, h6, h12, `h1-6 Avg`, `h1-12 Avg`)

# =============================================================================
# 7. PRINT REPORT OUT AND EXPORT FILES
# =============================================================================

cat("\n--- OVERALL SYSTEM SCALE-SAFE PROFILE MATRIX: RMSSE ---\n")
overall_rmsse_table %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>% print(n = Inf)
write.csv(overall_rmsse_table, file.path(dir_out, "RMSSE.csv"), row.names = FALSE)

cat("\n--- OVERALL SYSTEM MULTIVARIATE PROFILE MATRIX: ENERGY SCORE (ES) ---\n")
overall_es_table %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>% print(n = Inf)
write.csv(overall_es_table, file.path(dir_out, "ES.csv"), row.names = FALSE)

cat("\n--- APPENDIX HORIZON TIER TABLE: RMSE (EXCLUDING OVERALL) ---\n")
rmse_table_master %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>% print(n = Inf)
write.csv(rmse_table_master, file.path(dir_out, "RMSE.csv"), row.names = FALSE)

cat("\n--- APPENDIX HORIZON TIER TABLE: CRPS (EXCLUDING OVERALL) ---\n")
crps_table_master %>% dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>% print(n = Inf)
write.csv(crps_table_master, file.path(dir_out, "CRPS.csv"), row.names = FALSE)

message("\nAll milestone tables restructured and exported safely.")