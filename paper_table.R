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

# =============================================================================
# 1. LOAD RAW DATA FUNCTION
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
    uni    <- res$univariate   %>% dplyr::mutate(timestep = ts_num)
    multi  <- res$multivariate %>% dplyr::mutate(timestep = ts_num)
    list(uni = uni, multi = multi)
  }
  
  loaded_list <- lapply(acc_files, load_and_tag)
  all_uni     <- dplyr::bind_rows(lapply(loaded_list, `[[`, "uni"))
  all_multi   <- dplyr::bind_rows(lapply(loaded_list, `[[`, "multi"))
  
  # ── Per-level rows ────────────────────────────────────────────────────────
  uni_clean <- all_uni %>%
    dplyr::mutate(
      MAE   = ae,
      RMSE  = sqrt(se),
      CRPS  = crps,
      LogS  = ifelse(is.finite(logS), logS, NA_real_),
      Cov90 = covered_90,
      Cov95 = covered_95
    ) %>%
    dplyr::select(level, horizon, timestep, MAE, RMSE, CRPS, LogS, Cov90, Cov95)
  
  # ── Overall rows ──────────────────────────────────────────────────────────
  overall_clean <- all_uni %>%
    dplyr::group_by(horizon, timestep) %>%
    dplyr::summarise(
      MAE   = mean(ae,                          na.rm = TRUE),
      RMSE  = sqrt(mean(se,                     na.rm = TRUE)),
      CRPS  = mean(crps,                        na.rm = TRUE),
      LogS  = mean(ifelse(is.finite(logS), logS, NA_real_), na.rm = TRUE),
      Cov90 = mean(covered_90,                  na.rm = TRUE),
      Cov95 = mean(covered_95,                  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(level = "Overall")
  
  # ── Multivariate rows ─────────────────────────────────────────────────────
  multi_clean <- all_multi %>%
    dplyr::mutate(
      level        = "Overall",
      ES           = es,
      MAI_rp_total = abs(mai_rp_total),
      JSD_rp_total = jsd_rp_total
    ) %>%
    dplyr::select(level, horizon, timestep, ES, MAI_rp_total, JSD_rp_total)
  
  # ── Join and tag ──────────────────────────────────────────────────────────
  dplyr::bind_rows(uni_clean, overall_clean) %>%
    dplyr::left_join(multi_clean, by = c("level", "horizon", "timestep")) %>%
    dplyr::mutate(method = method_name)
}

# =============================================================================
# 2. INGEST DATA FROM ALL FOLDERS
# =============================================================================

all_data_list <- list()
for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  all_data_list[[i]] <- load_method(folders[i], methods[i])
}

all_data <- dplyr::bind_rows(all_data_list)

if (nrow(all_data) == 0) {
  stop("Critical Error: No data was loaded. Check your folder paths and RDS filenames.")
}

# =============================================================================
# 3. COMPUTE SUMMARY PERFORMANCE 
# =============================================================================

univariate_overall <- all_data %>%
  dplyr::filter(level == "Overall") %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    MAE  = mean(MAE, na.rm = TRUE),  
    CRPS  = mean(CRPS, na.rm = TRUE),
    Cov90 = mean(Cov90, na.rm = TRUE),
    Cov95 = mean(Cov95, na.rm = TRUE),
    .groups = "drop"
  )

multivariate_overall <- all_data %>%
  dplyr::filter(level == "Overall") %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    MAI = mean(MAI_rp_total, na.rm = TRUE),
    JSD = mean(JSD_rp_total, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# 4. EXTRACT AND FORMAT SUMMARY TABLE (Overall Level)
# =============================================================================

summary_table <- dplyr::tibble(method = methods) %>%
  dplyr::left_join(univariate_overall, by = "method") %>%
  dplyr::left_join(multivariate_overall, by = "method") %>%
  dplyr::mutate(
    Cov90 = Cov90 * 100,
    Cov95 = Cov95 * 100
  )

# =============================================================================
# 5. OUTPUT OPTIONS
# =============================================================================

summary_table %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  print(n = Inf)

# Export option
write.csv(summary_table, file.path(dir_out, "performance_summary_table.csv"), row.names = FALSE)


# =============================================================================
# 6. EXTRACT AND FORMAT LEVEL-SPECIFIC 95% COVERAGE TABLE
# =============================================================================

# Define your desired level factor levels to enforce the hierarchical order
level_order <- c("Overall", "Total", "Purpose", 
                 "StatePurpose", "ZonePurpose", "RegionPurpose")

coverage_95_by_level <- all_data %>%
  # Filter out any rows missing coverage data
  dplyr::filter(!is.na(Cov95)) %>%
  # Group by both method and level to get granular metrics
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    Avg_Cov95 = mean(Cov95, na.rm = TRUE) * 100, # Convert to percentage
    .groups = "drop"
  ) %>%
  # Convert level to a factor so columns pivot in your preferred structural order
  dplyr::mutate(
    level = factor(level, levels = level_order),
    method = factor(method, levels = methods)
  ) %>%
  # Filter out any levels that didn't have data matching our factor levels
  dplyr::filter(!is.na(level)) %>%
  # Pivot from long to wide format (Levels become columns)
  tidyr::pivot_wider(
    names_from = level, 
    values_from = Avg_Cov95
  ) %>%
  # Ensure rows stay sorted in the original custom methods order
  dplyr::arrange(method) %>%
  select(method, Overall, Total, Purpose, StatePurpose, ZonePurpose, RegionPurpose)

# =============================================================================
# 7. OUTPUT LEVEL-SPECIFIC COVERAGE TABLE
# =============================================================================

cat("\n--- 95% EMPIRICAL COVERAGE BY HIERARCHICAL LEVEL (%) ---\n")
coverage_95_by_level %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 2))) %>% # Coverage is best read at 2 decimals
  print(n = Inf)

# Export to CSV
write.csv(coverage_95_by_level, 
          file.path(dir_out, "coverage_95_by_level_summary.csv"), 
          row.names = FALSE)