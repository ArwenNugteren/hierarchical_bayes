library(dplyr)
library(tidyr)
library(stringr)
library(tsutils)

# =============================================================================
# 0. CONFIGURATION & STRATEGIC PATHS
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
  "accuracy_misspecified_independent_v2",
  "accuracy_misspecified_independent_mint_v2",
  "accuracy_misspecified_hierarchical_v2",
  "accuracy_misspecified_hierarchical_mint_v2",
  "accuracy_misspecified_coherency1_v2",
  "accuracy_misspecified_coherency0.707_v2",
  "accuracy_misspecified_coherency0.707_topweighted_v2"
)

stopifnot(length(methods) == length(folders))

dir_out   <- "nemenyi_misspecified_results"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

horizon     <- 5L   
base_levels <- c("Overall", "top", "mid", "bottom")
level_labels <- c(
  "Overall" = "Overall System",
  "top"     = "Top Level",
  "mid"     = "Middle Level",
  "bottom"  = "Bottom Level"
)

# =============================================================================
# 1. LOAD GRANULAR WINDOW DATA (Preserves Series and Rolls/Timesteps)
# =============================================================================
load_granular_method <- function(folder, method_name) {
  acc_files <- list.files(folder, pattern = "^accuracy_(roll|timestep)_[0-9]+\\.rds$", full.names = TRUE)
  if (length(acc_files) == 0) {
    warning("No accuracy files found in: ", folder)
    return(NULL)
  }
  
  load_and_tag <- function(f) {
    window_num <- as.integer(stringr::str_extract(basename(f), "[0-9]+"))
    res        <- readRDS(f)
    uni        <- res$univariate %>% dplyr::mutate(window = window_num)
    uni
  }
  
  all_uni <- dplyr::bind_rows(lapply(acc_files, load_and_tag))
  
  uni_clean <- all_uni %>%
    dplyr::mutate(se = ifelse(is.finite(se), se, NA_real_),
                  crps = ifelse(is.finite(crps), crps, NA_real_)) %>%
    dplyr::select(level, series, horizon, window, se, crps) %>%
    dplyr::mutate(method = method_name)
  
  return(uni_clean)
}

# =============================================================================
# 2. INGEST DATA PIPELINE
# =============================================================================
all_data_list <- list()
for (i in seq_along(folders)) {
  message("Extracting metrics from: ", folders[i])
  res <- load_granular_method(folders[i], methods[i])
  if (!is.null(res)) {
    all_data_list[[i]] <- res
  }
}
all_raw_data <- dplyr::bind_rows(all_data_list)

if (nrow(all_raw_data) == 0) {
  stop("Critical Error: No data was loaded. Check folder names and RDS paths.")
}

# =============================================================================
# 3. CONSTRUCT AGGREGATIONS FOR THE RANK MATRIX
# =============================================================================
message("\n--- Running Node Aggregations ---")

level_window_metrics <- all_raw_data %>%
  dplyr::group_by(method, level, horizon, window) %>%
  dplyr::summarise(
    RMSE = sqrt(mean(se, na.rm = TRUE)), # Structural rolling RMSE
    CRPS = mean(crps, na.rm = TRUE),
    .groups = "drop"
  )

overall_window_metrics <- all_raw_data %>%
  dplyr::group_by(method, horizon, window) %>%
  dplyr::summarise(
    RMSE = sqrt(mean(se, na.rm = TRUE)),
    CRPS = mean(crps, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

master_window_data <- dplyr::bind_rows(level_window_metrics, overall_window_metrics)

# =============================================================================
# 4. NEMENYI PAGE GENERATOR FUNCTION
# =============================================================================
append_nemenyi_page <- function(df, target_level, metric_name, target_h = NULL) {
  df_sub <- df %>% dplyr::filter(level == target_level)
  
  if (is.null(target_h)) {
    df_sub <- df_sub %>% dplyr::mutate(block_id = paste0("h", horizon, "_w", window))
    title_suffix <- " (Horizons Pooled 1-5)"
  } else {
    df_sub <- df_sub %>% dplyr::filter(horizon == target_h) %>% dplyr::mutate(block_id = paste0("w", window))
    title_suffix <- paste0(" (Horizon ", target_h, ")")
  }
  
  matrix_prep <- df_sub %>%
    dplyr::select(block_id, method, value = dplyr::all_of(metric_name)) %>%
    tidyr::pivot_wider(names_from = method, values_from = value, values_fn = mean)
  
  mat <- as.matrix(matrix_prep[, -1, drop = FALSE])
  rownames(mat) <- matrix_prep$block_id
  mat <- mat[rowSums(!is.na(mat)) == ncol(mat), , drop = FALSE]
  
  if (nrow(mat) < 2) {
    message("⚠️ Insufficient rows for Nemenyi: Level ", target_level, " | Metric ", metric_name)
    return(NULL)
  }
  
  # Structural side-by-side plot layout on the active PDF device page
  old_par <- par(mfrow = c(1, 2), mar = c(5, 14, 4, 2), cex.axis = 0.9)
  
  # Left Plot: Critical Difference (MCB)
  tsutils::nemenyi(mat, conf.int = 0.95, plottype = "vmcb")
  title(main = paste0(metric_name, ": ", level_labels[target_level], title_suffix), line = 2.5, adj = 0)
  
  # Right Plot: P-Value Matrix Significance Heatmap
  par(mar = c(10, 11, 4, 2))
  tsutils::nemenyi(mat, conf.int = 0.95, plottype = "matrix")
  title(main = "Pairwise Significance Matrix")
  
  par(old_par)
}

# =============================================================================
# 5. EXECUTION PIPELINE - ONE GRAPHICS DEVICE PER METRIC
# =============================================================================
target_metrics <- c("RMSE", "CRPS")

for (metric in target_metrics) {
  pdf_path <- file.path(dir_out, paste0("nemenyi_comprehensive_matrix_", metric, ".pdf"))
  message("Opening target device for: ", pdf_path)
  
  # Open a multi-page PDF document layout (Landscape orientation optimizes side-by-side plots)
  pdf(pdf_path, width = 14, height = 7.5)
  
  for (lev in base_levels) {
    message("  Adding layouts for Level: ", lev)
    
    # Page Variant A: Horizons Pooled 1-5
    append_nemenyi_page(df = master_window_data, target_level = lev, metric_name = metric, target_h = NULL)
    
    # Page Variant B-F: Horizons Evaluated Separately
    for (h in seq_len(horizon)) {
      append_nemenyi_page(df = master_window_data, target_level = lev, metric_name = metric, target_h = h)
    }
  }
  
  dev.off()
  message("Completed device print export for metric: ", metric)
}

message("\n--- Structural Nemenyi Unified Multi-page PDFs Generated Successfully ---")