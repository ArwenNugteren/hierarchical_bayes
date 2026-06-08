library(dplyr)
library(purrr)
library(tidyr)

# =============================================================================
# 1. SETTINGS & DIRECTORIES
# =============================================================================
job_name <- "accuracy_misspecified_coherency0.707"

dir_output_base <- file.path(job_name, "post_hoc_es_base")
dir_output_mint <- file.path(job_name, "post_hoc_es_mint")

# Find all successfully saved RDS result files
base_files <- list.files(dir_output_base, pattern = "es_base_roll_.*\\.rds$", full.names = TRUE)
mint_files <- list.files(dir_output_mint, pattern = "es_mint_roll_.*\\.rds$", full.names = TRUE)

if (length(base_files) == 0 || length(mint_files) == 0) {
  stop("No metric files found. Ensure your rolling window processing loop has run successfully.")
}

# =============================================================================
# 2. READ AND AGGREGATE PERFORMANCE METRICS
# =============================================================================
cat(sprintf("Found %d base files and %d MinT files. Aggregating...\n", length(base_files), length(mint_files)))

# Efficiently bind all window frames together
df_base <- map_dfr(base_files, readRDS) %>% rename(es_base = energy_score)
df_mint <- map_dfr(mint_files, readRDS) %>% rename(es_mint = energy_score)

# Join tables by window (roll) and horizon step
df_joint <- inner_join(df_base, df_mint, by = c("roll", "horizon"))

# =============================================================================
# 3. COMPUTE HORIZON-SPECIFIC SUMMARY STATISTICS
# =============================================================================
summary_table <- df_joint %>%
  group_by(horizon) %>%
  summarise(
    n_windows = n(),
    mean_es_base = mean(es_base, na.rm = TRUE),
    mean_es_mint = mean(es_mint, na.rm = TRUE),
    pct_improvement = (mean_es_base - mean_es_mint) / mean_es_base * 100,
    .groups = "drop"
  ) %>%
  # Convert horizon to character so it can bind with the summary label
  mutate(horizon = as.character(horizon))

# Add a global summary row tracking overall performance across all horizons
overall_summary <- df_joint %>%
  summarise(
    horizon = "Overall (Average)",
    n_windows = n(),
    mean_es_base = mean(es_base, na.rm = TRUE),
    mean_es_mint = mean(es_mint, na.rm = TRUE),
    pct_improvement = (mean_es_base - mean_es_mint) / mean_es_base * 100
  )

final_table <- bind_rows(summary_table, overall_summary)
write.csv(final_table, file = "ES_summary")

