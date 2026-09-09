library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(tsutils)   # For nemenyi() matrix & CD graphs
library(png)

# =============================================================================
# 0. CONFIGURATION & EXACT DIRECTORY PATHS
# =============================================================================

folders <- c(
  "accuracy_purposes_independent_v3/post_hoc_es_base",
  "accuracy_purposes_independent_MinT_v3/post_hoc_es_mint",
  "accuracy_purposes_hierarchy_v2/post_hoc_es_base",
  "accuracy_purposes_hierarchy_MinT_v2/post_hoc_es_mint",
  "accuracy_purposes_hierarchy_coherency0.005v3/es_post_hoc",
  "accuracy_purposes_hierarchy_coherency0.002v2/es_post_hoc",
  "accuracy_purposes_hierarchy_coherency0.002_totalweight_nrp/es_post_hoc",
  "accuracy_purposes_hierarchy_coherency0.002_purposesweight_nregions/es_post_hoc"
)

method_labels <- c(
  "Independent Bayes",
  "Independent Bayes MinT",
  "Hierarchical Bayes",
  "Hierarchical Bayes MinT",
  "Hierarchical Bayes softer coherency",  
  "Hierarchical Bayes tighter coherency", 
  "Hierarchical Bayes total-weighted",
  "Hierarchical Bayes purposes-weighted"
)

folder_mapping <- setNames(method_labels, folders)

dir_out <- "es_diagnostic_results"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

horizon   <- 12L
alpha     <- 0.05
plot_res  <- 150

# =============================================================================
# 1. UNIFORM DATA INGESTION LOADER
# =============================================================================

load_folder_es <- function(folder_path, display_name) {
  message("Ingesting Energy Score from directory: ", folder_path)
  
  post_hoc_files  <- list.files(folder_path, pattern = "^es_(base|mint)_origin_[0-9]+\\.rds$", full.names = TRUE)
  accuracy_files  <- list.files(folder_path, pattern = "^accuracy_timestep_[0-9]+\\.rds$", full.names = TRUE)
  
  if (length(post_hoc_files) > 0) {
    df <- dplyr::bind_rows(lapply(post_hoc_files, readRDS))
    if ("energy_score" %in% names(df)) {
      df <- df %>% dplyr::rename(es = energy_score)
    }
  } else if (length(accuracy_files) > 0) {
    df <- dplyr::bind_rows(lapply(accuracy_files, function(f) {
      res <- readRDS(f)
      if (is.list(res) && "multivariate" %in% names(res)) {
        return(res$multivariate %>% dplyr::select(origin, horizon, es))
      } else if ("es" %in% names(res)) {
        return(res %>% dplyr::select(origin, horizon, es))
      }
      NULL
    }))
  } else {
    warning("No valid RDS files located in target folder: ", folder_path)
    return(NULL)
  }
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>% 
    dplyr::mutate(method = display_name) %>% 
    dplyr::select(method, origin, horizon, es)
}

# =============================================================================
# 2. RUN PIPELINE INGESTION & MASTER FRAMING
# =============================================================================

raw_es_list <- list()
for (f in folders) {
  lbl <- folder_mapping[[f]]
  res <- load_folder_es(f, lbl)
  if (!is.null(res)) raw_es_list[[f]] <- res
}

master_es_data <- dplyr::bind_rows(raw_es_list) %>%
  dplyr::mutate(method = factor(method, levels = method_labels))

if (nrow(master_es_data) == 0) {
  stop("Execution halted: No valid Energy Score data could be aggregated from the specified folders.")
}

# =============================================================================
# 3. GENERATE COMPREHENSIVE HORIZON SUMMARY TABLE (EXPORTS CSV)
# =============================================================================

message("\nComputing comprehensive horizon summaries...")

horizon_summary_table <- master_es_data %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(mean_es = mean(es, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(horizon_lbl = paste0("h", horizon)) %>%
  dplyr::select(-horizon) %>%
  tidyr::pivot_wider(names_from = horizon_lbl, values_from = mean_es)

macro_averages <- master_es_data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    `h1-6 Avg`  = mean(es[horizon %in% 1:6], na.rm = TRUE),
    `h1-12 Avg` = mean(es[horizon %in% 1:12], na.rm = TRUE),
    .groups     = "drop"
  )

final_summary_table <- horizon_summary_table %>%
  dplyr::left_join(macro_averages, by = "method") %>%
  dplyr::arrange(method)

table_out_path <- file.path(dir_out, "es_comprehensive_summary_table.csv")
write.csv(final_summary_table, table_out_path, row.names = FALSE)
message("Saved summary layout spreadsheet to: ", table_out_path)

print(as.data.frame(final_summary_table), digits = 4)

# =============================================================================
# 4. REVISED NONPARAMETRIC MATRIX GENERATOR (WITH STRING WRAPPING)
# =============================================================================

message("\nBuilding Nonparametric Nemenyi Test Matrices (Averaged Over Horizons)...")

collapsed_es_data <- master_es_data %>%
  dplyr::group_by(method, origin) %>%
  dplyr::summarise(es_horizon_avg = mean(es, na.rm = TRUE), .groups = "drop") %>%
  # FIX: Apply clean string wrapping for long method labels to prevent canvas clipping
  dplyr::mutate(method = stringr::str_wrap(method, width = 22))

# Extract unique wrapped names in their designated order to preserve matrix tracking
wrapped_method_levels <- stringr::str_wrap(method_labels, width = 22)

wide_nemenyi <- collapsed_es_data %>%
  dplyr::mutate(method = factor(method, levels = wrapped_method_levels)) %>%
  dplyr::select(origin, method, es_horizon_avg) %>%
  tidyr::pivot_wider(names_from = method, values_from = es_horizon_avg)

matrix_splicer <- as.matrix(wide_nemenyi[, -1, drop = FALSE])
rownames(matrix_splicer) <- paste0("origin_", wide_nemenyi$origin)

matrix_splicer <- matrix_splicer[rowSums(!is.na(matrix_splicer)) == ncol(matrix_splicer), , drop = FALSE]

if (nrow(matrix_splicer) < 2) {
  stop("Matrix slicing error: Insufficient paired observations remain across all trials.")
}

# =============================================================================
# 5. RENDER DIAGNOSTIC GRAPHICS CANVAS PAGE
# =============================================================================

pdf_out_path <- file.path(dir_out, "energy_score_nemenyi_test.pdf")
message("Compiling Nemenyi framework profile page: ", pdf_out_path)

tmp_mcb <- tempfile(fileext = ".png")
tmp_mat <- tempfile(fileext = ".png")

png(tmp_mcb, width = 1000, height = 1400, res = plot_res)
par(mar = c(5, 12, 4, 1), cex.axis = 1.0)
tsutils::nemenyi(matrix_splicer, conf.int = 1 - alpha, plottype = "vmcb")
dev.off()

png(tmp_mat, width = 1400, height = 1400, res = plot_res)
par(mar = c(10, 11, 4, 2), pty = "s", cex.axis = 1.0)
tsutils::nemenyi(matrix_splicer, conf.int = 1 - alpha, plottype = "matrix")
dev.off()

g_mcb <- ggplot() + annotation_raster(png::readPNG(tmp_mcb), xmin=0, xmax=1, ymin=0, ymax=1) + coord_fixed() + theme_void()
g_mat <- ggplot() + annotation_raster(png::readPNG(tmp_mat), xmin=0, xmax=1, ymin=0, ymax=1) + coord_fixed() + theme_void()

grob_panel <- cowplot::plot_grid(g_mcb, g_mat, ncol = 2, rel_widths = c(0.40, 0.60))

title_header <- ggdraw() + 
  draw_label("Nemenyi Framework Matrix: System-Wide Overall Energy Score (ES) | All Horizons Pooled", 
             fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

final_grob_page <- plot_grid(title_header, grob_panel, ncol = 1, rel_heights = c(0.05, 0.95))

pdf(pdf_out_path, width = 15, height = 9)
print(final_grob_page)
dev.off()

unlink(c(tmp_mcb, tmp_mat))
message("Execution engine completed diagnostic tracking run successfully.")