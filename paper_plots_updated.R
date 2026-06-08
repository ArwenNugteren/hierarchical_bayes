library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(stringr)
library(tsutils)   # nemenyi() for CD + blue p-value matrix

# =============================================================================
# 0.  CONFIGURATION & HORIZON EXPANSION
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

horizon         <- 12L
dir_out         <- "comparison_plots"
baseline_method <- "Independent Bayes"
alpha           <- 0.05

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Structural hierarchy tiers
base_levels <- c("Overall", "Total", "Purpose", "StatePurpose", "ZonePurpose", "RegionPurpose")
# Dynamic expansion to include horizons 1-12 as distinct evaluation groups
horizon_levels <- paste0("h", seq_len(horizon))
level_order    <- c(base_levels, horizon_levels)

level_labels <- c(
  "Overall"       = "Overall System",
  "Total"         = "Australia",
  "Purpose"       = "Purposes",
  "StatePurpose"  = "States",
  "ZonePurpose"   = "Zones",
  "RegionPurpose" = "Regions"
)

# Append programmatic titles for the standalone horizon groups to handle plot titles safely
for(h in seq_len(horizon)) {
  level_labels[paste0("h", h)] <- paste0("Horizon ", h)
}

# =============================================================================
# 1.  LOAD RAW DATA FUNCTION (Preserves Granular Windows For Correct Math Progression)
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
  
  # Extract long granular structures with timestamps intact
  uni_clean <- all_uni %>%
    dplyr::mutate(LogS = ifelse(is.finite(logS), logS, NA_real_)) %>%
    dplyr::select(level, series, horizon, timestep, ae, se, crps, LogS, covered_90, covered_95)
  
  # System-wide multivariate blocks; safely mean across configurations here
  multi_clean <- all_multi %>%
    dplyr::group_by(horizon, timestep) %>%
    dplyr::summarise(
      ES           = mean(es, na.rm = TRUE),
      MAI_rp_total = mean(abs(mai_rp_total), na.rm = TRUE),   
      JSD_rp_total = mean(jsd_rp_total, na.rm = TRUE),
      .groups      = "drop"
    )
  
  list(uni = uni_clean %>% dplyr::mutate(method = method_name), 
       multi = multi_clean %>% dplyr::mutate(method = method_name))
}

# =============================================================================
# 2.  DATA INGESTION PIPELINE
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
  stop("Critical Error: No data was loaded. Check folder names and RDS paths.")
}

# =============================================================================
# 3.  CRITICAL MATHEMATICAL INTERSECT: TIME-SERIES ROLLING WINDOW AGGREGATION
# =============================================================================

# Step A: Compute accurate summary metric vectors for each asset over all 121 timesteps
series_horizon_metrics <- all_uni_raw %>%
  dplyr::group_by(method, level, series, horizon) %>%
  dplyr::summarise(
    Series_MAE   = mean(ae, na.rm = TRUE),
    Series_RMSE  = sqrt(mean(se, na.rm = TRUE)), # Correct: Root-mean squared error across windows
    Series_CRPS  = mean(crps, na.rm = TRUE),
    Series_LogS  = mean(LogS, na.rm = TRUE),
    Series_Cov90 = mean(covered_90, na.rm = TRUE) * 100,
    Series_Cov95 = mean(covered_95, na.rm = TRUE) * 100,
    .groups      = "drop"
  )

# Step B: Average metrics across distinct individual nodes inside their tier levels
level_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, level, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), # Simple mean over asset nodes
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  )

# Step C: Formulate systemic "Overall" layer metrics by compressing the entire layout 
overall_horizon_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), # Mean across all system assets
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = "Overall")

# Step D: Construct standalone horizon macro frames (h1-h12 rows)
horizon_macro_summary <- series_horizon_metrics %>%
  dplyr::group_by(method, horizon) %>%
  dplyr::summarise(
    MAE   = mean(Series_MAE, na.rm = TRUE),
    RMSE  = mean(Series_RMSE, na.rm = TRUE), 
    CRPS  = mean(Series_CRPS, na.rm = TRUE),
    LogS  = mean(Series_LogS, na.rm = TRUE),
    Cov90 = mean(Series_Cov90, na.rm = TRUE),
    Cov95 = mean(Series_Cov95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(level = paste0("h", horizon))

# Stitch all elements to compose the unified plotting frame matching the original `all_data` structure
all_data <- dplyr::bind_rows(level_horizon_summary, overall_horizon_summary, horizon_macro_summary) %>%
  dplyr::left_join(
    all_multi_raw %>% 
      dplyr::group_by(method, horizon) %>% 
      dplyr::summarise(across(c(ES, MAI_rp_total, JSD_rp_total), ~ mean(.x, na.rm = TRUE)), .groups = "drop"),
    by = c("method", "horizon")
  ) %>%
  dplyr::mutate(method = factor(method, levels = methods))

# =============================================================================
# 4.  COMPUTE SUMMARY_LM (Averaged Over Horizons for Dynamic Summaries)
# =============================================================================

summary_lm <- all_data %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE          = mean(MAE, na.rm = TRUE), 
    RMSE         = mean(RMSE, na.rm = TRUE), # Linear mean across horizons of structural rolling RMSE
    CRPS         = mean(CRPS, na.rm = TRUE),
    Cov90        = mean(Cov90, na.rm = TRUE),
    Cov95        = mean(Cov95, na.rm = TRUE),
    MAI_rp_total = mean(MAI_rp_total, na.rm = TRUE),
    JSD_rp_total = mean(JSD_rp_total, na.rm = TRUE),
    ES           = mean(ES, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  dplyr::mutate(method = factor(method, levels = methods))

# =============================================================================
# 5.  SHARED HELPERS & COLOR DEFINITIONS
# =============================================================================

okabe_ito <- c(
  "Independent Bayes"                 = "#000000",
  "Independent Bayes MinT"              = "#E69F00",
  "Hierarchical Bayes"                  = "#A8D8EA",
  "Hierarchical Bayes MinT"             = "#00B957",
  "Hierarchical Bayes softer coherency" = "#F0E442",
  "Hierarchical Bayes tighter coherency"= "#0072B2",
  "Hierarchical Bayes total-weighted"   = "#E03000",
  "Hierarchical Bayes purposes-weighted"= "#CC79A7"
)

cb_shapes <- c(
  "Independent Bayes"                 = 16,  
  "Independent Bayes MinT"              = 17,  
  "Hierarchical Bayes"                  = 15,  
  "Hierarchical Bayes MinT"             = 18,  
  "Hierarchical Bayes softer coherency" = 3,   
  "Hierarchical Bayes tighter coherency"= 4,   
  "Hierarchical Bayes total-weighted"   = 8,   
  "Hierarchical Bayes purposes-weighted"= 25  
)

shared_fill_scale   <- function() scale_fill_manual(values = okabe_ito, name = NULL, drop = FALSE)
shared_colour_scale <- function() scale_color_manual(values = okabe_ito, name = NULL, drop = FALSE)

shared_legend <- get_legend(
  ggplot(all_data %>% dplyr::distinct(method),
         aes(x = method, y = 1, colour = method, shape = method)) +
    geom_point(size = 4.5, stroke = 1.2) + 
    shared_colour_scale() +
    scale_shape_manual(values = cb_shapes, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                 label.theme = element_text(size = 11))) +
    theme_minimal() +
    theme(
      legend.position = "bottom", 
      legend.key.size = unit(0.8, "cm"),
      legend.text     = element_text(margin = margin(r = 10, l = 2))
    )
)

# =============================================================================
# 6.  METRIC SPECIFICATION (Uncomment the ones you wish to export to PDF plots)
# =============================================================================

metrics <- list(
  list(col = "MAE",          label = "MAE",              lower_better = TRUE,
       fix_floor = FALSE, coverage_target = NULL,
       levels_with_data = level_order,   filename = "MAE"),
  list(col = "RMSE",         label = "RMSE",             lower_better = TRUE,
       fix_floor = FALSE, coverage_target = NULL,
       levels_with_data = level_order,   filename = "RMSE"),
  list(col = "CRPS",         label = "CRPS",             lower_better = TRUE,
       fix_floor = FALSE, coverage_target = NULL,
       levels_with_data = level_order,   filename = "CRPS")
)

# =============================================================================
# 7.  BAR PANEL BUILDER  
# =============================================================================

make_bar_panel <- function(data, metric_col, metric_label, lev,
                           fix_floor, coverage_target = NULL) {
  
  df <- data %>%
    dplyr::filter(level == lev, !is.na(.data[[metric_col]])) %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(val = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::rename(!!metric_col := val)
  
  if (nrow(df) == 0)
    return(ggplot() + theme_void() + labs(title = paste(lev, "- no data")))
  
  vals    <- df[[metric_col]]
  y_min   <- min(vals, na.rm = TRUE)
  y_max   <- max(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (fix_floor) {
    y_lo <- 0; y_hi <- y_max * 1.08
  } else if (y_range < 1e-10) {
    y_lo <- 0; y_hi <- y_max * 1.08
  } else {
    y_lo <- max(2 * y_min - y_max, y_min - 0.10 * abs(y_min))
    y_hi <- y_max + 0.05 * y_range  # <-- Fixed typo from ':=' to '<-' here
  }
  
  ggplot(df, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    shared_fill_scale() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    labs(title = level_labels[lev], x = NULL, y = metric_label) +
    theme_bw(base_size = 13) +
    theme(legend.position    = "none",
          plot.title         = element_text(size = 13, face = "bold"),
          axis.text.x        = element_text(size = 11, angle = 35, hjust = 1),
          axis.text.y        = element_text(size = 11),
          axis.title.y       = element_text(size = 12),
          panel.grid.major.x = element_blank())
}

# =============================================================================
# 8.  AUDITED NEMENYI HELPERS (Constructs Matrix Slices for Nonparametric Testing)
# =============================================================================

build_nemenyi_matrix <- function(df, metric_col, lower_better) {
  df <- df %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::mutate(method = stringr::str_wrap(method, width = 22))
  
  # Map directly across available observations
  df <- df %>% 
    dplyr::mutate(obs_id = paste("h", horizon, sep = "_"))
  
  wide <- df %>%
    dplyr::select(obs_id, method, value = dplyr::all_of(metric_col)) %>%
    tidyr::pivot_wider(names_from  = method,
                       values_from = value,
                       values_fn   = mean) 
  
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$obs_id
  
  mat <- mat[rowSums(!is.na(mat)) == ncol(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  
  if (!lower_better) mat <- -mat
  mat
}

capture_nemenyi_plot <- function(mat, alpha_val = 0.05, res = 150) {
  if (is.null(mat)) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for Nemenyi test", size = 5, hjust = 0.5))
  }
  
  tmp_mcb <- tempfile(fileext = ".png")
  png(tmp_mcb, width = 1000, height = 1400, res = res)
  old_par_mcb <- par(mar = c(5, 12, 4, 1), cex.axis = 1.0)
  tryCatch({ tsutils::nemenyi(mat, conf.int = 1 - alpha_val, plottype = "vmcb") }, 
           error = function(e) { plot.new(); title(paste("MCB Error:", conditionMessage(e))) })
  par(old_par_mcb); dev.off()
  
  tmp_mat <- tempfile(fileext = ".png")
  png(tmp_mat, width = 1400, height = 1400, res = res)
  old_par_mat <- par(mar = c(10, 11, 4, 2), pty = "s", cex.axis = 1.0)
  tryCatch({ tsutils::nemenyi(mat, conf.int = 1 - alpha_val, plottype = "matrix") }, 
           error = function(e) { plot.new(); title(paste("Matrix Error:", conditionMessage(e))) })
  par(old_par_mat); dev.off()
  
  g_mcb <- ggplot() + annotation_raster(png::readPNG(tmp_mcb), xmin=0, xmax=1, ymin=0, ymax=1) + coord_fixed() + theme_void()
  g_mat <- ggplot() + annotation_raster(png::readPNG(tmp_mat), xmin=0, xmax=1, ymin=0, ymax=1) + coord_fixed() + theme_void()
  
  cowplot::plot_grid(g_mcb, g_mat, ncol = 2, rel_widths = c(0.38, 0.62))
}

make_nemenyi_page <- function(df, metric_col, lower_better, title_str, res = 150) {
  mat <- build_nemenyi_matrix(df, metric_col, lower_better)
  if (is.null(mat)) {
    return(list(grob = ggplot() + theme_void() + labs(title = paste(title_str, "- insufficient data"))))
  }
  grob <- capture_nemenyi_plot(mat, alpha_val = alpha, res = res)
  titled <- plot_grid(
    ggdraw() + draw_label(title_str, fontface = "bold", size = 14, x = 0.5, hjust = 0.5),
    grob, ncol = 1, rel_heights = c(0.05, 0.95)
  )
  list(grob = titled)
}

# =============================================================================
# 9.  RELATIVE LINE-PLOT
# =============================================================================

make_rel_lineplot <- function(df_avg, metric_col, metric_label, lower_better, coverage_target = NULL) {
  
  # Fetch baseline values for scaling
  baseline_vals <- df_avg %>%
    dplyr::filter(method == baseline_method) %>%
    dplyr::select(level, baseline = dplyr::all_of(metric_col))
  
  rel_df <- df_avg %>%
    dplyr::filter(
      method != baseline_method, 
      !is.na(.data[[metric_col]]),
      level %in% base_levels 
    ) %>%
    dplyr::left_join(baseline_vals, by = "level") %>%
    dplyr::filter(!is.na(baseline), abs(baseline) > .Machine$double.eps) %>%
    dplyr::mutate(
      rel_value = .data[[metric_col]] / baseline,
      level     = factor(level, levels = base_levels) 
    )
  
  baseline_df <- df_avg %>%
    dplyr::filter(
      method == baseline_method, 
      !is.na(.data[[metric_col]]),
      level %in% base_levels
    ) %>%
    dplyr::mutate(
      rel_value = 1, 
      level     = factor(level, levels = base_levels)
    )
  
  if (nrow(rel_df) == 0) return(ggplot() + theme_void() + labs(title = "No relative data"))
  
  direction_note <- if (lower_better) "(< 1 is better)" else "(> 1 is better)"
  
  ggplot(mapping = aes(x = level, y = rel_value, colour = method, group = method, shape = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_line(data  = rel_df,      linewidth = 0.9) +
    geom_line(data  = baseline_df, linewidth = 0.9) +
    geom_point(data = rel_df,      size = 3.0) + 
    geom_point(data = baseline_df, size = 3.0) +
    shared_colour_scale() +
    scale_shape_manual(values = cb_shapes, name = NULL, drop = FALSE) +
    scale_x_discrete(limits = base_levels, 
                     labels = function(x) {
                       clean_names <- ifelse(x %in% names(level_labels), level_labels[x], x)
                       str_wrap(clean_names, width = 12)
                     }) +
    labs(
      title    = paste(metric_label, "relative to", baseline_method), 
      subtitle = direction_note,
      x        = "Structural Tiers", 
      y        = paste0("Ratio vs. ", baseline_method)
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position    = "none", 
      plot.title         = element_text(size = 12, face = "bold"),
      axis.text.x        = element_text(size = 11, angle = 0, hjust = 0.5), 
      panel.grid.major.x = element_blank()
    )
}

# =============================================================================
# 10. MAIN EXECUTION LOOP
# =============================================================================

for (m in metrics) {
  avail <- all_data %>% filter(!is.na(.data[[m$col]])) %>% nrow()
  if (avail == 0) { message("Skipping ", m$filename, " - no data"); next }
  
  direction_note    = if (m$lower_better) "(lower is better)" else "(higher is better)"
  applicable_levels = intersect(level_order, as.character(unique(all_data$level[!is.na(all_data[[m$col]])])))
  
  pdf_path <- file.path(dir_out, paste0(m$filename, ".pdf"))
  message("Building ", pdf_path, " ...")
  pdf(pdf_path, width = 15, height = 9) 
  
  # ── Pages 1...H: bar panels per horizon ──────────────────────────────────
  for (h in seq_len(horizon)) {
    df_h   <- all_data %>% filter(horizon == h)
    panels <- lapply(applicable_levels, function(lev) make_bar_panel(df_h, m$col, m$label, lev, m$fix_floor, m$coverage_target))
    if (length(panels) == 0) next
    print(plot_grid(ggplotGrob(wrap_plots(panels, ncol = 3)), shared_legend, ncol = 1, rel_heights = c(0.85, 0.15)))
  }
  
  # ── Page H+1: averaged over horizons, bar panels ─────────────────────────
  df_avg_level <- summary_lm %>% filter(level %in% applicable_levels)
  panels_avg   <- lapply(applicable_levels, function(lev) make_bar_panel(df_avg_level, m$col, m$label, lev, m$fix_floor, m$coverage_target))
  if (length(panels_avg) > 0) {
    print(plot_grid(ggplotGrob(wrap_plots(panels_avg, ncol = 3)), shared_legend, ncol = 1, rel_heights = c(0.85, 0.15)))
  }
  
  # ── Nemenyi Framework: Evaluated for Macro Tiers ──────────────────────────
  target_nemenyi_levels <- c("Overall", "Total", "Purpose")
  for (lev in target_nemenyi_levels) {
    if (lev %in% applicable_levels) {
      df_lev <- all_data %>% dplyr::filter(level == lev)
      pg_all <- make_nemenyi_page(df_lev, m$col, m$lower_better, title_str = sprintf("Nemenyi Framework: %s | Level: %s | All Horizons Pooled", m$label, lev))
      print(pg_all$grob)
    }
  }
  
  # ── Final page: relative line plot ────────────────────────────────────────
  rel_plot <- make_rel_lineplot(summary_lm, m$col, m$label, m$lower_better, m$coverage_target)
  print(plot_grid(rel_plot, shared_legend, ncol = 1, rel_heights = c(0.80, 0.20)))
  
  dev.off()
  message("  Saved: ", pdf_path)
}