library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(stringr)
library(tsutils)   # nemenyi() for CD + blue p-value matrix

# =============================================================================
# 0.  CONFIGURATION
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

level_order <- c("Overall", "Total", "Purpose",
                 "StatePurpose", "ZonePurpose", "RegionPurpose")

# =============================================================================
# 1.  LOAD RAW DATA FUNCTION
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
  
  # ── Per-level rows (one row per timestep × horizon × series) ──────────────
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
# 1b. INGEST DATA FROM ALL FOLDERS (Fixes "all_data" missing error)
# =============================================================================

all_data_list <- list()
for (i in seq_along(folders)) {
  message("Processing folder: ", folders[i])
  all_data_list[[i]] <- load_method(folders[i], methods[i])
}

all_data <- dplyr::bind_rows(all_data_list)

if (nrow(all_data) == 0) {
  stop("Critical Error: No data was loaded. Check folder names and RDS paths.")
}

# Enforce your custom structural order for the legend and plots
all_data <- all_data %>%
  dplyr::mutate(method = factor(method, levels = methods))

# =============================================================================
# 1c. COMPUTE SUMMARY_LM (Option A: Simple Linear Mean to Match Line Plot)
# =============================================================================

# =============================================================================
# 1c. COMPUTE SUMMARY_LM (Option B: True Root of Mean for RMSE + MAE Added)
# =============================================================================

summary_lm <- all_data %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise(
    MAE          = mean(MAE, na.rm = TRUE), 
    RMSE         = sqrt(mean(RMSE^2, na.rm = TRUE)), 
    CRPS         = mean(CRPS, na.rm = TRUE),
    Cov90        = mean(Cov90, na.rm = TRUE),
    Cov95        = mean(Cov95, na.rm = TRUE),
    MAI_rp_total = mean(MAI_rp_total, na.rm = TRUE),
    JSD_rp_total = mean(JSD_rp_total, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  # Preserve factor sorting in the summary dataframe as well
  dplyr::mutate(method = factor(method, levels = methods))

# =============================================================================
# 2.  SHARED HELPERS  — Okabe-Ito colour scale
# =============================================================================

okabe_ito <- c(
  "Independent Bayes"                   = "#000000",  # black
  "Independent Bayes MinT"              = "#E69F00",  # orange
  "Hierarchical Bayes"                  = "#A8D8EA",  # lighter sky blue
  "Hierarchical Bayes MinT"             = "#00B957",  # more vivid green
  "Hierarchical Bayes softer coherency" = "#F0E442",  # yellow
  "Hierarchical Bayes tighter coherency"= "#0072B2",  # blue
  "Hierarchical Bayes total-weighted"   = "#E03000",  # redder vermillion
  "Hierarchical Bayes purposes-weighted"= "#CC79A7"   # reddish purple
)

shared_fill_scale <- function()
  scale_fill_manual(values = okabe_ito, name = NULL, drop = FALSE)

shared_colour_scale <- function()
  scale_color_manual(values = okabe_ito, name = NULL, drop = FALSE)

shared_legend <- get_legend(
  ggplot(all_data %>% dplyr::distinct(method),
         aes(x = method, y = 1, fill = method)) +
    geom_col() +
    shared_fill_scale() +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                               label.theme = element_text(size = 11))) +
    theme(legend.position = "bottom", legend.key.size = unit(0.6, "cm"))
)

# =============================================================================
# 3.  METRIC SPEC
# =============================================================================

metrics <- list(
  # list(col = "MAE",          label = "MAE",              lower_better = TRUE,
  #      fix_floor = FALSE, coverage_target = NULL,
  #      levels_with_data = level_order,   filename = "MAE"),
  # list(col = "RMSE",         label = "RMSE",             lower_better = TRUE,
  #      fix_floor = FALSE, coverage_target = NULL,
  #      levels_with_data = level_order,   filename = "RMSE")
  list(col = "CRPS",         label = "CRPS",             lower_better = TRUE,
       fix_floor = FALSE, coverage_target = NULL,
       levels_with_data = level_order,   filename = "CRPS")
  # list(col = "LogS",         label = "Log Score",        lower_better = TRUE,
  #      fix_floor = FALSE, coverage_target = NULL,
  #      levels_with_data = level_order,   filename = "LogScore"),
  # list(col = "Cov90",        label = "90% Coverage",     lower_better = FALSE,
  #      fix_floor = FALSE, coverage_target = 0.90,
  #      levels_with_data = level_order,   filename = "Coverage90"),
  # list(col = "Cov95",        label = "95% Coverage",     lower_better = FALSE,
  #      fix_floor = FALSE, coverage_target = 0.95,
  #      levels_with_data = level_order,   filename = "Coverage95"),
  # list(col = "ES",           label = "Energy Score",     lower_better = TRUE,
  #      fix_floor = FALSE, coverage_target = NULL,
  #      levels_with_data = "Overall",     filename = "EnergyScore"),
  # list(col = "MAI_rp_total", label = "MAI (rp->total)",  lower_better = TRUE,
  #      fix_floor = TRUE,  coverage_target = NULL,
  #      levels_with_data = "Overall",     filename = "MAI"),
  # list(col = "JSD_rp_total", label = "JSD (rp->total)",  lower_better = TRUE,
  #      fix_floor = TRUE,  coverage_target = NULL,
  #      levels_with_data = "Overall",     filename = "JSD")
)

# =============================================================================
# 4.  BAR PANEL BUILDER  
# =============================================================================

make_bar_panel <- function(data, metric_col, metric_label, lev,
                           fix_floor, coverage_target = NULL) {
  
  df <- data %>%
    dplyr::filter(level == lev, !is.na(.data[[metric_col]])) %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      val = mean(.data[[metric_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!metric_col := val)
  
  if (nrow(df) == 0)
    return(ggplot() + theme_void() + labs(title = paste(lev, "- no data")))
  
  vals    <- df[[metric_col]]
  y_min   <- min(vals, na.rm = TRUE)
  y_max   <- max(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (fix_floor) {
    y_lo <- 0
    y_hi <- y_max * 1.08
  } else if (y_range < 1e-10) {
    y_lo <- 0
    y_hi <- y_max * 1.08
  } else {
    y_lo <- max(2 * y_min - y_max, y_min - 0.10 * abs(y_min))
    y_hi <- y_max + 0.05 * y_range
  }
  
  p <- ggplot(df, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    shared_fill_scale() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    labs(title = lev, x = NULL, y = metric_label) +
    theme_bw(base_size = 13) +
    theme(legend.position    = "none",
          plot.title         = element_text(size = 13, face = "bold"),
          axis.text.x        = element_text(size = 11, angle = 35, hjust = 1),
          axis.text.y        = element_text(size = 11),
          axis.title.y       = element_text(size = 12),
          panel.grid.major.x = element_blank())
  
  if (!is.null(coverage_target))
    p <- p +
    geom_hline(yintercept = coverage_target, linetype = "dashed",
               colour = "black", linewidth = 0.7) +
    annotate("text", x = Inf, y = coverage_target,
             label = paste0(" ", coverage_target * 100, "%"),
             hjust = 1.05, vjust = -0.4, size = 3.5, colour = "black")
  p
}

# =============================================================================
# 5.  NEMENYI HELPERS 
# =============================================================================

build_nemenyi_matrix <- function(df, metric_col, lower_better) {
  df <- df %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::mutate(method = stringr::str_wrap(method, width = 22))
  
  df <- df %>% 
    dplyr::group_by(method, level, horizon) %>% 
    dplyr::mutate(obs_idx = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(obs_id = paste(level, "h", horizon, "ts", timestep, "idx", obs_idx, sep = "_"))
  
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
             annotate("text", x = 0.5, y = 0.5,
                      label = "Insufficient data for Nemenyi test",
                      size = 5, hjust = 0.5))
  }
  
  tmp_mcb <- tempfile(fileext = ".png")
  png(tmp_mcb, width = 1000, height = 1400, res = res)
  old_par_mcb <- par(mar = c(5, 12, 4, 1), cex.axis = 1.0)
  
  tryCatch({
    tsutils::nemenyi(mat, conf.int = 1 - alpha_val, plottype = "vmcb")
  }, error = function(e) {
    plot.new()
    title(paste("MCB Plot error:", conditionMessage(e)))
  })
  par(old_par_mcb)
  dev.off()
  
  tmp_mat <- tempfile(fileext = ".png")
  png(tmp_mat, width = 1400, height = 1400, res = res)
  old_par_mat <- par(mar = c(10, 11, 4, 2), pty = "s", cex.axis = 1.0)
  
  tryCatch({
    tsutils::nemenyi(mat, conf.int = 1 - alpha_val, plottype = "matrix")
  }, error = function(e) {
    plot.new()
    title(paste("Matrix Plot error:", conditionMessage(e)))
  })
  par(old_par_mat)
  dev.off()
  
  img_mcb <- png::readPNG(tmp_mcb)
  g_mcb <- ggplot() +
    annotation_raster(img_mcb, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
  
  img_mat <- png::readPNG(tmp_mat)
  g_mat <- ggplot() +
    annotation_raster(img_mat, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
  
  cowplot::plot_grid(g_mcb, g_mat, ncol = 2, rel_widths = c(0.38, 0.62))
}

make_nemenyi_page <- function(df, metric_col, lower_better, title_str, res = 150) {
  mat <- build_nemenyi_matrix(df, metric_col, lower_better)
  
  if (is.null(mat)) {
    g <- ggplot() + theme_void() +
      labs(title = paste(title_str, "- insufficient data"))
    return(list(grob = g, result = NULL))
  }
  
  grob <- capture_nemenyi_plot(mat, alpha_val = alpha, res = res)
  
  titled <- plot_grid(
    ggdraw() + draw_label(title_str, fontface = "bold", size = 14, x = 0.5, hjust = 0.5),
    grob,
    ncol = 1, rel_heights = c(0.05, 0.95)
  )
  
  list(grob = titled)
}

# =============================================================================
# 6.  RELATIVE LINE-PLOT  
# =============================================================================

make_rel_lineplot <- function(df_avg, metric_col, metric_label, lower_better,
                              coverage_target = NULL) {
  
  baseline_vals <- df_avg %>%
    dplyr::filter(method == baseline_method) %>%
    dplyr::select(level, baseline = dplyr::all_of(metric_col))
  
  rel_df <- df_avg %>%
    dplyr::filter(method != baseline_method,
                  !is.na(.data[[metric_col]])) %>%
    dplyr::left_join(baseline_vals, by = "level") %>%
    dplyr::filter(!is.na(baseline), abs(baseline) > .Machine$double.eps) %>%
    dplyr::mutate(rel_value = .data[[metric_col]] / baseline,
                  level     = factor(level, levels = level_order))
  
  baseline_df <- df_avg %>%
    dplyr::filter(method == baseline_method,
                  !is.na(.data[[metric_col]])) %>%
    dplyr::mutate(rel_value = 1,
                  level     = factor(level, levels = level_order))
  
  if (nrow(rel_df) == 0)
    return(ggplot() + theme_void() +
             labs(title = paste(metric_label, "- no relative data")))
  
  direction_note <- if (lower_better) "(< 1 is better)" else "(> 1 is better)"
  
  ggplot(mapping = aes(x = level, y = rel_value,
                       colour = method, group = method)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey60", linewidth = 0.5) +
    geom_line(data  = rel_df,     linewidth = 0.9) +
    geom_point(data = rel_df,     size = 2.5) +
    geom_line(data  = baseline_df, linewidth = 0.9) +
    geom_point(data = baseline_df, size = 2.5) +
    shared_colour_scale() +
    scale_x_discrete(limits = level_order,
                     labels = function(x) str_wrap(x, width = 12)) +
    labs(
      title    = paste(metric_label, "relative to", baseline_method),
      subtitle = direction_note,
      x        = "Level  (coarsest -> finest)",
      y        = paste0("Ratio vs. ", baseline_method)
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position    = "none",
          plot.title         = element_text(size = 12, face = "bold"),
          plot.subtitle      = element_text(size = 10, colour = "grey40"),
          axis.text.x        = element_text(size = 10, angle = 25, hjust = 1),
          panel.grid.major.x = element_blank())
}

# =============================================================================
# 7.  MAIN LOOP  
# =============================================================================

for (m in metrics) {
  
  avail <- all_data %>% filter(!is.na(.data[[m$col]])) %>% nrow()
  if (avail == 0) { message("Skipping ", m$filename, " - no data"); next }
  
  direction_note    <- if (m$lower_better) "(lower is better)" else "(higher is better)"
  applicable_levels <- intersect(
    if (identical(m$levels_with_data, level_order)) level_order    else as.character(m$levels_with_data),
    as.character(unique(all_data$level[!is.na(all_data[[m$col]])]))
  )
  
  nemenyi_per_horizon <- length(applicable_levels) > 1
  
  pdf_path <- file.path(dir_out, paste0(m$filename, ".pdf"))
  message("Building ", pdf_path, " ...")
  pdf(pdf_path, width = 12, height = 8)
  
  # ── Pages 1...H: bar panels per horizon ──────────────────────────────────
  for (h in seq_len(horizon)) {
    df_h   <- all_data %>% filter(horizon == h)
    panels <- lapply(applicable_levels, function(lev)
      make_bar_panel(df_h, m$col, m$label, lev, m$fix_floor, m$coverage_target))
    if (length(panels) == 0) next
    n_cols   <- min(length(panels), 3L)
    combined <- wrap_plots(panels, ncol = n_cols) +
      plot_annotation(
        title = sprintf("%s  %s  -  Horizon h = %d", m$label, direction_note, h),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    print(plot_grid(ggplotGrob(combined), shared_legend,
                    ncol = 1, rel_heights = c(0.82, 0.18)))
  }
  
  # ── Page H+1: averaged over horizons, bar panels ─────────────────────────
  df_avg_level <- summary_lm %>% filter(level %in% applicable_levels)
  panels_avg   <- lapply(applicable_levels, function(lev)
    make_bar_panel(df_avg_level, m$col, m$label, lev, m$fix_floor, m$coverage_target))
  if (length(panels_avg) > 0) {
    n_cols       <- min(length(panels_avg), 3L)
    combined_avg <- wrap_plots(panels_avg, ncol = n_cols) +
      plot_annotation(
        title = sprintf("%s  %s  -  Averaged over all horizons", m$label, direction_note),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    print(plot_grid(ggplotGrob(combined_avg), shared_legend,
                    ncol = 1, rel_heights = c(0.82, 0.18)))
  }
  
  # ── Nemenyi Framework: Evaluated for Overall, Total, and Purpose Tiers ─────
  target_nemenyi_levels <- c("Overall", "Total", "Purpose")
  
  for (lev in target_nemenyi_levels) {
    if (lev %in% applicable_levels) {
      df_lev <- all_data %>% dplyr::filter(level == lev)
      
      # 1. Pooled across all horizons for this level
      pg_all <- make_nemenyi_page(
        df_lev, m$col, m$lower_better,
        title_str = sprintf("Nemenyi Framework: %s | Level: %s | All Horizons Pooled", m$label, lev)
      )
      print(pg_all$grob)
      
      # 2. Individual horizon breakdown for this level
      if (nemenyi_per_horizon) {
        for (h in seq_len(horizon)) {
          pg_h <- make_nemenyi_page(
            df_lev %>% dplyr::filter(horizon == h),
            m$col, m$lower_better,
            title_str = sprintf("Nemenyi Framework: %s | Level: %s | Horizon h = %d", m$label, lev, h)
          )
          print(pg_h$grob)
        }
      }
    }
  }
  
  # ── Final page: relative line plot ───────────────────────────────────────
  rel_plot <- make_rel_lineplot(summary_lm, m$col, m$label,
                                m$lower_better, m$coverage_target)
  print(plot_grid(rel_plot, shared_legend,
                  ncol = 1, rel_heights = c(0.80, 0.20)))
  
  dev.off()
  message("  Saved: ", pdf_path)
}

message("All PDFs written to ", dir_out, "/")