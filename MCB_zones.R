# ══════════════════════════════════════════════════════════════════════════════
# FIXED MCB Test By Horizon: Correctly unpacking nested RDS structures
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(stringr)
library(tsutils)

# ── 0. Configuration ───────────────────────────────────────────────────────────
bayesian_runs <- list(
  "Bayes default"        = "accuracy_zones",
  "Bayes total_lambda=2" = "accuracy_zones_top_lambda2",
  "Bayes state_lambda=2" = "accuracy_zones_state_lambda2"
)

fable_rdata_path <- "fable_zone_hierarchy_121_origins.RData"
output_dir       <- "mcb_results_zones_by_horizon"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

horizon_max <- 12L

# ── 1. Load Bayesian Results (FIXED UNPACKING) ─────────────────────────────────
load_bayesian_data <- function(folder, method_name) {
  files <- list.files(folder, pattern = "^accuracy_timestep_[0-9]+\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop(sprintf("No files found in '%s'", folder))
  
  uni_list   <- list()
  multi_list <- list()
  
  for (i in seq_along(files)) {
    obj <- readRDS(files[i])
    
    # Explicitly pull from the nested univariate data frame
    if ("univariate" %in% names(obj)) {
      uni_list[[i]] <- obj$univariate %>%
        transmute(
          origin, horizon, series, actual, fc_mean, ae, se,
          level = case_when(
            level == "Total" ~ "Total",
            level == "State" ~ "State",
            level == "Zone"  ~ "Zone",
            TRUE             ~ level
          ),
          method = method_name
        )
    }
    
    # Explicitly pull from the nested multivariate data frame
    if ("multivariate" %in% names(obj)) {
      multi_list[[i]] <- obj$multivariate %>%
        transmute(origin, horizon, es, level = "Overall", method = method_name)
    }
  }
  
  list(uni = bind_rows(uni_list), multi = bind_rows(multi_list))
}

bayes_raw <- mapply(load_bayesian_data, 
                    folder = bayesian_runs, 
                    method_name = names(bayesian_runs), 
                    SIMPLIFY = FALSE)

bayes_uni   <- bind_rows(lapply(bayes_raw, `[[`, "uni"))
bayes_multi <- bind_rows(lapply(bayes_raw, `[[`, "multi"))

# ── 2. Load Fable Results ──────────────────────────────────────────────────────
load(fable_rdata_path) 

prep_fable_uni <- function(df, method_name) {
  df %>%
    transmute(
      origin, horizon, series, actual, fc_mean = forecast,
      ae = abs_error, se = sq_error,
      level = case_when(
        series == "Total"             ~ "Total",
        str_detect(series, "^[A-G]$")   ~ "State",
        str_detect(series, "^[A-G][A-Z]$") ~ "Zone",
        TRUE                          ~ "Other"
      ),
      method = method_name
    )
}

fable_uni <- bind_rows(
  prep_fable_uni(all_base_df, "ETS (base)"),
  prep_fable_uni(all_mint_df, "ETS (MinT)")
)

# ── 3. Core Plotting Engine ────────────────────────────────────────────────────
plot_mcb_panel <- function(wide_df, title) {
  mat <- wide_df %>% select(-uid) %>% as.matrix()
  mat <- mat[complete.cases(mat), , drop = FALSE]
  
  if (nrow(mat) >= 3 && ncol(mat) >= 2) {
    tsutils::nemenyi(mat, conf.level = 0.95, plottype = "vmcb", main = title)
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient Data\nor Single Method", cex = 0.8, col = "grey50")
    title(main = title)
  }
}

# ── 4. Run Pipelines for Univariate Metrics (MAE & RMSE) ──────────────────────
uni_metrics <- c("MAE", "RMSE")

for (metric in uni_metrics) {
  pdf_path <- file.path(output_dir, sprintf("MCB_%s.pdf", metric))
  pdf(pdf_path, width = 16, height = 4.5)
  
  for (h in 1:horizon_max) {
    par(mfrow = c(1, 4), mar = c(3, 4, 4, 1))
    
    # Combine true un-nested univariate frames
    h_uni <- filter(bind_rows(bayes_uni, fable_uni), horizon == h)
    
    # 1. Overall
    wide_overall <- h_uni %>% group_by(origin, method) %>%
      summarise(val = if(metric == "MAE") mean(ae) else sqrt(mean(se)), .groups = "drop") %>%
      mutate(uid = as.character(origin)) %>% pivot_wider(names_from = method, values_from = val)
    
    # 2. Total
    wide_total <- h_uni %>% filter(level == "Total") %>% group_by(origin, method) %>%
      summarise(val = if(metric == "MAE") mean(ae) else sqrt(mean(se)), .groups = "drop") %>%
      mutate(uid = as.character(origin)) %>% pivot_wider(names_from = method, values_from = val)
    
    # 3. State
    wide_state <- h_uni %>% filter(level == "State") %>% group_by(origin, method) %>%
      summarise(val = if(metric == "MAE") mean(ae) else sqrt(mean(se)), .groups = "drop") %>%
      mutate(uid = as.character(origin)) %>% pivot_wider(names_from = method, values_from = val)
    
    # 4. Zone
    wide_zone <- h_uni %>% filter(level == "Zone") %>% group_by(origin, method) %>%
      summarise(val = if(metric == "MAE") mean(ae) else sqrt(mean(se)), .groups = "drop") %>%
      mutate(uid = as.character(origin)) %>% pivot_wider(names_from = method, values_from = val)
    
    plot_mcb_panel(wide_overall, sprintf("%s - H %d\nOverall System", metric, h))
    plot_mcb_panel(wide_total,   sprintf("%s - H %d\nTotal Level", metric, h))
    plot_mcb_panel(wide_state,   sprintf("%s - H %d\nStates", metric, h))
    plot_mcb_panel(wide_zone,    sprintf("%s - H %d\nZones", metric, h))
  }
  dev.off()
}

# ── 5. Run Pipeline for Energy Score (ES) ──────────────────────────────────────
pdf_path_es <- file.path(output_dir, "MCB_ES.pdf")
pdf(pdf_path_es, width = 6, height = 5)

for (h in 1:horizon_max) {
  par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
  h_multi <- filter(bayes_multi, horizon == h)
  
  if (nrow(h_multi) > 0) {
    wide_es <- h_multi %>%
      group_by(origin, method) %>%
      summarise(val = mean(es), .groups = "drop") %>%
      mutate(uid = as.character(origin)) %>%
      pivot_wider(names_from = method, values_from = val)
    
    # Note: Fable won't show up here unless you generate sample paths from Fable
    plot_mcb_panel(wide_es, sprintf("Energy Score - Horizon %d\n(Bayesian Joint Comparison)", h))
  } else {
    plot.new()
    title(main = sprintf("ES - Horizon %d [No Data]", h))
  }
}
dev.off()
