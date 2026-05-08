library(fable)
library(dplyr)
library(stringr)
library(tidyr)
library(tsibble)
library(scoringRules)

### load in data
tourism_bottom <- read.csv("tourism_data.csv")

series_cols <- setdiff(names(tourism_bottom), c("Year", "Month"))

sum_cols <- function(df, pattern) {
  cols <- grep(pattern, names(df), value = TRUE)
  rowSums(df[, cols, drop = FALSE])
}

tourism <- tourism_bottom %>%
  mutate(
    Total    = rowSums(across(all_of(series_cols))),
    TotalHol = sum_cols(., "Hol$"),
    TotalVis = sum_cols(., "Vis$"),
    TotalBus = sum_cols(., "Bus$"),
    TotalOth = sum_cols(., "Oth$"),
  ) %>%
  {
    df <- .
    states <- unique(str_extract(series_cols, "^[A-Z]"))
    for (state in states)
      for (purpose in c("Hol", "Vis", "Bus", "Oth")) {
        df[[paste0(state, purpose)]] <- sum_cols(tourism_bottom, paste0("^", state, "..", purpose, "$"))
      }
    df
  } %>%
  {
    df <- .
    zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
    for (zone in zones)
      for (purpose in c("Hol", "Vis", "Bus", "Oth")) {
        df[[paste0(zone, purpose)]] <- sum_cols(tourism_bottom, paste0("^", zone, ".", purpose, "$"))
      }
    df
  }

states <- unique(str_extract(series_cols, "^[A-Z]"))
for (state in states)
  tourism[[state]] <- sum_cols(tourism_bottom, paste0("^", state, "..", "(Hol|Vis|Bus|Oth)$"))

zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
for (zone in zones)
  tourism[[zone]] <- sum_cols(tourism_bottom, paste0("^", zone, ".", "(Hol|Vis|Bus|Oth)$"))

regions <- unique(str_extract(series_cols, "^[A-Z]{3}"))
for (region in regions)
  tourism[[region]] <- sum_cols(tourism_bottom, paste0("^", region, "(Hol|Vis|Bus|Oth)$"))

purposes     <- c("Hol", "Vis", "Bus", "Oth")
lvl0         <- "Total"
lvl1         <- paste0("Total", purposes)
lvl2         <- as.vector(sapply(states,  function(s) c(s, paste0(s, purposes))))
lvl3         <- as.vector(sapply(zones,   function(z) c(z, paste0(z, purposes))))
lvl4         <- as.vector(sapply(regions, function(r) c(r, paste0(r, purposes))))
ordered_cols <- c("Year", "Month", lvl0, lvl1, lvl2, lvl3, lvl4)
tourism      <- tourism[, ordered_cols]

tourism <- tourism %>%
  fill(Year, .direction = "down") %>%
  mutate(Date = yearmonth(paste(Year, Month))) %>%
  select(-Year, -Month) %>%
  as_tsibble(index = Date)

# ── Setup ──────────────────────────────────────────────────────────────────────

state_codes  <- sort(unique(str_extract(series_cols, "^[A-Z]")))
window_size  <- 96
horizon      <- 12
test_set     <- 121
origins      <- window_size:(window_size + test_set - 1)

# ── Helper: compute accuracy from forecast and actuals ─────────────────────────
# fc_df must contain: series, Date, forecast, q05, q95, crps

compute_accuracy <- function(fc_df, actuals_df) {
  fc_df %>%
    left_join(actuals_df, by = c("series", "Date")) %>%
    group_by(series) %>%
    mutate(
      horizon   = row_number(),
      error     = forecast - actual,
      abs_error = abs(error),
      sq_error  = error^2,
      covered   = as.integer(actual >= q05 & actual <= q95)
    ) %>%
    ungroup()
}

# ── Helper: extract forecast metrics from a fable forecast object ──────────────

extract_fc_metrics <- function(fc_fable, series_label_fn) {
  # Extract parameters of the normal predictive distribution
  fc_fable %>%
    as_tibble() %>%
    mutate(
      series   = series_label_fn(state),
      forecast = .mean,
      # fable ETS predictive distribution is normal — extract mean and sd
      fc_sd    = sqrt(distributional::variance(value)),
      q05      = qnorm(0.05, mean = .mean, sd = fc_sd),
      q95      = qnorm(0.95, mean = .mean, sd = fc_sd),
      crps     = scoringRules::crps_norm(
        y    = NA_real_,   # placeholder, filled after join
        mean = .mean,
        sd   = fc_sd
      )
    ) %>%
    select(series, Date, forecast, fc_sd, q05, q95)
}

# ── Rolling window loop ────────────────────────────────────────────────────────

all_base_results <- vector("list", length(origins))
all_mint_results <- vector("list", length(origins))

series_label <- function(state_col) {
  ifelse(is_aggregated(state_col), "Total", as.character(state_col))
}

for (idx in seq_along(origins)) {
  
  origin      <- origins[idx]
  train_dates <- tourism$Date[(origin - window_size + 1):origin]
  test_dates  <- tourism$Date[(origin + 1):(origin + horizon)]
  
  tourism_hts <- tourism %>%
    as_tibble() %>%
    select(Date, all_of(state_codes)) %>%
    pivot_longer(-Date, names_to = "state", values_to = "value") %>%
    as_tsibble(index = Date, key = state) %>%
    aggregate_key(state, value = sum(value))
  
  train_hts <- tourism_hts %>% filter(Date %in% train_dates)
  
  actuals_df <- tourism_hts %>%
    filter(Date %in% test_dates) %>%
    as_tibble() %>%
    mutate(series = series_label(state)) %>%
    select(series, Date, actual = value)
  
  fit_hts <- train_hts %>%
    model(ets = ETS(value))
  
  # ── Base forecasts ───────────────────────────────────────────────────────────
  
  fc_base_raw <- fit_hts %>% forecast(h = horizon)
  
  fc_base <- fc_base_raw %>%
    as_tibble() %>%
    mutate(
      series   = series_label(state),
      forecast = .mean,
      fc_sd    = sqrt(distributional::variance(value)),
      q05      = qnorm(0.05, mean = .mean, sd = fc_sd),
      q95      = qnorm(0.95, mean = .mean, sd = fc_sd)
    ) %>%
    select(series, Date, forecast, fc_sd, q05, q95)
  
  base_acc <- fc_base %>%
    left_join(actuals_df, by = c("series", "Date")) %>%
    group_by(series) %>%
    mutate(
      horizon   = row_number(),
      error     = forecast - actual,
      abs_error = abs(error),
      sq_error  = error^2,
      crps      = scoringRules::crps_norm(actual, mean = forecast, sd = fc_sd),
      covered   = as.integer(actual >= q05 & actual <= q95)
    ) %>%
    ungroup() %>%
    mutate(origin = origin)
  
  # ── MinT reconciled forecasts ────────────────────────────────────────────────
  
  fc_mint_raw <- fit_hts %>%
    reconcile(mint = min_trace(ets, method = "mint_shrink")) %>%
    forecast(h = horizon) %>%
    filter(.model == "mint")
  
  fc_mint <- fc_mint_raw %>%
    as_tibble() %>%
    mutate(
      series   = series_label(state),
      forecast = .mean,
      fc_sd    = sqrt(distributional::variance(value)),
      q05      = qnorm(0.05, mean = .mean, sd = fc_sd),
      q95      = qnorm(0.95, mean = .mean, sd = fc_sd)
    ) %>%
    select(series, Date, forecast, fc_sd, q05, q95)
  
  mint_acc <- fc_mint %>%
    left_join(actuals_df, by = c("series", "Date")) %>%
    group_by(series) %>%
    mutate(
      horizon   = row_number(),
      error     = forecast - actual,
      abs_error = abs(error),
      sq_error  = error^2,
      crps      = scoringRules::crps_norm(actual, mean = forecast, sd = fc_sd),
      covered   = as.integer(actual >= q05 & actual <= q95)
    ) %>%
    ungroup() %>%
    mutate(origin = origin)
  
  all_base_results[[idx]] <- base_acc
  all_mint_results[[idx]] <- mint_acc
}

all_base_df <- bind_rows(all_base_results)
all_mint_df <- bind_rows(all_mint_results)

# ── Coherence helper ───────────────────────────────────────────────────────────

compute_coherence <- function(df) {
  df %>%
    group_by(origin, horizon) %>%
    summarise(
      coherence_error = forecast[series == "Total"] - sum(forecast[series != "Total"]),
      .groups = "drop"
    )
}

summarise_coherence_by_horizon <- function(coh_df) {
  coh_df %>%
    group_by(horizon) %>%
    summarise(
      ME_coherence   = mean(coherence_error),
      MAE_coherence  = mean(abs(coherence_error)),
      RMSE_coherence = sqrt(mean(coherence_error^2)),
      .groups = "drop"
    )
}

summarise_coherence_overall <- function(coh_df) {
  coh_df %>%
    summarise(
      ME_coherence   = mean(coherence_error),
      MAE_coherence  = mean(abs(coherence_error)),
      RMSE_coherence = sqrt(mean(coherence_error^2))
    )
}

# ── Updated summary functions ──────────────────────────────────────────────────

summarise_by_series <- function(df) {
  coh_df  <- compute_coherence(df)
  coh_ovr <- summarise_coherence_overall(coh_df)
  
  df %>%
    group_by(series) %>%
    summarise(MAE      = mean(abs_error),
              RMSE     = sqrt(mean(sq_error)),
              CRPS     = mean(crps),
              Coverage = mean(covered),
              .groups  = "drop") %>%
    bind_rows(
      df %>%
        summarise(series   = "Overall",
                  MAE      = mean(abs_error),
                  RMSE     = sqrt(mean(sq_error)),
                  CRPS     = mean(crps),
                  Coverage = mean(covered))
    ) %>%
    mutate(
      ME_coherence   = if_else(series == "Overall", coh_ovr$ME_coherence,   NA_real_),
      MAE_coherence  = if_else(series == "Overall", coh_ovr$MAE_coherence,  NA_real_),
      RMSE_coherence = if_else(series == "Overall", coh_ovr$RMSE_coherence, NA_real_)
    )
}

summarise_by_horizon <- function(df) {
  coh_by_h <- compute_coherence(df) %>% summarise_coherence_by_horizon()
  
  df %>%
    group_by(horizon) %>%
    summarise(MAE      = mean(abs_error),
              RMSE     = sqrt(mean(sq_error)),
              CRPS     = mean(crps),
              Coverage = mean(covered),
              .groups  = "drop") %>%
    left_join(coh_by_h, by = "horizon")
}

summarise_by_series_and_horizon <- function(df) {
  coh_by_h <- compute_coherence(df) %>% summarise_coherence_by_horizon()
  
  df %>%
    group_by(series, horizon) %>%
    summarise(MAE      = mean(abs_error),
              RMSE     = sqrt(mean(sq_error)),
              CRPS     = mean(crps),
              Coverage = mean(covered),
              .groups  = "drop") %>%
    left_join(coh_by_h, by = "horizon")
  # Same note as before: coherence repeats across series within a horizon —
  # it reflects aggregate mismatch, not a per-series contribution
}

base_series_accuracy         <- summarise_by_series(all_base_df)
base_horizon_accuracy        <- summarise_by_horizon(all_base_df)
base_series_horizon_accuracy <- summarise_by_series_and_horizon(all_base_df)

mint_series_accuracy         <- summarise_by_series(all_mint_df)
mint_horizon_accuracy        <- summarise_by_horizon(all_mint_df)
mint_series_horizon_accuracy <- summarise_by_series_and_horizon(all_mint_df)

save(base_series_accuracy, base_horizon_accuracy, base_series_horizon_accuracy,
     mint_series_accuracy, mint_horizon_accuracy, mint_series_horizon_accuracy,
     all_base_df, all_mint_df,
     file = "fable_state_hierarchy_121_origins.RData")

print(base_series_accuracy)
print(mint_series_accuracy)
print(base_horizon_accuracy)
print(mint_horizon_accuracy)
