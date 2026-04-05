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
      for (purpose in c("Hol", "Vis", "Bus", "Oth"))
        df[[paste0(state, purpose)]] <- sum_cols(tourism_bottom, paste0("^", state, "..", purpose, "$"))
    df
  } %>%
  {
    df <- .
    zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
    for (zone in zones)
      for (purpose in c("Hol", "Vis", "Bus", "Oth"))
        df[[paste0(zone, purpose)]] <- sum_cols(tourism_bottom, paste0("^", zone, ".", purpose, "$"))
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

state_codes <- sort(unique(str_extract(series_cols, "^[A-Z]")))
zone_codes  <- sort(unique(str_extract(series_cols, "^[A-Z]{2}")))

window_size <- 96
horizon     <- 12
test_set    <- 121
origins     <- window_size:(window_size + test_set - 1)

# ── Series label helper ────────────────────────────────────────────────────────

series_label <- function(state_col, zone_col) {
  case_when(
    is_aggregated(state_col) & is_aggregated(zone_col) ~ "Total",
    is_aggregated(zone_col)                            ~ as.character(state_col),
    TRUE                                               ~ as.character(zone_col)
  )
}

# ── Rolling window loop ────────────────────────────────────────────────────────

all_base_results <- vector("list", length(origins))
all_mint_results <- vector("list", length(origins))

for (idx in seq_along(origins)) {
  
  origin      <- origins[idx]
  train_dates <- tourism$Date[(origin - window_size + 1):origin]
  test_dates  <- tourism$Date[(origin + 1):(origin + horizon)]
  
  # ── Build two-level hierarchical tsibble ────────────────────────────────────
  # zone is nested within state
  
  tourism_hts <- tourism %>%
    as_tibble() %>%
    select(Date, all_of(zone_codes)) %>%
    pivot_longer(-Date, names_to = "zone", values_to = "value") %>%
    mutate(state = substr(zone, 1, 1)) %>%
    as_tsibble(index = Date, key = c(state, zone)) %>%
    aggregate_key(state / zone, value = sum(value))
  
  train_hts <- tourism_hts %>% filter(Date %in% train_dates)
  
  # ── Actuals ────────────────────────────────────────────────────────────────
  
  actuals_df <- tourism_hts %>%
    filter(Date %in% test_dates) %>%
    as_tibble() %>%
    mutate(series = series_label(state, zone)) %>%
    select(series, Date, actual = value)
  
  # ── Fit ────────────────────────────────────────────────────────────────────
  
  fit_hts <- train_hts %>%
    model(ets = ETS(value))
  
  # ── Base forecasts ──────────────────────────────────────────────────────────
  
  fc_base_raw <- fit_hts %>% forecast(h = horizon)
  
  fc_base <- fc_base_raw %>%
    as_tibble() %>%
    mutate(
      series   = series_label(state, zone),
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
  
  # ── MinT reconciled forecasts ───────────────────────────────────────────────
  
  fc_mint_raw <- fit_hts %>%
    reconcile(mint = min_trace(ets, method = "mint_shrink")) %>%
    forecast(h = horizon) %>%
    filter(.model == "mint")
  
  fc_mint <- fc_mint_raw %>%
    as_tibble() %>%
    mutate(
      series   = series_label(state, zone),
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

# ── Add level classification ───────────────────────────────────────────────────

add_level <- function(df) {
  df %>%
    mutate(level = case_when(
      series == "Total"   ~ "Total",
      nchar(series) == 1  ~ "State",
      nchar(series) == 2  ~ "Zone",
      TRUE                ~ "Other"
    ))
}

all_base_df <- add_level(all_base_df)
all_mint_df <- add_level(all_mint_df)

# ── Summary functions ──────────────────────────────────────────────────────────

summarise_by_series <- function(df) {
  df %>%
    group_by(level, series) %>%
    summarise(MAE      = mean(abs_error),
              RMSE     = sqrt(mean(sq_error)),
              CRPS     = mean(crps),
              Coverage = mean(covered),
              .groups  = "drop") %>%
    bind_rows(
      df %>%
        group_by(level) %>%
        summarise(series   = "Overall (level)",
                  MAE      = mean(abs_error),
                  RMSE     = sqrt(mean(sq_error)),
                  CRPS     = mean(crps),
                  Coverage = mean(covered),
                  .groups  = "drop")
    ) %>%
    bind_rows(
      df %>%
        summarise(level    = "All",
                  series   = "Overall",
                  MAE      = mean(abs_error),
                  RMSE     = sqrt(mean(sq_error)),
                  CRPS     = mean(crps),
                  Coverage = mean(covered))
    )
}

summarise_by_horizon <- function(df) {
  bind_rows(
    df %>%
      group_by(level, horizon) %>%
      summarise(MAE      = mean(abs_error),
                RMSE     = sqrt(mean(sq_error)),
                CRPS     = mean(crps),
                Coverage = mean(covered),
                .groups  = "drop"),
    df %>%
      group_by(horizon) %>%
      summarise(level    = "All",
                MAE      = mean(abs_error),
                RMSE     = sqrt(mean(sq_error)),
                CRPS     = mean(crps),
                Coverage = mean(covered),
                .groups  = "drop")
  )
}

summarise_by_series_and_horizon <- function(df) {
  df %>%
    group_by(level, series, horizon) %>%
    summarise(MAE      = mean(abs_error),
              RMSE     = sqrt(mean(sq_error)),
              CRPS     = mean(crps),
              Coverage = mean(covered),
              .groups  = "drop")
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
     file = "fable_zone_hierarchy_121_origins.RData")

print(base_series_accuracy)
print(mint_series_accuracy)
print(base_horizon_accuracy)
print(mint_horizon_accuracy)