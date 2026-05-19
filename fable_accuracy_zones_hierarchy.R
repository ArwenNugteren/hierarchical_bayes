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

# ── Helper for Multivariate Energy Score via Simulation ────────────────────────
compute_energy_score_row <- function(means, sds, actuals, n_sim = 500) {
  n_dims <- length(means)
  # Draw each dimension separately, then bind into n_sim × n_dims matrix
  sim_samples <- vapply(
    seq_len(n_dims),
    function(i) rnorm(n_sim, mean = means[i], sd = sds[i]),
    numeric(n_sim)
  )
  # sim_samples is now n_sim × n_dims — correct layout for es_sample
  scoringRules::es_sample(y = actuals, dat = t(sim_samples))
}
# ── Rolling window loop ────────────────────────────────────────────────────────

all_base_results  <- vector("list", length(origins))
all_mint_results  <- vector("list", length(origins))
all_base_es_list  <- vector("list", length(origins))
all_mint_es_list  <- vector("list", length(origins))

for (idx in seq_along(origins)) {
  
  origin      <- origins[idx]
  train_dates <- tourism$Date[(origin - window_size + 1):origin]
  test_dates  <- tourism$Date[(origin + 1):(origin + horizon)]
  
  # ── Build two-level hierarchical tsibble ────────────────────────────────────
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
  fit_hts <- train_hts %>% model(ets = ETS(value))
  
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
  
  # ── Joint Multivariate Energy Score Pipeline (ADDED) ────────────────────────
  # Calculate joint metric step by step across each predictive step horizon (1 to 12)
  # ── Joint Multivariate Energy Score Pipeline (NO-JOIN VERSION) ──────────────
  base_es_v <- numeric(horizon)
  mint_es_v <- numeric(horizon)
  
  for (h_step in 1:horizon) {
    target_date <- test_dates[h_step]
    
    # 1. Filter the raw fable forecast objects for this specific month step
    fc_base_slice <- fc_base_raw %>% filter(Date == target_date)
    fc_mint_slice <- fc_mint_raw %>% filter(Date == target_date)
    
    # 2. Extract the forecast distributions directly from the fable objects
    # This preserves fable's internal ordering perfectly
    base_means <- fc_base_slice$.mean
    base_sds   <- sqrt(distributional::variance(fc_base_slice$value))
    
    mint_means <- fc_mint_slice$.mean
    mint_sds   <- sqrt(distributional::variance(fc_mint_slice$value))
    
    # 3. Pull actuals directly from tourism_hts matching the exact order of keys
    # by using the key structures present in the forecast slice
    actuals_slice <- tourism_hts %>% 
      filter(Date == target_date) %>%
      right_join(as_tibble(fc_base_slice) %>% select(state, zone), by = c("state", "zone"))
    
    actuals_v <- actuals_slice$value
    
    # 4. Final safety guard: If anything dropped out, verify lengths match
    if (length(actuals_v) == length(base_means)) {
      base_es_v[h_step] <- compute_energy_score_row(
        means   = base_means,
        sds     = base_sds,
        actuals = actuals_v
      )
    } else {
      base_es_v[h_step] <- NA
    }
    
    if (length(actuals_v) == length(mint_means)) {
      mint_es_v[h_step] <- compute_energy_score_row(
        means   = mint_means,
        sds     = mint_sds,
        actuals = actuals_v
      )
    } else {
      mint_es_v[h_step] <- NA
    }
  }
  
  all_base_es_list[[idx]] <- tibble(origin = origin, horizon = 1:horizon, es = base_es_v)
  all_mint_es_list[[idx]] <- tibble(origin = origin, horizon = 1:horizon, es = mint_es_v)
  
  all_base_results[[idx]] <- base_acc
  all_mint_results[[idx]] <- mint_acc
}

all_base_df <- bind_rows(all_base_results)
all_mint_df <- bind_rows(all_mint_results)

fable_base_es <- bind_rows(all_base_es_list) %>% mutate(level = "Overall")
fable_mint_es <- bind_rows(all_mint_es_list) %>% mutate(level = "Overall")

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

# Save everything including the explicit multivariate Energy Score frames
save(base_series_accuracy, base_horizon_accuracy, base_series_horizon_accuracy,
     mint_series_accuracy, mint_horizon_accuracy, mint_series_horizon_accuracy,
     all_base_df, all_mint_df, fable_base_es, fable_mint_es,
     file = "fable_zone_hierarchy_121_origins.RData")

print(base_series_accuracy)
print(mint_series_accuracy)