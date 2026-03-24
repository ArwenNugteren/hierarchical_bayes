library(fable)
library(dplyr)
library(stringr)
library(tidyr)
library(tsibble)

### load in data
tourism_bottom <- read.csv("tourism_data.csv")

# Get all bottom-level series (excluding Year and Month)
series_cols <- setdiff(names(tourism_bottom), c("Year", "Month"))

# Helper to sum columns matching a pattern
sum_cols <- function(df, pattern) {
  cols <- grep(pattern, names(df), value = TRUE)
  rowSums(df[, cols, drop = FALSE])
}

tourism <- tourism_bottom %>%
  mutate(
    Total = rowSums(across(all_of(series_cols))),
    TotalHol = sum_cols(., "Hol$"),
    TotalVis = sum_cols(., "Vis$"),
    TotalBus = sum_cols(., "Bus$"),
    TotalOth = sum_cols(., "Oth$"),
  ) %>%
  { 
    df <- .
    states <- unique(str_extract(series_cols, "^[A-Z]"))
    for (state in states) {
      for (purpose in c("Hol", "Vis", "Bus", "Oth")) {
        col_name <- paste0(state, purpose)
        pattern  <- paste0("^", state, "..", purpose, "$")
        df[[col_name]] <- sum_cols(tourism_bottom, pattern)
      }
    }
    df
  } %>%
  {
    df <- .
    zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
    for (zone in zones) {
      for (purpose in c("Hol", "Vis", "Bus", "Oth")) {
        col_name <- paste0(zone, purpose)
        pattern  <- paste0("^", zone, ".", purpose, "$")
        df[[col_name]] <- sum_cols(tourism_bottom, pattern)
      }
    }
    df
  }

states <- unique(str_extract(series_cols, "^[A-Z]"))
for (state in states) {
  col_name <- state
  pattern  <- paste0("^", state, "..", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
for (zone in zones) {
  col_name <- zone
  pattern  <- paste0("^", zone, ".", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

regions <- unique(str_extract(series_cols, "^[A-Z]{3}"))
for (region in regions) {
  col_name <- region
  pattern  <- paste0("^", region, "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

purposes <- c("Hol", "Vis", "Bus", "Oth")
lvl0 <- "Total"
lvl1 <- paste0("Total", purposes)
lvl2 <- as.vector(sapply(states,  function(s) c(s, paste0(s, purposes))))
lvl3 <- as.vector(sapply(zones,   function(z) c(z, paste0(z, purposes))))
lvl4 <- as.vector(sapply(regions, function(r) c(r, paste0(r, purposes))))
ordered_cols <- c("Year", "Month", lvl0, lvl1, lvl2, lvl3, lvl4)
tourism <- tourism[, ordered_cols]

tourism <- tourism %>%
  fill(Year, .direction = "down") %>%
  mutate(Date = yearmonth(paste(Year, Month))) %>%
  select(-Year, -Month) %>%
  as_tsibble(index = Date)

# ── Select series of interest ──────────────────────────────────────────────────

state_codes <- sort(unique(str_extract(series_cols, "^[A-Z]")))  # A through G
series_of_interest <- c("Total", state_codes)

training_end <- 96
horizon      <- 12

# ── Pivot to long format for fable ────────────────────────────────────────────

tourism_long <- tourism %>%
  as_tibble() %>%
  select(Date, all_of(series_of_interest)) %>%
  pivot_longer(-Date, names_to = "series", values_to = "value") %>%
  as_tsibble(index = Date, key = series)

# ── Fit ETS models (auto-selected) on training data ───────────────────────────

fit <- tourism_long %>%
  filter_index(. ~ as.character(tourism$Date[training_end])) %>%
  model(ets = ETS(value))

report(fit)

# ── Information criteria ───────────────────────────────────────────────────────

ic_df <- glance(fit) %>%
  select(series, AIC, AICc, BIC)

# Sum across all 8 series for comparison with Bayesian model
ic_summary <- ic_df %>%
  summarise(
    AIC  = sum(AIC),
    AICc = sum(AICc),
    BIC  = sum(BIC)
  )

cat("\nPer-series information criteria:\n")
print(ic_df)
cat("\nSummed across all series (comparable to Bayesian):\n")
print(ic_summary)

# ── Forecasts ─────────────────────────────────────────────────────────────────

fc <- fit %>%
  forecast(h = horizon)

# ── Actuals ───────────────────────────────────────────────────────────────────

actuals_long <- tourism_long %>%
  filter_index(as.character(tourism$Date[training_end + 1]) ~
                 as.character(tourism$Date[training_end + horizon]))

# ── Accuracy ──────────────────────────────────────────────────────────────────

accuracy_df <- fc %>%
  as_tibble() %>%
  select(series, Date, .mean) %>%
  rename(forecast = .mean) %>%
  left_join(actuals_long %>% as_tibble() %>% select(series, Date, value),
            by = c("series", "Date")) %>%
  rename(actual = value) %>%
  group_by(series) %>%
  mutate(horizon = row_number(),
         error     = forecast - actual,
         abs_error = abs(error),
         sq_error  = error^2) %>%
  ungroup()

summary_accuracy <- accuracy_df %>%
  group_by(series) %>%
  summarise(MAE  = mean(abs_error),
            RMSE = sqrt(mean(sq_error)),
            .groups = "drop") %>%
  bind_rows(
    accuracy_df %>%
      summarise(series = "Overall",
                MAE    = mean(abs_error),
                RMSE   = sqrt(mean(sq_error)))
  )

cat("\nPer-series MAE and RMSE:\n")
print(summary_accuracy)