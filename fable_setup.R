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
    # ── Level 0: Grand total ──────────────────────────────────────────
    Total = rowSums(across(all_of(series_cols))),
    
    # ── Level 1: Purpose of travel (across all states) ───────────────
    TotalHol = sum_cols(., "Hol$"),
    TotalVis = sum_cols(., "Vis$"),
    TotalBus = sum_cols(., "Bus$"),
    TotalOth = sum_cols(., "Oth$"),
    
    # ── Level 2: State × Purpose ──────────────────────────────────────
    # Extract unique state letters (1st character of the 3-letter code)
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
    # ── Level 3: State × Zone × Purpose ──────────────────────────────
    # Extract unique state+zone pairs (first 2 characters)
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
# ── Level 2b: State totals (across all purposes) ─────────────────────────────
for (state in states) {
  col_name <- state
  pattern  <- paste0("^", state, "..", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
# ── Level 3b: Zone totals (across all purposes) ───────────────────────────────
for (zone in zones) {
  col_name <- zone
  pattern  <- paste0("^", zone, ".", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

# ── Level 4b: Region totals (across all purposes) ────────────────────────────
regions <- unique(str_extract(series_cols, "^[A-Z]{3}"))
for (region in regions) {
  col_name <- region
  pattern  <- paste0("^", region, "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

purposes <- c("Hol", "Vis", "Bus", "Oth")

lvl0 <- "Total"
lvl1 <- paste0("Total", purposes)

# For each level, put the aggregate total first, then the four purpose splits
lvl2 <- as.vector(sapply(states,  function(s) c(s, paste0(s, purposes))))
lvl3 <- as.vector(sapply(zones,   function(z) c(z, paste0(z, purposes))))
lvl4 <- as.vector(sapply(regions, function(r) c(r, paste0(r, purposes))))

ordered_cols <- c("Year", "Month", lvl0, lvl1, lvl2, lvl3, lvl4)

tourism <- tourism[, ordered_cols]

# replace NA years by back-filling and make a single time column
tourism <- tourism %>%
  fill(Year, .direction = "down") %>%
  # Convert the strings directly into a Year-Month objects
  mutate(Date = yearmonth(paste(Year, Month))) %>% 
  select(-Year, -Month) %>% 
  as_tsibble(index = Date)

#### qld tourism
# qld_tourism <- tourism %>%
#   select(C)
# 
# qld_tourism <- qld_tourism[1:96,]
# 
# qld_model <- qld_tourism %>%
#   model(ETS(C ~ error("A") + season("A")))
# report(qld_model)
# 
# # ── Forecast ──────────────────────────────────────────────────────────────────
# 
# fc <- qld_model %>%
#   forecast(h = 12)
# 
# # ── Actuals ───────────────────────────────────────────────────────────────────
# 
# actuals <- tourism %>%
#   select(C) %>%
#   slice(97:108)
# 
# # ── Accuracy ─────────────────────────────────────────────────────────────────
# 
# accuracy_df <- fc %>%
#   as_tibble() %>%
#   select(Date, .mean) %>%
#   rename(forecast = .mean) %>%
#   bind_cols(actual = actuals$C) %>%
#   mutate(
#     horizon   = 1:12,
#     error     = forecast - actual,
#     abs_error = abs(error),
#     sq_error  = error^2
#   )
# 
# mae  <- mean(accuracy_df$abs_error)
# rmse <- sqrt(mean(accuracy_df$sq_error))
# 
# cat(sprintf("MAE:  %.4f\n", mae))
# cat(sprintf("RMSE: %.4f\n", rmse))
# print(accuracy_df)



### aus total tourism

aus_tourism <- tourism %>%
  select(Total)

Y_mean <- mean(aus_tourism$Total)
Y_sd   <- sd(aus_tourism$Total)
# Y_mean <- 0
# Y_sd <- 1

aus_tourism <- aus_tourism[1:96,] %>%
  mutate(scaled_total = (Total - Y_mean)/Y_sd)

aus_model <- aus_tourism %>%
  model(ETS(scaled_total ~ error("A") + season("A")))
report(aus_model)
# aus_model <- aus_tourism %>%
#   model(ETS(scaled_total))
# report(aus_model)

# ── Forecast ──────────────────────────────────────────────────────────────────

fc <- aus_model %>%
  forecast(h = 12)

# ── Actuals ───────────────────────────────────────────────────────────────────

actuals <- tourism %>%
  select(Total) %>%
  slice(97:108)

# ── Accuracy ─────────────────────────────────────────────────────────────────

accuracy_df <- fc %>%
  as_tibble() %>%
  select(Date, .mean) %>%
  rename(forecast = .mean) %>%
  mutate(forecast = forecast*Y_sd + Y_mean) %>%
  bind_cols(actual = actuals$Total) %>%
  mutate(
    horizon   = 1:12,
    error     = forecast - actual,
    abs_error = abs(error),
    sq_error  = error^2
  )

mae  <- mean(accuracy_df$abs_error)
rmse <- sqrt(mean(accuracy_df$sq_error))

cat(sprintf("MAE:  %.4f\n", mae))
cat(sprintf("RMSE: %.4f\n", rmse))
print(accuracy_df)
