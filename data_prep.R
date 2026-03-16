library(fable)
library(dplyr)
library(stringr)

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

# ── Reorder columns ───────────────────────────────────────────────────────────
purposes <- c("Hol", "Vis", "Bus", "Oth")

# Level 0 & 1
lvl0 <- "Total"
lvl1 <- paste0("Total", purposes)

# Level 2: single-letter state codes × purpose
states <- unique(str_extract(series_cols, "^[A-Z]"))
lvl2   <- as.vector(outer(states, purposes, paste0))

# Level 3: two-letter zone codes × purpose
zones <- unique(str_extract(series_cols, "^[A-Z]{2}"))
lvl3  <- as.vector(outer(zones, purposes, paste0))

# Level 4: three-letter region codes × purpose (original bottom-level columns)
regions <- unique(str_extract(series_cols, "^[A-Z]{3}"))
lvl4    <- as.vector(outer(regions, purposes, paste0))

ordered_cols <- c("Year", "Month", lvl0, lvl1, lvl2, lvl3, lvl4)

tourism <- tourism[, ordered_cols]

# replace NA years by back-filling and make a single time column
tourism <- tourism %>%
  fill(Year, .direction="down") %>%
  mutate(Month = paste0(Year, "_", Month)) %>%
  rename(Date = Month) %>%
  select(-Year)

save(tourism, file="tourism.Rdata")
