library(fable)
library(dplyr)
library(stringr)
library(tidyr)

### load in data
tourism_bottom <- read.csv("data/tourism_data.csv")

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

states  <- unique(str_extract(series_cols, "^[A-Z]"))
zones   <- unique(str_extract(series_cols, "^[A-Z]{2}"))
regions <- unique(str_extract(series_cols, "^[A-Z]{3}"))

# ── Level 2b: State totals (across all purposes) ─────────────────────────────
for (state in states) {
  col_name <- state
  pattern  <- paste0("^", state, "..", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

# ── Level 3b: Zone totals (across all purposes) ───────────────────────────────
for (zone in zones) {
  col_name <- zone
  pattern  <- paste0("^", zone, ".", "(Hol|Vis|Bus|Oth)$")
  tourism[[col_name]] <- sum_cols(tourism_bottom, pattern)
}

# ── Level 4b: Region totals (across all purposes) ────────────────────────────
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

# replace NA years (as year is only listed with the january column) 
# and make a single time column
tourism <- tourism %>%
  fill(Year, .direction="down") %>%
  mutate(Month = paste0(Year, "_", Month)) %>%
  rename(Date = Month) %>%
  select(-Year)

save(tourism, file="data/tourism.Rdata")
