# read in original data
synthetic_data <- read.csv("data/synthetic_hierarchy.csv")

# add noise to all levels/series so that the noise is roughly within +/- 5%
# of each series' range. We treat +/- 3 SD as corresponding to that 5%
# window, i.e. sd = (0.05 * range) / 3, computed separately per column.
# A fresh noise draw is generated for every row (time point) and every
# column, so the added noise varies over time as well as across series.

set.seed(2025)

noisy_data <- synthetic_data

for (col in names(synthetic_data)) {
  series <- synthetic_data[[col]]
  
  col_range <- max(series, na.rm = TRUE) - min(series, na.rm = TRUE)
  noise_sd <- (0.05 * col_range) / 3
  
  noise <- rnorm(n = length(series), mean = 0, sd = noise_sd)
  
  noisy_data[[col]] <- series + noise
}

# save the resulting (now incoherent, since noise is added independently
# per series/row rather than respecting the aggregation structure) dataset
write.csv(noisy_data, "data/incoherent_synthetic_hierarchy.csv", row.names = FALSE)