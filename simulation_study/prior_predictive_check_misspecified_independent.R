# =============================================================================
# Prior predictive check — sampling + plotting
# =============================================================================
# Reuses prior_predictive_check.stan, which has no likelihood: sampling from
# it draws directly from the prior, and generated quantities forward-simulates
# a full length-T path per replication for every series in the hierarchy.

library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
# 1. LOAD DATA (same convention as the main fitting script)
# =============================================================================

synthetic_data <- read.csv("synthetic_hierarchy.csv")
data_mat       <- as.matrix(synthetic_data)   # columns: [Total, A, B, AA, AB, BA, BB]

IDX <- list(total = 1L, A = 2L, B = 3L, bot = 4:7)
bot_names <- c("AA", "AB", "BA", "BB")

# Length of the check: match this to whatever window you want to validate,
# e.g. TRAIN_SIZE from the main script, or the full observed history.
T_CHECK <- min(50L, nrow(data_mat))

# =============================================================================
# 2. SAMPLE FROM THE PRIOR
# =============================================================================
# No likelihood means no real posterior geometry to fight — this is cheap.
# Multiple chains x moderate iterations already gives a few thousand draws.


model_string <- "
// =============================================================================
// Prior predictive check for the Hierarchical ETS(A,A) Model
// =============================================================================
// Omits observed time-series data and likelihood evaluations. 
// Parameters and hyperpriors are sampled in model{}, and 
// generated quantities{} forward-simulates the length-T in-sample paths 
// and length-H out-of-sample forecasts.

// =============================================================================
// Prior Predictive Check for Hierarchical ETS(A,A) Model
// =============================================================================
// Generates prior sample paths for in-sample length T.
// No observational data, likelihood, or out-of-sample forecasting included.

data {
  int<lower=1> T;  // length of the simulated check series
}

parameters {
  // ── Independent Smoothing Parameters ──────────────────────────────────────
  // Bottom-level parameters: [AA, AB, BA, BB]
  vector<lower=0, upper=1>[4] alpha_bot;
  vector<lower=0, upper=1>[4] beta_bot;

  // Mid-level parameters
  real<lower=0, upper=1> alpha_A;
  real<lower=0, upper=1> beta_A;

  real<lower=0, upper=1> alpha_B;
  real<lower=0, upper=1> beta_B;

  // Top-level parameters
  real<lower=0, upper=1> alpha_top;
  real<lower=0, upper=1> beta_top;

  // ── Laplace noise scales (independent per series) ─────────────────────────
  vector<lower=0>[4] sigma_bottom;
  real<lower=0> sigma_A;
  real<lower=0> sigma_B;
  real<lower=0> sigma_top;

  // ── Initial State parameters ──────────────────────────────────────────────
  vector[4] l0_bot;
  vector[4] b0_bot;
  real l0_A; real b0_A;
  real l0_B; real b0_B;
  real l0_top; real b0_top;
}

model {
  // ── Priors only — this block IS the target distribution ──────────────────
  alpha_bot ~ uniform(0, 1);
  beta_bot  ~ uniform(0, 1);

  alpha_A   ~ uniform(0, 1);
  beta_A    ~ uniform(0, 1);

  alpha_B   ~ uniform(0, 1);
  beta_B    ~ uniform(0, 1);

  alpha_top ~ uniform(0, 1);
  beta_top  ~ uniform(0, 1);

  l0_bot ~ normal(0, 10); b0_bot ~ normal(0, 5);
  l0_A   ~ normal(0, 10); b0_A   ~ normal(0, 5);
  l0_B   ~ normal(0, 10); b0_B   ~ normal(0, 5);
  l0_top ~ normal(0, 10); b0_top ~ normal(0, 5);

  sigma_bottom ~ cauchy(0, 1);
  sigma_A      ~ cauchy(0, 1);
  sigma_B      ~ cauchy(0, 1);
  sigma_top    ~ cauchy(0, 1);
}

generated quantities {
  matrix[T, 4] Y_bottom_sim;
  vector[T]    Y_A_sim;
  vector[T]    Y_B_sim;
  vector[T]    Y_total_sim;

  {
    matrix[T, 4] l_bot_sim;
    matrix[T, 4] b_bot_sim;
    vector[T] l_A_sim;   vector[T] b_A_sim;
    vector[T] l_B_sim;   vector[T] b_B_sim;
    vector[T] l_top_sim; vector[T] b_top_sim;

    // ── Bottom Level Simulation ──────────────────────────────────────────────
    for (j in 1:4) {
      real e_1 = double_exponential_rng(0, sigma_bottom[j]);
      Y_bottom_sim[1, j] = l0_bot[j] + b0_bot[j] + e_1;
      l_bot_sim[1, j]    = l0_bot[j] + b0_bot[j] + alpha_bot[j] * e_1;
      b_bot_sim[1, j]    = b0_bot[j] + alpha_bot[j] * beta_bot[j] * e_1;

      for (t in 2:T) {
        real e_t = double_exponential_rng(0, sigma_bottom[j]);
        Y_bottom_sim[t, j] = l_bot_sim[t-1, j] + b_bot_sim[t-1, j] + e_t;
        l_bot_sim[t, j]    = l_bot_sim[t-1, j] + b_bot_sim[t-1, j] + alpha_bot[j] * e_t;
        b_bot_sim[t, j]    = b_bot_sim[t-1, j] + alpha_bot[j] * beta_bot[j] * e_t;
      }
    }

    // ── Mid Level A Simulation ───────────────────────────────────────────────
    {
      real e_1 = double_exponential_rng(0, sigma_A);
      Y_A_sim[1] = l0_A + b0_A + e_1;
      l_A_sim[1] = l0_A + b0_A + alpha_A * e_1;
      b_A_sim[1] = b0_A + alpha_A * beta_A * e_1;

      for (t in 2:T) {
        real e_t = double_exponential_rng(0, sigma_A);
        Y_A_sim[t] = l_A_sim[t-1] + b_A_sim[t-1] + e_t;
        l_A_sim[t] = l_A_sim[t-1] + b_A_sim[t-1] + alpha_A * e_t;
        b_A_sim[t] = b_A_sim[t-1] + alpha_A * beta_A * e_t;
      }
    }

    // ── Mid Level B Simulation ───────────────────────────────────────────────
    {
      real e_1 = double_exponential_rng(0, sigma_B);
      Y_B_sim[1] = l0_B + b0_B + e_1;
      l_B_sim[1] = l0_B + b0_B + alpha_B * e_1;
      b_B_sim[1] = b0_B + alpha_B * beta_B * e_1;

      for (t in 2:T) {
        real e_t = double_exponential_rng(0, sigma_B);
        Y_B_sim[t] = l_B_sim[t-1] + b_B_sim[t-1] + e_t;
        l_B_sim[t] = l_B_sim[t-1] + b_B_sim[t-1] + alpha_B * e_t;
        b_B_sim[t] = b_B_sim[t-1] + alpha_B * beta_B * e_t;
      }
    }

    // ── Top Level Simulation ─────────────────────────────────────────────────
    {
      real e_1 = double_exponential_rng(0, sigma_top);
      Y_total_sim[1] = l0_top + b0_top + e_1;
      l_top_sim[1]   = l0_top + b0_top + alpha_top * e_1;
      b_top_sim[1]   = b0_top + alpha_top * beta_top * e_1;

      for (t in 2:T) {
        real e_t = double_exponential_rng(0, sigma_top);
        Y_total_sim[t] = l_top_sim[t-1] + b_top_sim[t-1] + e_t;
        l_top_sim[t]   = l_top_sim[t-1] + b_top_sim[t-1] + alpha_top * e_t;
        b_top_sim[t]   = b_top_sim[t-1] + alpha_top * beta_top * e_t;
      }
    }
  }
}
"


ppc_model <- stan_model(model_code = model_string)

ppc_fit <- sampling(
  object  = ppc_model,
  data    = list(T = T_CHECK),
  seed    = 2025,
  chains  = 5,
  iter    = 2500,
  warmup  = 500,
  cores   = 5,
  refresh = 0
)

# =============================================================================
# 3. EXTRACT REPLICATIONS
# =============================================================================

Y_bottom_sim <- rstan::extract(ppc_fit, "Y_bottom_sim")$Y_bottom_sim  # [draws, T, 4]
Y_A_sim      <- rstan::extract(ppc_fit, "Y_A_sim")$Y_A_sim            # [draws, T]
Y_B_sim      <- rstan::extract(ppc_fit, "Y_B_sim")$Y_B_sim            # [draws, T]
Y_total_sim  <- rstan::extract(ppc_fit, "Y_total_sim")$Y_total_sim    # [draws, T]

n_draws <- dim(Y_bottom_sim)[1]
message(sprintf("Extracted %d prior predictive replications of length %d", n_draws, T_CHECK))

# =============================================================================
# 4. BUILD QUANTILE RIBBONS PER TIMEPOINT
# =============================================================================

# Include quantiles for 99% (0.5% - 99.5%), 95% (2.5% - 97.5%), 90% (5% - 95%), and 50% (25% - 75%)
qs <- c(0.005, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.995)

ribbon_from_draws <- function(draws_mat, series_name) {
  # draws_mat: [n_draws x T]
  apply(draws_mat, 2, quantile, probs = qs) %>%
    t() %>%
    as.data.frame() %>%
    setNames(c("q005", "q025", "q05", "q25", "q50", "q75", "q95", "q975", "q995")) %>%
    mutate(t = row_number(), series = series_name)
}

ribbon_AA    <- ribbon_from_draws(Y_bottom_sim[, , 1], "AA (bottom, representative)")
ribbon_A     <- ribbon_from_draws(Y_A_sim,             "A (mid, representative)")
ribbon_total <- ribbon_from_draws(Y_total_sim,         "Total (top)")

# Define explicit top-to-bottom factor levels for facet ordering
series_levels <- c(
  "Total (top)",
  "A (mid, representative)",
  "AA (bottom, representative)"
)

ribbons <- bind_rows(ribbon_AA, ribbon_A, ribbon_total) %>%
  mutate(series = factor(series, levels = series_levels))

# Matching observed series, aligned to the same t = 1..T_CHECK index
actual <- bind_rows(
  data.frame(t = seq_len(T_CHECK), y = data_mat[seq_len(T_CHECK), IDX$bot[1]],
             series = "AA (bottom, representative)"),
  data.frame(t = seq_len(T_CHECK), y = data_mat[seq_len(T_CHECK), IDX$A],
             series = "A (mid, representative)"),
  data.frame(t = seq_len(T_CHECK), y = data_mat[seq_len(T_CHECK), IDX$total],
             series = "Total (top)")
) %>%
  mutate(series = factor(series, levels = series_levels))

# =============================================================================
# 5. FAN CHART: replications vs. observed data (with interval legend)
# =============================================================================

# Define manual fill scale palette mapping interval labels to colors/alphas
interval_colors <- c(
  "99% Interval" = alpha("steelblue", 0.15),
  "95% Interval" = alpha("steelblue", 0.25),
  "90% Interval" = alpha("steelblue", 0.40),
  "50% Interval" = alpha("steelblue", 0.60)
)

fan_plot <- ggplot(ribbons, aes(x = t)) +
  # Map fill aesthetic to strings so ggplot builds the discrete legend
  geom_ribbon(aes(ymin = q005, ymax = q995, fill = "99% Interval")) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = "95% Interval")) +
  geom_ribbon(aes(ymin = q05,  ymax = q95,  fill = "90% Interval")) +
  geom_ribbon(aes(ymin = q25,  ymax = q75,  fill = "50% Interval")) +
  # Line layers with legend key for the actual path
  geom_line(aes(y = q50, colour = "Prior Median"), linewidth = 0.6) +
  geom_line(data = actual, aes(x = t, y = y, colour = "Observed Data"), linewidth = 0.6) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  scale_fill_manual(
    name = "Predictive Intervals",
    values = interval_colors,
    breaks = c("99% Interval", "95% Interval", "90% Interval", "50% Interval")
  ) +
  scale_colour_manual(
    name = "Series",
    values = c("Prior Median" = "steelblue4", "Observed Data" = "black")
  ) +
  labs(
    title    = "Prior predictive check: Replications vs. observed data",
    subtitle = sprintf("%d prior predictive replications, length %d", n_draws, T_CHECK),
    x = "t", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10, face = "bold"),
    legend.margin = margin(t = -5)
  )

ggsave("prior_predictive_misspecified_independent.png", fan_plot, width = 10, height = 10, dpi = 150)


# =============================================================================
# 6. CROSS-LEVEL COHERENCE CHECK
# =============================================================================
# The model does not structurally enforce Total = A + B = sum(bottom) at the
# generative level (each level has independent priors). This checks whether
# the priors at least imply comparable *scale* across levels — a mismatch
# here means the reconciliation step (MinT) will be pulling badly-matched
# distributions into agreement, not a modest correction.

sum_bottom_draws <- apply(Y_bottom_sim, c(1, 2), sum)   # [draws x T], sum over the 4 bottom series
sum_mid_draws    <- Y_A_sim + Y_B_sim                   # [draws x T]

coherence_df <- bind_rows(
  data.frame(value = as.vector(sum_bottom_draws[, T_CHECK]), group = "sum(bottom) at final t"),
  data.frame(value = as.vector(sum_mid_draws[, T_CHECK]),    group = "A + B at final t"),
  data.frame(value = as.vector(Y_total_sim[, T_CHECK]),      group = "Total at final t")
)

coherence_plot <- ggplot(coherence_df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Cross-level scale consistency under the prior",
    subtitle = "If these are badly mismatched, priors imply an incoherent hierarchy before any reconciliation",
    x = NULL, y = "density", fill = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave("prior_predictive_coherence.png", coherence_plot, width = 8, height = 5, dpi = 150)

# Quick numeric summary alongside the plot
coherence_summary <- coherence_df %>%
  group_by(group) %>%
  summarise(
    median = median(value),
    iqr    = IQR(value),
    p05    = quantile(value, 0.05),
    p95    = quantile(value, 0.95),
    .groups = "drop"
  )

print(coherence_summary)

message("Saved: prior_predictive_fanchart.png, prior_predictive_coherence.png")