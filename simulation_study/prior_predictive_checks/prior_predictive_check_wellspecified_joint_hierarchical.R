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
// Prior predictive check for the Hierarchical Bayesian ARMA(2,2)
// =============================================================================
// This model retains the full level-0 and level-1 hyperprior structure while
// omitting the observed time-series data and likelihood evaluations.
// Sampling from model{} draws directly from the joint prior distribution.
// generated quantities{} then forward-simulates a length-T path for every series.

data {
  int<lower=3> T;  // length of the simulated check series (e.g., TRAIN_SIZE)
}

parameters {
  // ── Grand-mean hyperpriors (level 0) ──────────────────────────────────────
  real<lower=-1, upper=1> mu_phi2;
  real<lower=-1, upper=1> mu_phi1_raw;
  real<lower=-1, upper=1> mu_theta2;
  real<lower=-1, upper=1> mu_theta1_raw;

  // ── Group-mean hyperpriors (level 1) ──────────────────────────────────────
  real<lower=-1, upper=1> mu_phi2_A;
  real<lower=-1, upper=1> mu_phi1_raw_A;
  real<lower=-1, upper=1> mu_theta2_A;
  real<lower=-1, upper=1> mu_theta1_raw_A;

  real<lower=-1, upper=1> mu_phi2_B;
  real<lower=-1, upper=1> mu_phi1_raw_B;
  real<lower=-1, upper=1> mu_theta2_B;
  real<lower=-1, upper=1> mu_theta1_raw_B;

  // ── Bottom-level parameters: [AA, AB, BA, BB] ─────────────────────────────
  vector<lower=-1, upper=1>[4] phi2_bot;
  vector<lower=-1, upper=1>[4] phi1_raw_bot;
  vector<lower=-1, upper=1>[4] theta2_bot;
  vector<lower=-1, upper=1>[4] theta1_raw_bot;

  // ── Mid-level parameters ──────────────────────────────────────────────────
  real<lower=-1, upper=1> phi2_A;
  real<lower=-1, upper=1> phi1_raw_A;
  real<lower=-1, upper=1> theta2_A;
  real<lower=-1, upper=1> theta1_raw_A;

  real<lower=-1, upper=1> phi2_B;
  real<lower=-1, upper=1> phi1_raw_B;
  real<lower=-1, upper=1> theta2_B;
  real<lower=-1, upper=1> theta1_raw_B;

  // ── Top-level parameters ──────────────────────────────────────────────────
  real<lower=-1, upper=1> phi2_top;
  real<lower=-1, upper=1> phi1_raw_top;
  real<lower=-1, upper=1> theta2_top;
  real<lower=-1, upper=1> theta1_raw_top;

  // ── Bottom-level innovation covariance ────────────────────────────────────
  cholesky_factor_corr[4] L_Omega;
  vector<lower=0>[4]      tau;

  // ── Upper-level noise ─────────────────────────────────────────────────────
  real<lower=0> sigma_A;
  real<lower=0> sigma_B;
  real<lower=0> sigma_top;
}

transformed parameters {
  real mu_phi1   = mu_phi1_raw   * (1 - fabs(mu_phi2));
  real mu_theta1 = mu_theta1_raw * (1 - fabs(mu_theta2));

  real mu_phi1_A   = mu_phi1_raw_A   * (1 - fabs(mu_phi2_A));
  real mu_theta1_A = mu_theta1_raw_A * (1 - fabs(mu_theta2_A));
  real mu_phi1_B   = mu_phi1_raw_B   * (1 - fabs(mu_phi2_B));
  real mu_theta1_B = mu_theta1_raw_B * (1 - fabs(mu_theta2_B));

  vector[4] phi1_bot;
  vector[4] theta1_bot;
  for (j in 1:4) {
    phi1_bot[j]   = phi1_raw_bot[j]   * (1 - fabs(phi2_bot[j]));
    theta1_bot[j] = theta1_raw_bot[j] * (1 - fabs(theta2_bot[j]));
  }

  real phi1_A   = phi1_raw_A   * (1 - fabs(phi2_A));
  real theta1_A = theta1_raw_A * (1 - fabs(theta2_A));
  real phi1_B   = phi1_raw_B   * (1 - fabs(phi2_B));
  real theta1_B = theta1_raw_B * (1 - fabs(theta2_B));

  real phi1_top   = phi1_raw_top   * (1 - fabs(phi2_top));
  real theta1_top = theta1_raw_top * (1 - fabs(theta2_top));

  matrix[4, 4] L_Sigma = diag_pre_multiply(tau, L_Omega);
}

model {
  // ── Level 0: Grand-mean hyperpriors ──────────────────────────────────────
  mu_phi2       ~ uniform(-1, 1);
  mu_phi1_raw   ~ uniform(-1, 1);
  mu_theta2     ~ uniform(-1, 1);
  mu_theta1_raw ~ uniform(-1, 1);

  // ── Level 1: Group-mean hyperpriors ──────────────────────────────────────
  mu_phi2_A       ~ normal(mu_phi2,       0.1);
  mu_phi1_raw_A   ~ normal(mu_phi1_raw,   0.1);
  mu_theta2_A     ~ normal(mu_theta2,     0.1);
  mu_theta1_raw_A ~ normal(mu_theta1_raw, 0.1);

  mu_phi2_B       ~ normal(mu_phi2,       0.1);
  mu_phi1_raw_B   ~ normal(mu_phi1_raw,   0.1);
  mu_theta2_B     ~ normal(mu_theta2,     0.1);
  mu_theta1_raw_B ~ normal(mu_theta1_raw, 0.1);

  // ── Level 2: Series-specific parameters ──────────────────────────────────
  // Bottom level (AA, AB -> Group A; BA, BB -> Group B)
  phi2_bot[1]       ~ normal(mu_phi2_A,       0.1);
  phi1_raw_bot[1]   ~ normal(mu_phi1_raw_A,   0.1);
  theta2_bot[1]     ~ normal(mu_theta2_A,     0.1);
  theta1_raw_bot[1] ~ normal(mu_theta1_raw_A, 0.1);

  phi2_bot[2]       ~ normal(mu_phi2_A,       0.1);
  phi1_raw_bot[2]   ~ normal(mu_phi1_raw_A,   0.1);
  theta2_bot[2]     ~ normal(mu_theta2_A,     0.1);
  theta1_raw_bot[2] ~ normal(mu_theta1_raw_A, 0.1);

  phi2_bot[3]       ~ normal(mu_phi2_B,       0.1);
  phi1_raw_bot[3]   ~ normal(mu_phi1_raw_B,   0.1);
  theta2_bot[3]     ~ normal(mu_theta2_B,     0.1);
  theta1_raw_bot[3] ~ normal(mu_theta1_raw_B, 0.1);

  phi2_bot[4]       ~ normal(mu_phi2_B,       0.1);
  phi1_raw_bot[4]   ~ normal(mu_phi1_raw_B,   0.1);
  theta2_bot[4]     ~ normal(mu_theta2_B,     0.1);
  theta1_raw_bot[4] ~ normal(mu_theta1_raw_B, 0.1);

  // Mid level
  phi2_A       ~ normal(mu_phi2_A,       0.1);
  phi1_raw_A   ~ normal(mu_phi1_raw_A,   0.1);
  theta2_A     ~ normal(mu_theta2_A,     0.1);
  theta1_raw_A ~ normal(mu_theta1_raw_A, 0.1);

  phi2_B       ~ normal(mu_phi2_B,       0.1);
  phi1_raw_B   ~ normal(mu_phi1_raw_B,   0.1);
  theta2_B     ~ normal(mu_theta2_B,     0.1);
  theta1_raw_B ~ normal(mu_theta1_raw_B, 0.1);

  // Top level
  phi2_top       ~ normal(mu_phi2,       0.1);
  phi1_raw_top   ~ normal(mu_phi1_raw,   0.1);
  theta2_top     ~ normal(mu_theta2,     0.1);
  theta1_raw_top ~ normal(mu_theta1_raw, 0.1);

  // ── Variance & Covariance Components ─────────────────────────────────────
  L_Omega ~ lkj_corr_cholesky(2);
  tau     ~ cauchy(0, 1);

  sigma_A   ~ cauchy(0, 1);
  sigma_B   ~ cauchy(0, 1);
  sigma_top ~ cauchy(0, 1);
}

generated quantities {
  matrix[T, 4] Y_bottom_sim;
  vector[T]    Y_A_sim;
  vector[T]    Y_B_sim;
  vector[T]    Y_total_sim;

  {
    matrix[T, 4] e_bottom_sim;
    vector[T]    e_A_sim;
    vector[T]    e_B_sim;
    vector[T]    e_top_sim;

    // ── t = 1: e[1] = Y[1] ──────────────────────────────────────────────────
    vector[4] innov1_bot = multi_normal_cholesky_rng(rep_vector(0, 4), L_Sigma);
    e_bottom_sim[1] = innov1_bot';
    Y_bottom_sim[1] = innov1_bot';

    real eA1 = normal_rng(0, sigma_A);
    e_A_sim[1] = eA1;  Y_A_sim[1] = eA1;

    real eB1 = normal_rng(0, sigma_B);
    e_B_sim[1] = eB1;  Y_B_sim[1] = eB1;

    real et1 = normal_rng(0, sigma_top);
    e_top_sim[1] = et1;  Y_total_sim[1] = et1;

    // ── t = 2: e[2] = Y[2] - phi1*Y[1] ──────────────────────────────────────
    vector[4] innov2_bot = multi_normal_cholesky_rng(rep_vector(0, 4), L_Sigma);
    for (j in 1:4) {
      e_bottom_sim[2, j] = innov2_bot[j];
      Y_bottom_sim[2, j] = phi1_bot[j] * Y_bottom_sim[1, j] + innov2_bot[j];
    }

    real eA2 = normal_rng(0, sigma_A);
    e_A_sim[2] = eA2;  Y_A_sim[2] = phi1_A * Y_A_sim[1] + eA2;

    real eB2 = normal_rng(0, sigma_B);
    e_B_sim[2] = eB2;  Y_B_sim[2] = phi1_B * Y_B_sim[1] + eB2;

    real et2 = normal_rng(0, sigma_top);
    e_top_sim[2] = et2;  Y_total_sim[2] = phi1_top * Y_total_sim[1] + et2;

    // ── t = 3..T: Full ARMA(2,2) forward recursion ──────────────────────────
    for (t in 3:T) {
      vector[4] innov_bot = multi_normal_cholesky_rng(rep_vector(0, 4), L_Sigma);
      for (j in 1:4) {
        Y_bottom_sim[t, j] = phi1_bot[j]   * Y_bottom_sim[t-1, j]
                           + phi2_bot[j]   * Y_bottom_sim[t-2, j]
                           + innov_bot[j]
                           + theta1_bot[j] * e_bottom_sim[t-1, j]
                           + theta2_bot[j] * e_bottom_sim[t-2, j];
        e_bottom_sim[t, j] = innov_bot[j];
      }

      real innovA = normal_rng(0, sigma_A);
      Y_A_sim[t] = phi1_A * Y_A_sim[t-1] + phi2_A * Y_A_sim[t-2]
                 + innovA
                 + theta1_A * e_A_sim[t-1] + theta2_A * e_A_sim[t-2];
      e_A_sim[t] = innovA;

      real innovB = normal_rng(0, sigma_B);
      Y_B_sim[t] = phi1_B * Y_B_sim[t-1] + phi2_B * Y_B_sim[t-2]
                 + innovB
                 + theta1_B * e_B_sim[t-1] + theta2_B * e_B_sim[t-2];
      e_B_sim[t] = innovB;

      real innovTop = normal_rng(0, sigma_top);
      Y_total_sim[t] = phi1_top * Y_total_sim[t-1] + phi2_top * Y_total_sim[t-2]
                     + innovTop
                     + theta1_top * e_top_sim[t-1] + theta2_top * e_top_sim[t-2];
      e_top_sim[t] = innovTop;
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

ggsave("prior_predictive_wellspecified_hierarchical.png", fan_plot, width = 10, height = 10, dpi = 150)


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