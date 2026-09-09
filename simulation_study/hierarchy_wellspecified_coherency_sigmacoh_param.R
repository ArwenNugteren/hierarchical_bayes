library(rstan)
library(MASS)
library(scoringRules)
library(dplyr)
library(ggplot2)

# =============================================================================
# 1.  SETUP
# =============================================================================

synthetic_data <- read.csv("synthetic_hierarchy.csv")
data_mat       <- as.matrix(synthetic_data)  # [Total, A, B, AA, AB, BA, BB]

TRAIN_SIZE <- 50
HORIZON    <- 5

lambda_top <- 1
lambda_mid <- 1
lambda_bot <- 1

# sigma_coh prior: TruncNormal(sigma_coh_0, sigma_coh_0, lower=0)
# [1] = bottom <-> mid  [2] = mid <-> top
sigma_coh_0 <- sqrt(0.5)

IDX <- list(total = 1L, A = 2L, B = 3L, bot = 4:7)

# Single window: roll = 1
roll      <- 1L
train_idx <- roll:(roll + TRAIN_SIZE - 1)
test_idx  <- (roll + TRAIN_SIZE):(roll + TRAIN_SIZE + HORIZON - 1)
train     <- data_mat[train_idx, ]
test      <- data_mat[test_idx,  ]

# =============================================================================
# 2.  STAN MODEL
# =============================================================================

model_string <- "
data {
  int<lower=1> T;
  int<lower=1> H;
  matrix[T, 4] Y_bottom;
  vector[T]    Y_A;
  vector[T]    Y_B;
  vector[T]    Y_total;
  real<lower=0> lambda_bot;
  real<lower=0> lambda_mid;
  real<lower=0> lambda_top;
  real<lower=0> sigma_coh_0;   // shared prior mean and SD for both sigma_coh elements
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

  // ── Learnt coherence penalty SDs ─────────────────────────────────────────
  // [1] bottom <-> mid   [2] mid <-> top
  vector<lower=0>[2] sigma_coh;
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

  // ── ARMA residuals ────────────────────────────────────────────────────────
  matrix[T, 4] e_bottom;
  for (j in 1:4) {
    e_bottom[1, j] = Y_bottom[1, j];
    e_bottom[2, j] = Y_bottom[2, j] - phi1_bot[j] * Y_bottom[1, j];
    for (t in 3:T)
      e_bottom[t, j] = Y_bottom[t, j]
                     - phi1_bot[j]   * Y_bottom[t-1, j]
                     - phi2_bot[j]   * Y_bottom[t-2, j]
                     - theta1_bot[j] * e_bottom[t-1, j]
                     - theta2_bot[j] * e_bottom[t-2, j];
  }

  vector[T] e_A;
  e_A[1] = Y_A[1];
  e_A[2] = Y_A[2] - phi1_A * Y_A[1];
  for (t in 3:T)
    e_A[t] = Y_A[t]
           - phi1_A * Y_A[t-1] - phi2_A * Y_A[t-2]
           - theta1_A * e_A[t-1] - theta2_A * e_A[t-2];

  vector[T] e_B;
  e_B[1] = Y_B[1];
  e_B[2] = Y_B[2] - phi1_B * Y_B[1];
  for (t in 3:T)
    e_B[t] = Y_B[t]
           - phi1_B * Y_B[t-1] - phi2_B * Y_B[t-2]
           - theta1_B * e_B[t-1] - theta2_B * e_B[t-2];

  vector[T] e_top;
  e_top[1] = Y_total[1];
  e_top[2] = Y_total[2] - phi1_top * Y_total[1];
  for (t in 3:T)
    e_top[t] = Y_total[t]
             - phi1_top * Y_total[t-1] - phi2_top * Y_total[t-2]
             - theta1_top * e_top[t-1] - theta2_top * e_top[t-2];
}
model {
  // ── Hyperpriors ───────────────────────────────────────────────────────────
  mu_phi2       ~ uniform(-1, 1);
  mu_phi1_raw   ~ uniform(-1, 1);
  mu_theta2     ~ uniform(-1, 1);
  mu_theta1_raw ~ uniform(-1, 1);

  mu_phi2_A       ~ normal(mu_phi2,       0.1);
  mu_phi1_raw_A   ~ normal(mu_phi1_raw,   0.1);
  mu_theta2_A     ~ normal(mu_theta2,     0.1);
  mu_theta1_raw_A ~ normal(mu_theta1_raw, 0.1);

  mu_phi2_B       ~ normal(mu_phi2,       0.1);
  mu_phi1_raw_B   ~ normal(mu_phi1_raw,   0.1);
  mu_theta2_B     ~ normal(mu_theta2,     0.1);
  mu_theta1_raw_B ~ normal(mu_theta1_raw, 0.1);

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

  phi2_A       ~ normal(mu_phi2_A,       0.1);
  phi1_raw_A   ~ normal(mu_phi1_raw_A,   0.1);
  theta2_A     ~ normal(mu_theta2_A,     0.1);
  theta1_raw_A ~ normal(mu_theta1_raw_A, 0.1);

  phi2_B       ~ normal(mu_phi2_B,       0.1);
  phi1_raw_B   ~ normal(mu_phi1_raw_B,   0.1);
  theta2_B     ~ normal(mu_theta2_B,     0.1);
  theta1_raw_B ~ normal(mu_theta1_raw_B, 0.1);

  phi2_top       ~ normal(mu_phi2,       0.1);
  phi1_raw_top   ~ normal(mu_phi1_raw,   0.1);
  theta2_top     ~ normal(mu_theta2,     0.1);
  theta1_raw_top ~ normal(mu_theta1_raw, 0.1);

  L_Omega ~ lkj_corr_cholesky(2);
  tau     ~ cauchy(0, 1);

  sigma_A   ~ cauchy(0, 1);
  sigma_B   ~ cauchy(0, 1);
  sigma_top ~ cauchy(0, 1);

  // ── Learnt sigma_coh: TruncNormal(sigma_coh_0, sigma_coh_0, lower=0) ─────
  // [1] bottom <-> mid   [2] mid <-> top
  for (k in 1:2)
    sigma_coh[k] ~ normal(sigma_coh_0, sigma_coh_0) T[0, ];

  // ── Likelihoods ───────────────────────────────────────────────────────────
  for (t in 3:T)
    target += lambda_bot * multi_normal_cholesky_lpdf(
                e_bottom[t] | rep_vector(0, 4), L_Sigma);

  for (t in 3:T) {
    target += lambda_mid * normal_lpdf(e_A[t]   | 0, sigma_A);
    target += lambda_mid * normal_lpdf(e_B[t]   | 0, sigma_B);
    target += lambda_top * normal_lpdf(e_top[t] | 0, sigma_top);
  }

  // ── Coherence penalties on fitted values (level-specific sigma_coh) ───────
  for (t in 3:T) {
    real fitted_AA  = Y_bottom[t, 1] - e_bottom[t, 1];
    real fitted_AB  = Y_bottom[t, 2] - e_bottom[t, 2];
    real fitted_BA  = Y_bottom[t, 3] - e_bottom[t, 3];
    real fitted_BB  = Y_bottom[t, 4] - e_bottom[t, 4];
    real fitted_A   = Y_A[t]         - e_A[t];
    real fitted_B   = Y_B[t]         - e_B[t];
    real fitted_top = Y_total[t]     - e_top[t];

    real sum_bot = fitted_AA + fitted_AB + fitted_BA + fitted_BB;
    real sum_mid = fitted_A  + fitted_B;

    target += normal_lpdf(sum_bot - sum_mid    | 0, sigma_coh[1]);
    target += normal_lpdf(sum_mid - fitted_top | 0, sigma_coh[2]);
  }
}
generated quantities {
  matrix[H, 4] fc_bot;
  vector[H]    fc_A;
  vector[H]    fc_B;
  vector[H]    fc_total;

  {
    vector[4] y_tm1 = to_vector(Y_bottom[T]);
    vector[4] y_tm2 = to_vector(Y_bottom[T-1]);
    vector[4] e_tm1 = to_vector(e_bottom[T]);
    vector[4] e_tm2 = to_vector(e_bottom[T-1]);

    real yA_tm1 = Y_A[T];      real yA_tm2 = Y_A[T-1];
    real eA_tm1 = e_A[T];      real eA_tm2 = e_A[T-1];

    real yB_tm1 = Y_B[T];      real yB_tm2 = Y_B[T-1];
    real eB_tm1 = e_B[T];      real eB_tm2 = e_B[T-1];

    real yt_tm1 = Y_total[T];  real yt_tm2 = Y_total[T-1];
    real et_tm1 = e_top[T];    real et_tm2 = e_top[T-1];

    for (h in 1:H) {
      vector[4] innov_bot = multi_normal_cholesky_rng(rep_vector(0, 4), L_Sigma);

      vector[4] y_h;
      for (j in 1:4)
        y_h[j] = phi1_bot[j] * y_tm1[j] + phi2_bot[j] * y_tm2[j]
               + innov_bot[j]
               + theta1_bot[j] * e_tm1[j] + theta2_bot[j] * e_tm2[j];
      fc_bot[h] = y_h';
      e_tm2 = e_tm1;  e_tm1 = innov_bot;
      y_tm2 = y_tm1;  y_tm1 = y_h;

      real innov_A = normal_rng(0, sigma_A);
      real yA_h    = phi1_A * yA_tm1 + phi2_A * yA_tm2
                   + innov_A
                   + theta1_A * eA_tm1 + theta2_A * eA_tm2;
      fc_A[h] = yA_h;
      eA_tm2 = eA_tm1;  eA_tm1 = innov_A;
      yA_tm2 = yA_tm1;  yA_tm1 = yA_h;

      real innov_B = normal_rng(0, sigma_B);
      real yB_h    = phi1_B * yB_tm1 + phi2_B * yB_tm2
                   + innov_B
                   + theta1_B * eB_tm1 + theta2_B * eB_tm2;
      fc_B[h] = yB_h;
      eB_tm2 = eB_tm1;  eB_tm1 = innov_B;
      yB_tm2 = yB_tm1;  yB_tm1 = yB_h;

      real innov_top = normal_rng(0, sigma_top);
      real yt_h      = phi1_top * yt_tm1 + phi2_top * yt_tm2
                     + innov_top
                     + theta1_top * et_tm1 + theta2_top * et_tm2;
      fc_total[h] = yt_h;
      et_tm2 = et_tm1;  et_tm1 = innov_top;
      yt_tm2 = yt_tm1;  yt_tm1 = yt_h;
    }
  }
}
"

# =============================================================================
# 3.  COMPILE AND FIT
# =============================================================================

stanDso <- stan_model(model_code = model_string)

stan_data <- list(
  T           = TRAIN_SIZE,
  H           = HORIZON,
  Y_bottom    = train[, IDX$bot],
  Y_A         = train[, IDX$A],
  Y_B         = train[, IDX$B],
  Y_total     = train[, IDX$total],
  lambda_bot  = lambda_bot,
  lambda_mid  = lambda_mid,
  lambda_top  = lambda_top,
  sigma_coh_0 = sigma_coh_0
)

start_time <- Sys.time()

stanFit <- sampling(
  object  = stanDso,
  data    = stan_data,
  seed    = 2025,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  cores   = 4,
  refresh = 100, 
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

cat("Elapsed:", format(Sys.time() - start_time), "\n")
cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")

# =============================================================================
# 4.  EXTRACT sigma_coh POSTERIOR AND PLOT
# =============================================================================

sigma_coh_draws <- rstan::extract(stanFit, pars = "sigma_coh")$sigma_coh
# sigma_coh_draws: [n_draws x 2]
# [1] = bottom <-> mid   [2] = mid <-> top

level_pair_labels <- c(
  "Bottom (AA+AB+BA+BB) <-> Mid (A+B)",
  "Mid (A+B) <-> Top (Total)"
)
colnames(sigma_coh_draws) <- level_pair_labels

# Truncated normal density (lower = 0)
dtruncnorm <- function(x, mu, sigma) {
  nc <- 1 - pnorm(0, mean = mu, sd = sigma)
  ifelse(x < 0, 0, dnorm(x, mean = mu, sd = sigma) / nc)
}

pdf("sigma_coh_posterior_arma.pdf", width = 8, height = 5)

for (k in 1:2) {
  draws_k <- sigma_coh_draws[, k]
  
  x_max <- max(quantile(draws_k, 0.995), 3 * sigma_coh_0)
  x_seq <- seq(0, x_max, length.out = 1000)
  
  prior_df   <- data.frame(x = x_seq,
                           y = dtruncnorm(x_seq, sigma_coh_0, sigma_coh_0))
  post_df    <- data.frame(draws = draws_k)
  
  p <- ggplot() +
    geom_histogram(data = post_df, aes(x = draws, y = after_stat(density)),
                   bins = 60, fill = "#4C72B0", colour = "white", alpha = 0.75) +
    geom_density(data = post_df, aes(x = draws),
                 colour = "#1a3f6f", linewidth = 0.9) +
    geom_line(data = prior_df, aes(x = x, y = y),
              colour = "firebrick", linewidth = 0.9, linetype = "dashed") +
    geom_vline(xintercept = sigma_coh_0,
               colour = "firebrick", linewidth = 0.7, linetype = "dotted") +
    geom_vline(xintercept = median(draws_k),
               colour = "#1a3f6f", linewidth = 0.7, linetype = "dotted") +
    scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = paste0("sigma_coh posterior: ", level_pair_labels[k]),
      subtitle = sprintf(
        "Prior: TruncNormal(mu=%.3f, sd=%.3f)   |   Prior mean (dotted red): %.3f   |   Posterior median (dotted blue): %.3f",
        sigma_coh_0, sigma_coh_0, sigma_coh_0, median(draws_k)
      ),
      x = expression(sigma[coh]),
      y = "Density"
    ) +
    annotate("text", x = x_max * 0.98, y = Inf,
             label = sprintf(
               "Post. mean = %.4f\nPost. sd   = %.4f\nPost. 95%% CI = [%.4f, %.4f]",
               mean(draws_k), sd(draws_k),
               quantile(draws_k, 0.025), quantile(draws_k, 0.975)
             ),
             hjust = 1, vjust = 1.3, size = 3.2, colour = "grey20") +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 7.5, colour = "grey30")
    )
  
  print(p)
}

dev.off()
cat("sigma_coh plots saved to sigma_coh_posterior_arma.pdf\n")

# =============================================================================
# 5.  EXTRACT FORECASTS (single pass)
# =============================================================================

all_fc   <- rstan::extract(stanFit, pars = c("fc_bot", "fc_A", "fc_B", "fc_total"))
fc_bot   <- all_fc$fc_bot    # [draws, H, 4]
fc_A     <- all_fc$fc_A      # [draws, H]
fc_B     <- all_fc$fc_B      # [draws, H]
fc_total <- all_fc$fc_total  # [draws, H]
rm(all_fc)

n_draws <- nrow(fc_bot)

# =============================================================================
# 6.  ACCURACY
# =============================================================================

JSD_metric <- function(p_draws, q_draws, n_grid = 1000) {
  lo    <- min(c(p_draws, q_draws)) - 3 * sd(c(p_draws, q_draws))
  hi    <- max(c(p_draws, q_draws)) + 3 * sd(c(p_draws, q_draws))
  p_kde <- pmax(density(p_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  q_kde <- pmax(density(q_draws, from = lo, to = hi, n = n_grid)$y, 1e-300)
  m_kde <- 0.5 * (p_kde + q_kde)
  dx    <- (hi - lo) / (n_grid - 1)
  jsd   <- 0.5 * (sum(p_kde * log(p_kde / m_kde)) +
                    sum(q_kde * log(q_kde / m_kde))) * dx
  sqrt(max(0, min(log(2), jsd)))
}

# Panagiotelis-style: unscaled, two independent subsamples for pairwise term
compute_energy_score <- function(draws_mat, actuals, n_subsample = NULL) {
  if (!is.null(n_subsample)) {
    n    <- min(n_subsample, nrow(draws_mat))
    idx1 <- sample.int(nrow(draws_mat), size = n, replace = FALSE)
    idx2 <- sample.int(nrow(draws_mat), size = n, replace = FALSE)
    x    <- draws_mat[idx1, , drop = FALSE]
    xs   <- draws_mat[idx2, , drop = FALSE]
  } else {
    x    <- draws_mat
    xs   <- draws_mat
  }
  diffs <- sweep(x, 2, as.numeric(actuals), "-")
  t1    <- mean(sqrt(rowSums(diffs^2)))
  t2    <- mean(sqrt(rowSums((x - xs)^2)))
  t1 - 0.5 * t2
}

n_es_subsample <- 1000L
set.seed(roll)

bot_names <- c("AA", "AB", "BA", "BB")

make_row <- function(draws, actual, series, level) {
  data.frame(
    roll       = roll,
    series     = series,
    level      = level,
    actual     = actual,
    fc_mean    = mean(draws),
    ae         = abs(mean(draws) - actual),
    se         = (mean(draws) - actual)^2,
    crps       = scoringRules::crps_sample(actual, draws),
    logS       = scoringRules::logs_sample(actual, draws),
    covered_90 = as.integer(actual >= quantile(draws, 0.05) &
                              actual <= quantile(draws, 0.95)),
    covered_95 = as.integer(actual >= quantile(draws, 0.025) &
                              actual <= quantile(draws, 0.975)),
    pit        = mean(draws < actual)
  )
}

horizon_results <- lapply(seq_len(HORIZON), function(h) {
  
  actual_bot   <- test[h, IDX$bot]
  actual_A     <- test[h, IDX$A]
  actual_B     <- test[h, IDX$B]
  actual_total <- test[h, IDX$total]
  
  univariate_df <- bind_rows(
    lapply(1:4, function(j)
      make_row(fc_bot[, h, j], actual_bot[j], bot_names[j], "bottom")),
    make_row(fc_A[, h],     actual_A,     "A",     "mid"),
    make_row(fc_B[, h],     actual_B,     "B",     "mid"),
    make_row(fc_total[, h], actual_total, "Total", "top")
  ) %>% mutate(horizon = h, .after = roll)
  
  bot_mat     <- fc_bot[, h, ]
  sum_bot     <- rowSums(bot_mat)
  sum_mid     <- fc_A[, h] + fc_B[, h]
  draws_total <- fc_total[, h]
  
  joint_mat <- cbind(
    draws_total,
    fc_A[, h],
    fc_B[, h],
    bot_mat[, 1],
    bot_mat[, 2],
    bot_mat[, 3],
    bot_mat[, 4]
  )
  actuals_joint <- c(actual_total, actual_A, actual_B, actual_bot)
  
  multivariate_df <- data.frame(
    roll        = roll,
    horizon     = h,
    es          = compute_energy_score(joint_mat, actuals_joint,
                                       n_subsample = n_es_subsample),
    mai_bot_mid = mean(sum_bot) - mean(sum_mid),
    mai_mid_top = mean(sum_mid) - mean(draws_total),
    mai_bot_top = mean(sum_bot) - mean(draws_total),
    jsd_bot_mid = JSD_metric(sum_bot,  sum_mid),
    jsd_mid_top = JSD_metric(sum_mid,  draws_total),
    jsd_bot_top = JSD_metric(sum_bot,  draws_total)
  )
  
  list(univariate = univariate_df, multivariate = multivariate_df)
})

roll_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
roll_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))

# =============================================================================
# 7.  SUMMARY TABLE
# =============================================================================

make_level_row <- function(lev) {
  roll_univariate %>%
    filter(level == lev) %>%
    summarise(
      level      = lev,
      MAE        = mean(ae),
      RMSE       = sqrt(mean(se)),
      CRPS       = mean(crps),
      LogS       = mean(logS),
      Cover_90   = mean(covered_90),
      Cover_95   = mean(covered_95),
      ES          = NA_real_,
      MAI_bot_mid = NA_real_,
      MAI_mid_top = NA_real_,
      MAI_bot_top = NA_real_,
      JSD_bot_mid = NA_real_,
      JSD_mid_top = NA_real_,
      JSD_bot_top = NA_real_
    )
}

overall_row <- bind_cols(
  roll_univariate %>%
    summarise(
      level    = "Overall",
      MAE      = mean(ae),
      RMSE     = sqrt(mean(se)),
      CRPS     = mean(crps),
      LogS     = mean(logS),
      Cover_90 = mean(covered_90),
      Cover_95 = mean(covered_95)
    ),
  roll_multivariate %>%
    summarise(
      ES          = mean(es),
      MAI_bot_mid = mean(abs(mai_bot_mid)),
      MAI_mid_top = mean(abs(mai_mid_top)),
      MAI_bot_top = mean(abs(mai_bot_top)),
      JSD_bot_mid = mean(jsd_bot_mid, na.rm = TRUE),
      JSD_mid_top = mean(jsd_mid_top, na.rm = TRUE),
      JSD_bot_top = mean(jsd_bot_top, na.rm = TRUE)
    )
)

summary_table <- bind_rows(
  make_level_row("bottom"),
  make_level_row("mid"),
  make_level_row("top"),
  overall_row
)

print(summary_table)

# =============================================================================
# 8.  SAVE
# =============================================================================

save(stanFit, summary_table, roll_univariate, roll_multivariate,
     sigma_coh_draws, sigma_coh_0, level_pair_labels,
     file = "arma_sigma_coh_single_window.RData")

cat("Results saved to arma_sigma_coh_single_window.RData\n")