# =============================================================================
# Bayesian Hierarchical ARMA(2,2) — rolling window evaluation (HPC version)
# =============================================================================

library(rstan)
library(MASS)
library(scoringRules)
library(dplyr)

# =============================================================================
# 1.  SETUP
# =============================================================================

args     <- commandArgs(trailingOnly = TRUE)
job_name <- args[1]

dir_accuracy <- file.path(job_name, "accuracy")
dir_stanfits <- file.path(job_name, "stanfits")

for (d in c(dir_accuracy, dir_stanfits))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

task_id <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
# task_id <-1

# =============================================================================
# 2.  LOAD DATA
# =============================================================================

synthetic_data <- read.csv("synthetic_hierarchy.csv")
data_mat       <- as.matrix(synthetic_data)   # columns: [Total, A, B, AA, AB, BA, BB]

TRAIN_SIZE <- 50  
N_ROLL     <- 500
# N_ROLL <- 1
# task_id <- 1
HORIZON    <- 5

lambda_top <- 1
lambda_mid <- 1
sigma_coh <- 0.8

IDX <- list(total = 1L, A = 2L, B = 3L, bot = 4:7)

# =============================================================================
# 3.  STAN MODEL
# =============================================================================

model_string <- "
data {
  int<lower=1> T;
  int<lower=1> H;
  matrix[T, 4] Y_bottom;
  vector[T]    Y_A;
  vector[T]    Y_B;
  vector[T]    Y_total;
  real<lower=0> lambda_mid;
  real<lower=0> lambda_top;
  real<lower=0> sigma_coh;
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
  mu_phi2       ~ uniform(-1, 1);
  mu_phi1_raw   ~ uniform(-1, 1);
  mu_theta2     ~ uniform(-1, 1);
  mu_theta1_raw ~ uniform(-1, 1);

  mu_phi2_A       ~ normal(mu_phi2,       0.14);
  mu_phi1_raw_A   ~ normal(mu_phi1_raw,   0.14);
  mu_theta2_A     ~ normal(mu_theta2,     0.14);
  mu_theta1_raw_A ~ normal(mu_theta1_raw, 0.14);

  mu_phi2_B       ~ normal(mu_phi2,       0.14);
  mu_phi1_raw_B   ~ normal(mu_phi1_raw,   0.14);
  mu_theta2_B     ~ normal(mu_theta2,     0.14);
  mu_theta1_raw_B ~ normal(mu_theta1_raw, 0.14);

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

  for (t in 3:T)
    e_bottom[t] ~ multi_normal_cholesky(rep_vector(0, 4), L_Sigma);

  for (t in 3:T) {
    target += lambda_mid * normal_lpdf(e_A[t]   | 0, sigma_A);
    target += lambda_mid * normal_lpdf(e_B[t]   | 0, sigma_B);
    target += lambda_top * normal_lpdf(e_top[t] | 0, sigma_top);
  }
  
  // coherency penalisation
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

    target += normal_lpdf(sum_bot - sum_mid | 0, sigma_coh);
    target += normal_lpdf(sum_mid - fitted_top | 0, sigma_coh);
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
# 4.  COMPILE
# =============================================================================

stanDso <- stan_model(model_code = model_string)

# =============================================================================
# 5.  10 rolls per job
# =============================================================================

rolls_per_job <- 10L
roll_start    <- (task_id - 1L) * rolls_per_job + 1L
roll_end      <- min(roll_start + rolls_per_job - 1L, N_ROLL)

for (roll in roll_start:roll_end) {
  
  message(sprintf("Starting roll %d of %d (job task %d)",
                  roll, N_ROLL, task_id))
  
  train_idx <- roll:(roll + TRAIN_SIZE - 1)
  test_idx  <- (roll + TRAIN_SIZE):(roll + TRAIN_SIZE + HORIZON - 1)
  
  train <- data_mat[train_idx, ]
  test  <- data_mat[test_idx,  ]
  
  stan_data <- list(
    T          = TRAIN_SIZE,
    H          = HORIZON,
    Y_bottom   = train[, IDX$bot],
    Y_A        = train[, IDX$A],
    Y_B        = train[, IDX$B],
    Y_total    = train[, IDX$total],
    lambda_mid = lambda_mid,
    lambda_top = lambda_top, 
    sigma_coh  = sigma_coh
  )
  
  stanFit <- sampling(
    object  = stanDso,
    data    = stan_data,
    seed    = 2025,
    chains  = 10,
    iter    = 2000,
    warmup  = 1000,
    cores   = 10,
    refresh = 0
  )
  
  saveRDS(stanFit,
          file.path(dir_stanfits, sprintf("stanfit_roll_%03d.rds", roll)))
  
  # =============================================================================
  # 6.  EXTRACT FORECASTS
  # =============================================================================
  
  fc_bot   <- rstan::extract(stanFit, "fc_bot")$fc_bot     # [draws, H, 4]
  fc_A     <- rstan::extract(stanFit, "fc_A")$fc_A         # [draws, H]
  fc_B     <- rstan::extract(stanFit, "fc_B")$fc_B         # [draws, H]
  fc_total <- rstan::extract(stanFit, "fc_total")$fc_total # [draws, H]
  
  n_draws <- nrow(fc_bot)
  
  # =============================================================================
  # 7.  ACCURACY — full suite matching tourism scripts
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
  
  compute_energy_score <- function(draws_mat, actuals, sub_idx = NULL) {
    if (!is.null(sub_idx)) draws_mat <- draws_mat[sub_idx, , drop = FALSE]
    col_sds    <- apply(draws_mat, 2, sd)
    draws_sc   <- sweep(draws_mat, 2, col_sds, "/")
    actuals_sc <- as.numeric(actuals) / col_sds
    diffs1     <- sweep(draws_sc, 2, actuals_sc, "-")
    t1         <- mean(sqrt(rowSums(diffs1^2)))
    t2         <- mean(dist(draws_sc, method = "euclidean"))
    t1 - 0.5 * t2
  }
  
  n_es_subsample <- 1000L
  set.seed(roll)
  sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)
  
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
    
    # Multivariate coherence metrics
    bot_mat      <- fc_bot[, h, ]                        # [draws x 4]
    sum_bot      <- rowSums(bot_mat)
    sum_A        <- fc_A[, h] + 0                        # copy
    sum_B        <- fc_B[, h] + 0
    sum_mid      <- fc_A[, h] + fc_B[, h]
    draws_total  <- fc_total[, h]
    
    # Energy score over all series jointly
    joint_mat     <- cbind(bot_mat, fc_A[, h], fc_B[, h], draws_total)
    actuals_joint <- c(actual_bot, actual_A, actual_B, actual_total)
    es_val        <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
    
    multivariate_df <- data.frame(
      roll             = roll,
      horizon          = h,
      es               = es_val,
      mai_bot_mid      = mean(sum_bot)   - mean(sum_mid),
      mai_mid_top      = mean(sum_mid)   - mean(draws_total),
      mai_bot_top      = mean(sum_bot)   - mean(draws_total),
      jsd_bot_mid      = JSD_metric(sum_bot,  sum_mid),
      jsd_mid_top      = JSD_metric(sum_mid,  draws_total),
      jsd_bot_top      = JSD_metric(sum_bot,  draws_total)
    )
    
    list(univariate = univariate_df, multivariate = multivariate_df)
  })
  
  roll_univariate   <- bind_rows(lapply(horizon_results, `[[`, "univariate"))
  roll_multivariate <- bind_rows(lapply(horizon_results, `[[`, "multivariate"))
  
  saveRDS(
    list(univariate = roll_univariate, multivariate = roll_multivariate),
    file.path(dir_accuracy, sprintf("accuracy_roll_%03d.rds", roll))
  )
  
  rm(stanFit); gc()
}