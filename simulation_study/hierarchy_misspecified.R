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
dir_mint_accuracy <- file.path(job_name, "accuracy_mint")
dir_stanfits <- file.path(job_name, "stanfits")

for (d in c(dir_accuracy, dir_stanfits, dir_mint_accuracy))
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
}
parameters {
  // ── Grand-mean hyperpriors (level 0) ──────────────────────────────────────
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_beta;

  // ── Group-mean hyperpriors (level 1) ──────────────────────────────────────
  real<lower=0, upper=1> mu_alpha_A;
  real<lower=0, upper=1> mu_beta_A;
  
  real<lower=0, upper=1> mu_alpha_B;
  real<lower=0, upper=1> mu_beta_B;

  // ── Bottom-level parameters: [AA, AB, BA, BB] ─────────────────────────────
  vector<lower=0, upper=1>[4] alpha_bot;
  vector<lower=0, upper=1>[4] beta_bot;

  // ── Mid-level parameters ──────────────────────────────────────────────────
  real<lower=0, upper=1> alpha_A;
  real<lower=0, upper=1> beta_A;

  real<lower=0, upper=1> alpha_B;
  real<lower=0, upper=1> beta_B;

  // ── Top-level parameters ──────────────────────────────────────────────────
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
transformed parameters {
  // ── Bottom-level States & Residuals ───────────────────────────────────────
  matrix[T, 4] l_bot;
  matrix[T, 4] b_bot;
  matrix[T, 4] e_bot;

  for (j in 1:4) {
    // Initialization (t=1)
    e_bot[1, j] = Y_bottom[1, j] - (l0_bot[j] + b0_bot[j]);
    l_bot[1, j] = l0_bot[j] + b0_bot[j] + alpha_bot[j] * e_bot[1, j];
    b_bot[1, j] = b0_bot[j] + alpha_bot[j] * beta_bot[j] * e_bot[1, j]; // Standard ETS formulation

    for (t in 2:T) {
      e_bot[t, j] = Y_bottom[t, j] - (l_bot[t-1, j] + b_bot[t-1, j]);
      l_bot[t, j] = l_bot[t-1, j] + b_bot[t-1, j] + alpha_bot[j] * e_bot[t, j];
      b_bot[t, j] = b_bot[t-1, j] + alpha_bot[j] * beta_bot[j] * e_bot[t, j];
    }
  }

  // ── Mid-level A States & Residuals ────────────────────────────────────────
  vector[T] l_A; vector[T] b_A; vector[T] e_A;
  e_A[1] = Y_A[1] - (l0_A + b0_A);
  l_A[1] = l0_A + b0_A + alpha_A * e_A[1];
  b_A[1] = b0_A + alpha_A * beta_A * e_A[1];
  for (t in 2:T) {
    e_A[t] = Y_A[t] - (l_A[t-1] + b_A[t-1]);
    l_A[t] = l_A[t-1] + b_A[t-1] + alpha_A * e_A[t];
    b_A[t] = b_A[t-1] + alpha_A * beta_A * e_A[t];
  }

  // ── Mid-level B States & Residuals ────────────────────────────────────────
  vector[T] l_B; vector[T] b_B; vector[T] e_B;
  e_B[1] = Y_B[1] - (l0_B + b0_B);
  l_B[1] = l0_B + b0_B + alpha_B * e_B[1];
  b_B[1] = b0_B + alpha_B * beta_B * e_B[1];
  for (t in 2:T) {
    e_B[t] = Y_B[t] - (l_B[t-1] + b_B[t-1]);
    l_B[t] = l_B[t-1] + b_B[t-1] + alpha_B * e_B[t];
    b_B[t] = b_B[t-1] + alpha_B * beta_B * e_B[t];
  }

  // ── Top-level States & Residuals ──────────────────────────────────────────
  vector[T] l_top; vector[T] b_top; vector[T] e_top;
  e_top[1] = Y_total[1] - (l0_top + b0_top);
  l_top[1] = l0_top + b0_top + alpha_top * e_top[1];
  b_top[1] = b0_top + alpha_top * beta_top * e_top[1];
  for (t in 2:T) {
    e_top[t] = Y_total[t] - (l_top[t-1] + b_top[t-1]);
    l_top[t] = l_top[t-1] + b_top[t-1] + alpha_top * e_top[t];
    b_top[t] = b_top[t-1] + alpha_top * beta_top * e_top[t];
  }
}
model {
  // Grand-mean Hyperpriors
  mu_alpha ~ uniform(0, 1);
  mu_beta  ~ uniform(0, 1);

  // Group-mean Hyperpriors (Level 1)
  mu_alpha_A ~ normal(mu_alpha, 0.14);
  mu_beta_A  ~ normal(mu_beta,  0.14);
  mu_alpha_B ~ normal(mu_alpha, 0.14);
  mu_beta_B  ~ normal(mu_beta,  0.14);

  // Bottom parameters
  for (j in 1:2) {
    alpha_bot[j] ~ normal(mu_alpha_A, 0.1);
    beta_bot[j]  ~ normal(mu_beta_A,  0.1);
  }
  for (j in 3:4) {
    alpha_bot[j] ~ normal(mu_alpha_B, 0.1);
    beta_bot[j]  ~ normal(mu_beta_B,  0.1);
  }

  // Mid-level parameters
  alpha_A ~ normal(mu_alpha_A, 0.1);
  beta_A  ~ normal(mu_beta_A,  0.1);
  alpha_B ~ normal(mu_alpha_B, 0.1);
  beta_B  ~ normal(mu_beta_B,  0.1);

  // Top-level parameters
  alpha_top ~ normal(mu_alpha, 0.1);
  beta_top  ~ normal(mu_beta,  0.1);

  // Initial States Priors
  l0_bot ~ normal(0, 10); b0_bot ~ normal(0, 5);
  l0_A   ~ normal(0, 10); b0_A   ~ normal(0, 5);
  l0_B   ~ normal(0, 10); b0_B   ~ normal(0, 5);
  l0_top ~ normal(0, 10); b0_top ~ normal(0, 5);

  // Scale Priors
  sigma_bottom ~ cauchy(0, 1);
  sigma_A      ~ cauchy(0, 1);
  sigma_B      ~ cauchy(0, 1);
  sigma_top    ~ cauchy(0, 1);

  // ── Laplace (Double Exponential) Likelihood ───────────────────────────────
  for (t in 1:T) {
    for (j in 1:4) {
      target += double_exponential_lpdf(e_bot[t, j] | 0, sigma_bottom[j]);
    }
    target += lambda_mid * double_exponential_lpdf(e_A[t] | 0, sigma_A);
    target += lambda_mid * double_exponential_lpdf(e_B[t] | 0, sigma_B);
    target += lambda_top * double_exponential_lpdf(e_top[t] | 0, sigma_top);
  }
}
generated quantities {
  matrix[H, 4] fc_bot;
  vector[H]    fc_A;
  vector[H]    fc_B;
  vector[H]    fc_total;

  {
    // Retrieve the final states at time T
    vector[4] l_bot_t = to_vector(l_bot[T]);
    vector[4] b_bot_t = to_vector(b_bot[T]);
    real l_A_t = l_A[T]; real b_A_t = b_A[T];
    real l_B_t = l_B[T]; real b_B_t = b_B[T];
    real l_top_t = l_top[T]; real b_top_t = b_top[T];

    for (h in 1:H) {
      // Bottom level projections
      for (j in 1:4) {
        real innov_bot = double_exponential_rng(0, sigma_bottom[j]);
        fc_bot[h, j]   = l_bot_t[j] + b_bot_t[j] + innov_bot;
        
        // Update states for next step forecasting
        l_bot_t[j] = l_bot_t[j] + b_bot_t[j] + alpha_bot[j] * innov_bot;
        b_bot_t[j] = b_bot_t[j] + alpha_bot[j] * beta_bot[j] * innov_bot;
      }

      // Mid level A projections
      real innov_A = double_exponential_rng(0, sigma_A);
      fc_A[h] = l_A_t + b_A_t + innov_A;
      l_A_t   = l_A_t + b_A_t + alpha_A * innov_A;
      b_A_t   = b_A_t + alpha_A * beta_A * innov_A;

      // Mid level B projections
      real innov_B = double_exponential_rng(0, sigma_B);
      fc_B[h] = l_B_t + b_B_t + innov_B;
      l_B_t   = l_B_t + b_B_t + alpha_B * innov_B;
      b_B_t   = b_B_t + alpha_B * beta_B * innov_B;

      // Top level projections
      real innov_top = double_exponential_rng(0, sigma_top);
      fc_total[h]    = l_top_t + b_top_t + innov_top;
      l_top_t        = l_top_t + b_top_t + alpha_top * innov_top;
      b_top_t        = b_top_t + alpha_top * beta_top * innov_top;
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
    lambda_top = lambda_top
  )
  
  stanFit <- sampling(
    object  = stanDso,
    data    = stan_data,
    seed    = 2025,
    chains  = 10,
    iter    = 200,
    warmup  = 100,
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
  
  n_draws <- dim(fc_bot)[1]
  
  # =============================================================================
  # 6b. MinT(Shrink) RECONCILIATION
  # =============================================================================
  # S matrix: [7 x 4], ordering [Total, A, B, AA, AB, BA, BB]
  S_mat <- rbind(
    c(1, 1, 1, 1),
    c(1, 1, 0, 0),
    c(0, 0, 1, 1),
    c(1, 0, 0, 0),
    c(0, 1, 0, 0),
    c(0, 0, 1, 0),
    c(0, 0, 0, 1)
  )
  
  # ── Build in-sample residual matrix from Stan's transformed parameters ────────
  # e_bottom: [T x 4], e_A: [T], e_B: [T], e_top: [T] — posterior means used
  # for W estimation (using posterior means keeps this computationally tractable)
  
  e_bot_pm  <- summary(stanFit, pars = "e_bot")$summary[, "mean"]
  e_A_pm    <- summary(stanFit, pars = "e_A")$summary[,  "mean"]
  e_B_pm    <- summary(stanFit, pars = "e_B")$summary[,  "mean"]
  e_top_pm  <- summary(stanFit, pars = "e_top")$summary[, "mean"]
  
  # e_bot_pm is flattened [T*4] in column-major order from Stan
  e_bot_mat <- matrix(e_bot_pm, nrow = TRAIN_SIZE, ncol = 4)  # [T x 4]
  
  # Assemble full residual matrix: [T x 7], order matches S rows
  E_full <- cbind(
    e_top_pm,           # Total
    e_A_pm,             # A
    e_B_pm,             # B
    e_bot_mat           # AA, AB, BA, BB
  )
  # Drop first 2 rows (ARMA initialisation is unreliable there)
  E_full <- E_full[3:TRAIN_SIZE, ]
  T_eff  <- nrow(E_full)  # TRAIN_SIZE - 2
  
  # ── MinT(Shrink) covariance estimator (Ledoit-Wolf-style shrinkage) ───────────
  # Shrinks the sample covariance toward a scaled identity, following
  # Wickramasuriya et al. (2019) Appendix A.
  mint_shrink_W <- function(E) {
    m     <- ncol(E)
    T_e   <- nrow(E)
    S_cov <- cov(E)                          # [m x m] sample covariance
    
    # Target: diagonal of sample covariance (preserves scale, kills off-diagonal)
    diag_S <- diag(diag(S_cov))
    
    # Ledoit-Wolf optimal shrinkage intensity (Oracle-free analytical formula)
    # as used in the hts / FoReco implementation
    sum_sq_offdiag <- sum((S_cov - diag_S)^2)
    denom          <- sum(S_cov^2) - sum(diag(S_cov)^2)
    
    if (denom < .Machine$double.eps) {
      lambda_star <- 1   # already diagonal
    } else {
      lambda_star <- min(1, max(0, (sum_sq_offdiag / denom) * ((T_e + 2) / T_e)))
    }
    
    W_shrink <- (1 - lambda_star) * S_cov + lambda_star * diag_S
    list(W = W_shrink, lambda = lambda_star)
  }
  
  W_res      <- mint_shrink_W(E_full)
  W_hat      <- W_res$W
  lambda_lw  <- W_res$lambda
  message(sprintf("  MinT shrinkage lambda = %.4f", lambda_lw))
  
  # ── MinT projection matrix: P = (S'W⁻¹S)⁻¹ S'W⁻¹  ──────────────────────────
  # Reconciled forecasts: y_tilde = S P y_hat  where y_hat is [7 x 1]
  W_inv    <- solve(W_hat)
  StWinv   <- t(S_mat) %*% W_inv                        # [4 x 7]
  StWinvS  <- StWinv %*% S_mat                          # [4 x 4]
  P_mint   <- solve(StWinvS) %*% StWinv                 # [4 x 7]
  SP_mint  <- S_mat %*% P_mint                          # [7 x 7]  full projection
  
  # Sanity check: SP should map coherent vectors to themselves
  coherence_err <- max(abs(SP_mint %*% S_mat[, 1] - S_mat[, 1]))
  if (coherence_err > 1e-8)
    warning("MinT projection sanity check failed; check S_mat and W_hat.")
  
  # ── Apply P to every posterior draw at every horizon ─────────────────────────
  # fc_all_draws: list of length H, each element [n_draws x 7]
  reconcile_draws <- function(h) {
    # Assemble [n_draws x 7] base forecast matrix
    y_hat_draws <- cbind(
      fc_total[, h],       # Total
      fc_A[, h],           # A
      fc_B[, h],           # B
      fc_bot[, h, 1],      # AA
      fc_bot[, h, 2],      # AB
      fc_bot[, h, 3],      # BA
      fc_bot[, h, 4]       # BB
    )  # [n_draws x 7]
    
    # Apply SP: each row is a 7-vector of base forecasts
    # y_tilde = y_hat %*% t(SP)  (broadcasting across draws)
    y_tilde_draws <- y_hat_draws %*% t(SP_mint)   # [n_draws x 7]
    colnames(y_tilde_draws) <- c("Total", "A", "B", "AA", "AB", "BA", "BB")
    y_tilde_draws
  }
  
  fc_reconciled <- lapply(seq_len(HORIZON), reconcile_draws)
  # fc_reconciled[[h]] is [n_draws x 7], coherent by construction
  
  # =============================================================================
  # 7.  ACCURACY — base + reconciled, saved separately
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
    # Apply subsampling if requested for computational tractability
    if (!is.null(sub_idx)) {
      draws_mat <- draws_mat[sub_idx, , drop = FALSE]
    }
    
    # Ensure actuals is a numeric vector matching the columns of draws_mat
    actuals <- as.numeric(actuals)
    
    # Term 1: Mean Euclidean distance from each posterior draw to the actuals
    diffs1 <- sweep(draws_mat, 2, actuals, "-")
    t1     <- mean(sqrt(rowSums(diffs1^2)))
    
    # Term 2: Mean pairwise Euclidean distance between all posterior draws
    t2     <- mean(dist(draws_mat, method = "euclidean"))
    
    # Energy Score formula
    return(t1 - 0.5 * t2)
  }
  n_es_subsample <- 1000L
  set.seed(roll)
  sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)
  
  bot_names <- c("AA", "AB", "BA", "BB")
  
  # ---------------------------------------------------------------------------
  # Shared accuracy helper — accepts either base or reconciled forecast arrays.
  # rec_h_list: list of length H, each element a [n_draws x 7] matrix with
  #   columns ordered [Total, A, B, AA, AB, BA, BB].
  # ---------------------------------------------------------------------------
  make_accuracy <- function(fc_bot_arr, fc_A_mat, fc_B_mat, fc_total_mat,
                            label_roll) {
    
    make_row <- function(draws, actual, series, level) {
      data.frame(
        roll       = label_roll,
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
          make_row(fc_bot_arr[, h, j], actual_bot[j], bot_names[j], "bottom")),
        make_row(fc_A_mat[, h],     actual_A,     "A",     "mid"),
        make_row(fc_B_mat[, h],     actual_B,     "B",     "mid"),
        make_row(fc_total_mat[, h], actual_total, "Total", "top")
      ) %>% mutate(horizon = h, .after = roll)
      
      # Multivariate coherence metrics
      bot_mat      <- fc_bot_arr[, h, ]
      sum_bot      <- rowSums(bot_mat)
      sum_mid      <- fc_A_mat[, h] + fc_B_mat[, h]
      draws_total  <- fc_total_mat[, h]
      
      # Explicitly align columns to match S_mat: [Total, A, B, AA, AB, BA, BB]
      joint_mat <- cbind(
        draws_total,         # Column 1: Total
        fc_A_mat[, h],       # Column 2: A
        fc_B_mat[, h],       # Column 3: B
        bot_mat[, 1],        # Column 4: AA
        bot_mat[, 2],        # Column 5: AB
        bot_mat[, 3],        # Column 6: BA
        bot_mat[, 4]         # Column 7: BB
      )
      
      # Align target vector to match the same coordinate mapping space
      actuals_joint <- c(actual_total, actual_A, actual_B, actual_bot)
      es_val        <- compute_energy_score(joint_mat, actuals_joint, sub_idx)
      
      multivariate_df <- data.frame(
        roll             = label_roll,
        horizon          = h,
        es               = es_val,
        mai_bot_mid      = mean(sum_bot) - mean(sum_mid),
        mai_mid_top      = mean(sum_mid) - mean(draws_total),
        mai_bot_top      = mean(sum_bot) - mean(draws_total),
        jsd_bot_mid      = JSD_metric(sum_bot,  sum_mid),
        jsd_mid_top      = JSD_metric(sum_mid,  draws_total),
        jsd_bot_top      = JSD_metric(sum_bot,  draws_total)
      )
      
      list(univariate = univariate_df, multivariate = multivariate_df)
    })
    
    list(
      univariate   = bind_rows(lapply(horizon_results, `[[`, "univariate")),
      multivariate = bind_rows(lapply(horizon_results, `[[`, "multivariate"))
    )
  }
  
  # ---------------------------------------------------------------------------
  # Unpack reconciled draws back into separate arrays matching base structure
  # fc_reconciled[[h]] is [n_draws x 7]: [Total, A, B, AA, AB, BA, BB]
  # ---------------------------------------------------------------------------
  fc_bot_mint    <- array(NA_real_, dim = c(n_draws, HORIZON, 4))
  fc_A_mint      <- matrix(NA_real_, nrow = n_draws, ncol = HORIZON)
  fc_B_mint      <- matrix(NA_real_, nrow = n_draws, ncol = HORIZON)
  fc_total_mint  <- matrix(NA_real_, nrow = n_draws, ncol = HORIZON)
  
  for (h in seq_len(HORIZON)) {
    rec_h                <- fc_reconciled[[h]]   # [n_draws x 7]
    fc_total_mint[, h]   <- rec_h[, 1]
    fc_A_mint[, h]       <- rec_h[, 2]
    fc_B_mint[, h]       <- rec_h[, 3]
    fc_bot_mint[, h, ]   <- rec_h[, 4:7]
  }
  
  # ---------------------------------------------------------------------------
  # Compute accuracy for both variants and save to separate directories
  # ---------------------------------------------------------------------------
  base_acc <- make_accuracy(fc_bot,         fc_A,       fc_B,       fc_total,
                            label_roll = roll)
  
  mint_acc <- make_accuracy(fc_bot_mint,    fc_A_mint,  fc_B_mint,  fc_total_mint,
                            label_roll = roll)
  
  saveRDS(
    base_acc,
    file.path(dir_accuracy,      sprintf("accuracy_roll_%03d.rds", roll))
  )
  
  saveRDS(
    list(univariate   = mint_acc$univariate,
         multivariate = mint_acc$multivariate,
         mint_lambda  = lambda_lw),
    file.path(dir_mint_accuracy, sprintf("accuracy_roll_%03d.rds", roll))
  )
  
  rm(stanFit); gc()
}