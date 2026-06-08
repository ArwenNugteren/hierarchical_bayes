library(rstan)
library(dplyr)

# =============================================================================
# 1. SETUP & DIRECTORIES
# =============================================================================
synthetic_data <- read.csv("synthetic_hierarchy.csv")
data_mat       <- as.matrix(synthetic_data)   # columns: [Total, A, B, AA, AB, BA, BB]

job_name <- "synthetic_hierarchy_22208112"  # Change this to your target directory

dir_stanfits      <- file.path(job_name, "stanfits")
dir_output_base   <- file.path(job_name, "post_hoc_es_base")
dir_output_mint   <- file.path(job_name, "post_hoc_es_mint")

for (d in c(dir_output_base, dir_output_mint)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Changed to 1:500 to match your actual file names on disk
origins <- 1:500  

TRAIN_SIZE <- 50  
N_ROLL     <- 500
HORIZON    <- 5

IDX <- list(total = 1L, A = 2L, B = 3L, bot = 4:7)

n_es_subsample <- 1000L

# S matrix mapping: [Total, A, B, AA, AB, BA, BB]
S_mat <- rbind(
  c(1, 1, 1, 1), # Total
  c(1, 1, 0, 0), # A
  c(0, 0, 1, 1), # B
  c(1, 0, 0, 0), # AA
  c(0, 1, 0, 0), # AB
  c(0, 0, 1, 0), # BA
  c(0, 0, 0, 1)  # BB
)

# =============================================================================
# 2. COMPARABLE ENERGY SCORE FUNCTION
# =============================================================================
compute_comparable_energy_score <- function(draws_mat, actuals, scale_vec, sub_idx = NULL) {
  if (!is.null(sub_idx)) {
    draws_mat <- draws_mat[sub_idx, , drop = FALSE]
  }
  
  # Coordinate space transformation using fixed baseline scales
  draws_sc   <- sweep(draws_mat, 2, scale_vec, "/")
  actuals_sc <- as.numeric(actuals) / scale_vec
  
  # Term 1: Mean Euclidean distance to actuals
  diffs1 <- sweep(draws_sc, 2, actuals_sc, "-")
  t1     <- mean(sqrt(rowSums(diffs1^2)))
  
  # Term 2: Mean pairwise distance between draws
  pw     <- dist(draws_sc, method = "euclidean")
  t2     <- mean(pw)
  
  return(t1 - 0.5 * t2)
}

# =============================================================================
# 3. PROCESSING LOOP
# =============================================================================
for (roll in origins) {
  
  stanfit_path <- file.path(dir_stanfits, sprintf("stanfit_roll_%03d.rds", roll))
  
  if (!file.exists(stanfit_path)) {
    warning(sprintf("Stanfit file missing for roll index %d, skipping.", roll))
    next
  }
  
  # Back to original logic: roll = 1 reads data rows 1:50 and tests 51:55
  train_idx <- roll:(roll + TRAIN_SIZE - 1)
  test_idx  <- (roll + TRAIN_SIZE):(roll + TRAIN_SIZE + HORIZON - 1)
  
  train <- data_mat[train_idx, ]
  test  <- data_mat[test_idx,  ]
  
  cat(sprintf("Processing Roll Origin: %d (Data Rows: %d to %d)...\n", 
              roll, min(train_idx), max(test_idx)))
  stanFit <- readRDS(stanfit_path)
  
  # ── A. Calculate In-Sample Baseline Scale Vector ───────────────────────────
  e_bot_pm  <- summary(stanFit, pars = "e_bot")$summary[, "mean"]
  e_A_pm    <- summary(stanFit, pars = "e_A")$summary[,  "mean"]
  e_B_pm    <- summary(stanFit, pars = "e_B")$summary[,  "mean"]
  e_top_pm  <- summary(stanFit, pars = "e_top")$summary[, "mean"]
  
  e_bot_mat <- matrix(e_bot_pm, nrow = TRAIN_SIZE, ncol = 4)
  E_full    <- cbind(e_top_pm, e_A_pm, e_B_pm, e_bot_mat)
  E_full    <- E_full[3:TRAIN_SIZE, ] # Drop ARMA initialization steps
  
  # Compute historical standard deviations across all 7 nodes
  scale_vec <- apply(E_full, 2, sd)
  scale_vec[scale_vec == 0] <- 1e-6 
  
  # ── B. Extract Base Forecast Matrix Elements ────────────────────────────────
  fc_bot   <- rstan::extract(stanFit, "fc_bot")$fc_bot     # [draws, H, 4]
  fc_A     <- rstan::extract(stanFit, "fc_A")$fc_A         # [draws, H]
  fc_B     <- rstan::extract(stanFit, "fc_B")$fc_B         # [draws, H]
  fc_total <- rstan::extract(stanFit, "fc_total")$fc_total # [draws, H]
  
  n_draws <- dim(fc_bot)[1]
  
  # ── C. Reconstruct MinT Projection Elements ─────────────────────────────────
  W_hat     <- cov(E_full)
  diag_S    <- diag(diag(W_hat))
  sum_sq    <- sum((W_hat - diag_S)^2)
  denom     <- sum(W_hat^2) - sum(diag(W_hat)^2)
  T_e       <- nrow(E_full)
  
  lambda_star <- if (denom < .Machine$double.eps) 1 else min(1, max(0, (sum_sq / denom) * ((T_e + 2) / T_e)))
  W_shrink    <- (1 - lambda_star) * W_hat + lambda_star * diag_S
  
  SP_mint   <- S_mat %*% solve(t(S_mat) %*% solve(W_shrink) %*% S_mat) %*% t(S_mat) %*% solve(W_shrink)
  
  # Reproducible sub-sampling setup
  set.seed(roll)
  sub_idx <- sample.int(n_draws, size = min(n_es_subsample, n_draws), replace = FALSE)
  
  # Data frames to collect results across horizons
  base_es_list <- vector("list", HORIZON)
  mint_es_list <- vector("list", HORIZON)
  
  # ── D. Loop Over Horizons ──────────────────────────────────────────────────
  for (h in seq_len(HORIZON)) {
    actual_bot   <- test[h, IDX$bot]
    actual_A     <- test[h, IDX$A]
    actual_B     <- test[h, IDX$B]
    actual_total <- test[h, IDX$total]
    actuals_joint <- c(actual_total, actual_A, actual_B, actual_bot) 
    
    # Assemble Joint Matrix: [Base Forecasts]
    joint_mat_base <- cbind(
      fc_total[, h],
      fc_A[, h],
      fc_B[, h],
      fc_bot[, h, 1],
      fc_bot[, h, 2],
      fc_bot[, h, 3],
      fc_bot[, h, 4]
    )
    
    # Assemble Joint Matrix: [MinT Reconciled Forecasts]
    joint_mat_mint <- joint_mat_base %*% t(SP_mint)
    
    # Calculate Scores
    base_es_list[[h]] <- data.frame(
      roll = roll, horizon = h, 
      energy_score = compute_comparable_energy_score(joint_mat_base, actuals_joint, scale_vec, sub_idx)
    )
    
    mint_es_list[[h]] <- data.frame(
      roll = roll, horizon = h, 
      energy_score = compute_comparable_energy_score(joint_mat_mint, actuals_joint, scale_vec, sub_idx)
    )
  }
  
  # ── E. Save Results and Free RAM Immediately ───────────────────────────────
  saveRDS(bind_rows(base_es_list), file.path(dir_output_base, sprintf("es_base_roll_%03d.rds", roll)))
  saveRDS(bind_rows(mint_es_list), file.path(dir_output_mint, sprintf("es_mint_roll_%03d.rds", roll)))
  
  rm(stanFit, fc_bot, fc_A, fc_B, fc_total, joint_mat_base, joint_mat_mint)
  gc()
}

cat("All 500 post-hoc energy scores computed successfully!\n")