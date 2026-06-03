# =============================================================================
# Replication: Wickramasuriya, Athanasopoulos & Hyndman (2019)
# Section 3.3: Smoothing hierarchy DGP
# Generates a single dataset: 50 training + 500 rolling observations
#
# ARMA(2,2) coefficients drawn from a 3-level Gaussian hyperprior hierarchy:
#   Level 0: grand mean sampled uniformly from Table 1 ranges
#   Level 1: group means mu_A, mu_B ~ N(grand_mean, sigma2_hyper1)
#   Level 2: series coefficients    ~ N(mu_parent,  sigma2_hyper2)
#
# Parameters are sampled in the order that respects their dependence:
#   phi2, theta2 first (free ranges), then phi1, theta1 conditional on them.
# =============================================================================

library(MASS)

set.seed(2026)

# =============================================================================
# 1.  HIERARCHY STRUCTURE
# =============================================================================
#
#          Total (y)
#         /         \
#        A           B
#       / \         / \
#     AA  AB       BA  BB
#
#  m = 7 series, n = 4 bottom-level series
#  Ordering: [Total, A, B, AA, AB, BA, BB]

S <- rbind(
  c(1, 1, 1, 1),   # Total
  c(1, 1, 0, 0),   # A
  c(0, 0, 1, 1),   # B
  c(1, 0, 0, 0),   # AA
  c(0, 1, 0, 0),   # AB
  c(0, 0, 1, 0),   # BA
  c(0, 0, 0, 1)    # BB
)

m <- nrow(S)  # 7
n <- ncol(S)  # 4

# =============================================================================
# 2.  INNOVATION COVARIANCE (defined here so generate_smooth can receive it
#     explicitly rather than resolving it from the global environment)
# =============================================================================

Sigma_bottom <- matrix(
  data = c(5, 3, 2, 1,
           3, 4, 2, 1,
           2, 2, 5, 3,
           1, 1, 3, 4),
  nrow = 4, ncol = 4, byrow = TRUE
)

# =============================================================================
# 3.  3-LEVEL HYPERPRIOR SAMPLING FOR ARMA(2,2) COEFFICIENTS
# =============================================================================

sigma2_hyper1 <- 0.02   # grand mean  -> group means
sigma2_hyper2 <- 0.01   # group means -> series coefficients

# ---------------------------------------------------------------------------
# Stationarity (AR) and invertibility (MA) checks for ARMA(2,2)
# AR polynomial roots outside unit circle iff:
#   |phi2| < 1,  phi1 + phi2 < 1,  phi2 - phi1 < 1
# Same conditions apply to theta1, theta2 for invertibility.
# ---------------------------------------------------------------------------
is_valid_ar2 <- function(phi1, phi2) {
  abs(phi2) < 1 && (phi1 + phi2) < 1 && (phi2 - phi1) < 1
}

is_valid_ma2 <- function(theta1, theta2) {
  abs(theta2) < 1 && (theta1 + theta2) < 1 && (theta2 - theta1) < 1
}

# --- Level 0: sample grand means from Table 1 ranges ----------------------
sample_grand_means <- function() {
  phi2_gm   <- runif(1, 0.5, 0.7)
  phi1_gm   <- runif(1, phi2_gm - 0.9, 0.9 - phi2_gm)
  theta2_gm <- runif(1, 0.5, 0.7)
  theta1_gm <- runif(1, -(0.9 + theta2_gm) / 3.2, (0.9 + theta2_gm) / 3.2)
  c(phi1 = phi1_gm, phi2 = phi2_gm, theta1 = theta1_gm, theta2 = theta2_gm)
}

# --- Levels 1 & 2: cascade Gaussian hyperpriors ---------------------------
# Resamples until all four series satisfy stationarity and invertibility.
sample_coefficients <- function(grand_means, max_attempts = 10000) {
  
  series_names <- c("AA", "AB", "BA", "BB")
  
  for (attempt in seq_len(max_attempts)) {
    
    # Level 1: group means (2 groups x 4 parameters)
    group_means <- sapply(grand_means, function(gm) {
      rnorm(2, mean = gm, sd = sqrt(sigma2_hyper1))
    })
    # group_means: rows = [A, B], cols = [phi1, phi2, theta1, theta2]
    
    # Level 2: series coefficients (4 series x 4 parameters)
    series_coefs <- rbind(
      rnorm(4, mean = group_means[1, ], sd = sqrt(sigma2_hyper2)),  # AA (parent A)
      rnorm(4, mean = group_means[1, ], sd = sqrt(sigma2_hyper2)),  # AB (parent A)
      rnorm(4, mean = group_means[2, ], sd = sqrt(sigma2_hyper2)),  # BA (parent B)
      rnorm(4, mean = group_means[2, ], sd = sqrt(sigma2_hyper2))   # BB (parent B)
    )
    colnames(series_coefs) <- c("phi1", "phi2", "theta1", "theta2")
    rownames(series_coefs) <- series_names
    
    # Validate every series
    valid <- all(sapply(series_names, function(nm) {
      is_valid_ar2(series_coefs[nm, "phi1"], series_coefs[nm, "phi2"]) &&
        is_valid_ma2(series_coefs[nm, "theta1"], series_coefs[nm, "theta2"])
    }))
    
    if (valid) {
      return(list(
        grand_means  = grand_means,
        group_means  = group_means,
        series_coefs = series_coefs,
        attempts     = attempt
      ))
    }
  }
  
  stop(sprintf(
    "Could not draw valid ARMA coefficients after %d attempts. ",
    "Consider widening sigma2_hyper1/sigma2_hyper2 or tightening grand mean ranges.",
    max_attempts
  ))
}

# =============================================================================
# 4.  SECTION 3.3 DGP: smoothing hierarchy
# =============================================================================

arma22_filter <- function(innov, phi1, phi2, theta1, theta2, n_burnin = 200) {
  n_total <- length(innov) - n_burnin
  y <- numeric(length(innov))
  e <- innov
  for (t in 3:length(innov)) {
    y[t] <- phi1   * y[t-1] + phi2   * y[t-2] +
      e[t]   +
      theta1 * e[t-1] + theta2 * e[t-2]
  }
  tail(y, n_total)
}

generate_smooth <- function(n_total, series_coefs, Sigma_bottom,
                            sigma2_1 = 10, sigma2_2 = 6) {
  
  sigma1   <- sqrt(sigma2_1)
  sigma2   <- sqrt(sigma2_2)
  n_burnin <- 200
  n_sim    <- n_total + n_burnin
  
  # Joint ARMA innovations from bottom-level covariance
  joint_innov <- mvrnorm(n_sim, mu = rep(0, 4), Sigma = Sigma_bottom)
  
  series_names <- c("AA", "AB", "BA", "BB")
  z_mat <- matrix(NA, n_total, 4)
  for (j in seq_along(series_names)) {
    nm <- series_names[j]
    z_mat[, j] <- arma22_filter(
      innov  = joint_innov[, j],
      phi1   = series_coefs[nm, "phi1"],
      phi2   = series_coefs[nm, "phi2"],
      theta1 = series_coefs[nm, "theta1"],
      theta2 = series_coefs[nm, "theta2"],
      n_burnin = n_burnin
    )
  }
  
  # Wickramasuriya Section 3.3 smoothing noise:
  #   nu  cancels within groups:  AA-AB -> A,  BA-BB -> B
  #   om  cancels across groups:  A-B   -> Total
  # Net effect: noise-free at Total; nu-free at A and B; raw at bottom level.
  nu <- rnorm(n_total, 0, sigma1)
  om <- rnorm(n_total, 0, sigma2)
  
  bottom_mat <- cbind(
    AA = z_mat[, 1] - nu - 0.5 * om,
    AB = z_mat[, 2] + nu - 0.5 * om,
    BA = z_mat[, 3] - nu + 0.5 * om,
    BB = z_mat[, 4] + nu + 0.5 * om
  )
  
  full_mat <- t(S %*% t(bottom_mat))
  colnames(full_mat) <- c("Total", "A", "B", "AA", "AB", "BA", "BB")
  full_mat
}

# =============================================================================
# 5.  GENERATE DATA
# =============================================================================

TRAIN_SIZE <- 50
N_ROLL     <- 500
HORIZON    <- 5
N_TOTAL    <- TRAIN_SIZE + N_ROLL + HORIZON   # 555

grand_means  <- sample_grand_means()
coef_draw    <- sample_coefficients(grand_means)
series_coefs <- coef_draw$series_coefs
group_means  <- coef_draw$group_means

message(sprintf("Valid coefficients found after %d attempt(s).", coef_draw$attempts))

synthetic_data <- generate_smooth(
  n_total      = N_TOTAL,
  series_coefs = series_coefs,
  Sigma_bottom = Sigma_bottom
)

# =============================================================================
# 6.  SAVE
# =============================================================================

write.csv(synthetic_data, file = "synthetic_hierarchy.csv", row.names = FALSE)
save(series_coefs, group_means, grand_means, Sigma_bottom,
     file = "true_params.Rdata")