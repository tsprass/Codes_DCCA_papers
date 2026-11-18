# =============================================================================
# SCRIPT 04: SIMULATION - SHORT-MEMORY SCENARIOS 1 TO 4 (SUPPLEMENTARY MATERIAL)
# =============================================================================
# Paper: Prass, T. S. and Pumi, G. (2021). On the behavior of the DFA and DCCA
#        in trend-stationary processes. Journal of Multivariate Analysis,
#        182, 104703.
#
# Simulates 1,000 Monte Carlo replications (n = 2,000, m in {3, 5, ..., 101},
# nu = 0) across all four short-memory scenarios studied in Section 1 of the
# Supplementary Material:
#   - Scenario 1a (`_1.txt`):      Two independent i.i.d. N(0, 1) sequences
#   - Scenario 1b (`_1_2.txt`):    Independent AR(1) and MA(1) processes
#   - Scenario 2a (`_2.txt`):      Signal plus noise with sigma_eps^2 = 4
#   - Scenario 2b (`_2_2.txt`):    Signal plus noise with sigma_eps^2 = 64
#   - Scenario 2c (`_2_cor5.txt`): Bivariate white noise with rho = 0.5
#   - Scenario 2d (`_2_cor8.txt`): Bivariate white noise with rho = 0.8
#   - Scenario 3a (`_3.txt`):      Cross-correlation driven by MA(20) structure
#   - Scenario 3b (`_3_2.txt`):    Cross-correlation driven by AR(1) structure
#   - Scenario 4a (`_4.txt`):      Common AR(1) innovation with AR(1) filters
#   - Scenario 4b (`_4_2.txt`):    Common AR(1) innovation with MA(1) filters
#
# Output: Saved in `results/part2_scenarios/` (pre-computed results are already
#         included in the repository; run this script only to regenerate them).
# =============================================================================

library(DCCA)
library(MASS)

outdir <- "results/part2_scenarios"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m_min <- 3
m_max <- 101
ms <- seq(m_min, m_max, 2)
dm <- length(ms)
n <- 2000
n_rep <- 1000

save_rep <- function(suffix, dfa1, dfa2, dcca12, dccacor) {
  write.table(t(dfa1),    file = file.path(outdir, paste0("dfa1", suffix)),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dfa2),    file = file.path(outdir, paste0("dfa2", suffix)),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dcca12),  file = file.path(outdir, paste0("dcca", suffix)),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dccacor), file = file.path(outdir, paste0("rhodcca", suffix)), row.names = FALSE, col.names = FALSE, append = TRUE)
}

run_curve <- function(x, y) {
  dfa1 <- dfa2 <- dcca12 <- dccacor <- numeric(dm)
  for (m_idx in seq_along(ms)) {
    temp <- DCCA::rhodcca(x, y, m = ms[m_idx], nu = 0)
    dfa1[m_idx]    <- temp$F2dfa1
    dfa2[m_idx]    <- temp$F2dfa2
    dcca12[m_idx]  <- temp$Fdcca
    dccacor[m_idx] <- temp$rhodcca
  }
  list(dfa1 = dfa1, dfa2 = dfa2, dcca12 = dcca12, dccacor = dccacor)
}

# =============================================================================
# Scenario 1a: Two independent i.i.d. N(0, 1) sequences
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(n)
  y <- rnorm(n)
  res <- run_curve(x, y)
  save_rep("_1.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 1a - rep:", r, "\n")
}

# =============================================================================
# Scenario 1b: Independent AR(1) and MA(1) processes
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- y <- numeric(n + 10)
  eps1 <- rnorm(n + 10)
  eps2 <- rnorm(n + 10)
  x[1] <- eps1[1]
  y[1] <- eps2[1]
  for (i in 2:(n + 10)) {
    x[i] <- 0.6 * x[i - 1] + eps1[i]
    y[i] <- eps2[i] + 0.6 * eps2[i - 1]
  }
  res <- run_curve(x[-c(1:10)], y[-c(1:10)])
  save_rep("_1_2.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 1b - rep:", r, "\n")
}

# =============================================================================
# Scenario 2a: Signal plus noise (sigma_eps = 2, i.e., variance = 4)
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(n)
  y <- 2 * x + 3 + rnorm(n, sd = 2)
  res <- run_curve(x, y)
  save_rep("_2.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 2a - rep:", r, "\n")
}

# =============================================================================
# Scenario 2b: Signal plus noise (sigma_eps = 8, i.e., variance = 64)
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(n)
  y <- 2 * x + 3 + rnorm(n, sd = 8)
  res <- run_curve(x, y)
  save_rep("_2_2.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 2b - rep:", r, "\n")
}

# =============================================================================
# Scenario 2c: Bivariate Gaussian white noise with correlation rho = 0.5
# =============================================================================
SIG5 <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
for (r in 1:n_rep) {
  set.seed(100 + r)
  vec5 <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = SIG5)
  res <- run_curve(vec5[, 1], vec5[, 2])
  save_rep("_2_cor5.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 2c - rep:", r, "\n")
}

# =============================================================================
# Scenario 2d: Bivariate Gaussian white noise with correlation rho = 0.8
# =============================================================================
SIG8 <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
for (r in 1:n_rep) {
  set.seed(100 + r)
  vec8 <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = SIG8)
  res <- run_curve(vec8[, 1], vec8[, 2])
  save_rep("_2_cor8.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 2d - rep:", r, "\n")
}

# =============================================================================
# Scenario 3a: Cross-correlation induced by MA(20) filter
# =============================================================================
thetas <- (20:1) / 10
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(n + length(thetas))
  y <- numeric(n + length(thetas))
  for (i in (length(thetas) + 1):(n + length(thetas))) {
    y[i] <- x[i] + sum(thetas * x[(i - 1):(i - length(thetas))])
  }
  res <- run_curve(x[-(1:length(thetas))], y[-(1:length(thetas))])
  save_rep("_3.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 3a - rep:", r, "\n")
}

# =============================================================================
# Scenario 3b: Cross-correlation induced by AR(1) filter
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(n + 10)
  y <- numeric(n + 10)
  y[1] <- x[1]
  for (i in 2:(n + 10)) {
    y[i] <- 0.6 * y[i - 1] + x[i]
  }
  res <- run_curve(x[-(1:10)], y[-(1:10)])
  save_rep("_3_2.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 3b - rep:", r, "\n")
}

# =============================================================================
# Scenario 4a: Common AR(1) innovation driving two AR(1) filters (AR(2) marginals)
# =============================================================================
burnEps <- 50
burnXY  <- 10
for (r in 1:n_rep) {
  set.seed(100 + r)
  nt <- rnorm(n + burnEps + burnXY)
  et <- numeric(n + burnEps + burnXY)
  et[1] <- nt[1]
  for (i in 2:(n + burnEps + burnXY)) et[i] <- 0.7 * et[i - 1] + nt[i]
  et <- et[-c(1:burnEps)]

  X <- Y <- numeric(n + burnXY)
  X[1] <- et[1]
  Y[1] <- et[1]
  for (i in 2:(n + burnXY)) {
    X[i] <- 0.2 * X[i - 1] + et[i]
    Y[i] <- 0.6 * Y[i - 1] + et[i]
  }
  res <- run_curve(X[-c(1:burnXY)], Y[-c(1:burnXY)])
  save_rep("_4.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 4a - rep:", r, "\n")
}

# =============================================================================
# Scenario 4b: Common AR(1) innovation driving two MA(1) filters (ARMA(1,1) marginals)
# =============================================================================
for (r in 1:n_rep) {
  set.seed(100 + r)
  nt <- rnorm(n + burnEps + burnXY)
  et <- numeric(n + burnEps + burnXY)
  et[1] <- nt[1]
  for (i in 2:(n + burnEps + burnXY)) et[i] <- 0.7 * et[i - 1] + nt[i]
  et <- et[-c(1:burnEps)]

  X <- Y <- numeric(n + burnXY)
  X[1] <- et[1]
  Y[1] <- et[1]
  for (i in 2:(n + burnXY)) {
    X[i] <- et[i] + 0.2 * et[i - 1]
    Y[i] <- et[i] + 0.6 * et[i - 1]
  }
  res <- run_curve(X[-c(1:burnXY)], Y[-c(1:burnXY)])
  save_rep("_4_2.txt", res$dfa1, res$dfa2, res$dcca12, res$dccacor)
  cat("Scenario 4b - rep:", r, "\n")
}
