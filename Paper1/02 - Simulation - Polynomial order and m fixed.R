# =============================================================================
# SCRIPT 02: SIMULATION - POLYNOMIAL ORDER (nu = 0, 2, 5) AND FIXED m
# =============================================================================
# Paper: Prass, T. S. and Pumi, G. (2021). On the behavior of the DFA and DCCA
#        in trend-stationary processes. Journal of Multivariate Analysis,
#        182, 104703.
#
# Runs two Monte Carlo simulation experiments (1,000 replications each):
#   1. Effect of the polynomial detrending order nu in {0, 2, 5} for large
#      window sizes m in {1000, 1010, ..., 1100} (n = 2,000) under independent
#      AR(1) and MA(1) processes (Figure 4 in the main manuscript).
#   2. Consistency as sample size n grows from 500 to 10,000 for a fixed window
#      size m = 100 under (a) two independent i.i.d. series and (b) independent
#      AR(1) and MA(1) series (Figure S2 in the Supplementary Material).
#
# Output: Saved in `results/part1_nu_mfixed/` (pre-computed results are already
#         included in the repository; run this script only to regenerate them).
# =============================================================================

library(DCCA)

outdir <- "results/part1_nu_mfixed"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Part 1: Effect of polynomial order nu in {0, 2, 5} for m in {1000, ..., 1100}
#         X(t) = 0.6 * X(t-1) + eps1(t)
#         Y(t) = eps2(t) + 0.6 * eps2(t-1),   eps1 independent of eps2
# =============================================================================
m_min <- 1000
m_max <- 1100
n <- 2000
n_rep <- 1000
ms <- seq(m_min, m_max, 10)
dm <- length(ms)

dfa12 <- dfa11 <- dfa10 <- numeric(dm)
dfa22 <- dfa21 <- dfa20 <- numeric(dm)
dcca2 <- dcca1 <- dcca0 <- numeric(dm)
dccacor2 <- dccacor1 <- dccacor0 <- numeric(dm)

for (r in 1:n_rep) {
  set.seed(100 + r)

  x <- numeric(n + 10)
  y <- numeric(n + 10)
  eps1 <- rnorm(n + 10)
  eps2 <- rnorm(n + 10)
  x[1] <- eps1[1]
  y[1] <- eps2[1]
  for (i in 2:(n + 10)) {
    x[i] <- 0.6 * x[i - 1] + eps1[i]
    y[i] <- eps2[i] + 0.6 * eps2[i - 1]
  }
  x <- x[-c(1:10)]
  y <- y[-c(1:10)]

  for (m_idx in seq_along(ms)) {
    temp <- DCCA::rhodcca(x, y, m = ms[m_idx], nu = 0)
    dfa10[m_idx]    <- temp$F2dfa1
    dfa20[m_idx]    <- temp$F2dfa2
    dcca0[m_idx]    <- temp$Fdcca
    dccacor0[m_idx] <- temp$rhodcca

    temp <- DCCA::rhodcca(x, y, m = ms[m_idx], nu = 2)
    dfa11[m_idx]    <- temp$F2dfa1
    dfa21[m_idx]    <- temp$F2dfa2
    dcca1[m_idx]    <- temp$Fdcca
    dccacor1[m_idx] <- temp$rhodcca

    temp <- DCCA::rhodcca(x, y, m = ms[m_idx], nu = 5)
    dfa12[m_idx]    <- temp$F2dfa1
    dfa22[m_idx]    <- temp$F2dfa2
    dcca2[m_idx]    <- temp$Fdcca
    dccacor2[m_idx] <- temp$rhodcca

    cat("Part 1 - rep:", r, "m:", ms[m_idx], "\n")
  }

  write.table(t(dfa10),    file = file.path(outdir, "dfa10_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dfa20),    file = file.path(outdir, "dfa20_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dcca0),    file = file.path(outdir, "dcca0_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dccacor0), file = file.path(outdir, "rhodcca0_1_2.txt"), row.names = FALSE, col.names = FALSE, append = TRUE)

  write.table(t(dfa11),    file = file.path(outdir, "dfa11_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dfa21),    file = file.path(outdir, "dfa21_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dcca1),    file = file.path(outdir, "dcca1_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dccacor1), file = file.path(outdir, "rhodcca1_1_2.txt"), row.names = FALSE, col.names = FALSE, append = TRUE)

  write.table(t(dfa12),    file = file.path(outdir, "dfa12_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dfa22),    file = file.path(outdir, "dfa22_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dcca2),    file = file.path(outdir, "dcca2_1_2.txt"),    row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(t(dccacor2), file = file.path(outdir, "rhodcca2_1_2.txt"), row.names = FALSE, col.names = FALSE, append = TRUE)
}

# =============================================================================
# Part 2: Fixed window size m = 100, varying sample size n in {500, ..., 10000}
# =============================================================================
n_min <- 500
n_max <- 10000
ns <- seq(n_min, n_max, 100)
dn <- length(ns)
m <- 100

# --- Model 1: Two independent i.i.d. N(0, 1) sequences ---
dccacor <- numeric(dn)
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- rnorm(ns[dn])
  y <- rnorm(ns[dn])
  for (n_idx in 1:dn) {
    dccacor[n_idx] <- DCCA::rhodcca(x[1:ns[n_idx]], y[1:ns[n_idx]], m = m)$rhodcca
    cat("mfixed1 - rep:", r, "n:", ns[n_idx], "\n")
  }
  write.table(t(dccacor), file = file.path(outdir, "mfixed1.txt"), row.names = FALSE, col.names = FALSE, append = TRUE)
}

# --- Model 2: Independent AR(1) and MA(1) processes ---
dccacor <- numeric(dn)
for (r in 1:n_rep) {
  set.seed(100 + r)
  x <- numeric(ns[dn] + 10)
  y <- numeric(ns[dn] + 10)
  eps1 <- rnorm(ns[dn] + 10)
  eps2 <- rnorm(ns[dn] + 10)
  x[1] <- eps1[1]
  y[1] <- eps2[1]
  for (i in 2:(ns[dn] + 10)) {
    x[i] <- 0.6 * x[i - 1] + eps1[i]
    y[i] <- eps2[i] + 0.6 * eps2[i - 1]
  }
  x <- x[-c(1:10)]
  y <- y[-c(1:10)]

  for (n_idx in 1:dn) {
    dccacor[n_idx] <- DCCA::rhodcca(x[1:ns[n_idx]], y[1:ns[n_idx]], m = m)$rhodcca
    cat("mfixed2 - rep:", r, "n:", ns[n_idx], "\n")
  }
  write.table(t(dccacor), file = file.path(outdir, "mfixed2.txt"), row.names = FALSE, col.names = FALSE, append = TRUE)
}
