# Parametric bootstrap experiments: simulation stage (scenarios SC2 and SC2_null).
#
# For each scenario, sample size n and Monte Carlo replication r:
#   1. simulate a bivariate ARFIMA(1, d, 1) sample;
#   2. estimate (d_k, S_k, S_12) with estimate(type = "joint");
#   3. run the bivariate parametric bootstrap (boot.dS12), which re-estimates
#      all five parameters in every bootstrap sample;
#   4. compute percentile intervals for d_k, S_k, S_12, F2dfa_k(m), Fdcca(m)
#      and rho(m).
# All `nrep` replications of one (scenario, n) block are saved together in a
# single file, so an interrupted run can be resumed by running the script
# again: a block whose file already exists is skipped, an unfinished one is
# recomputed from its first replication.
#
# Reference values:
#   parameters: d_k, S_k = Psi_kk * sigma_k^2, S_12 = Psi_12 * rho
#   curves:     asymptotic limit (rhoLim) and exact expectation (rhoE)
#
# Scenarios:
#   SC2       d = (0.45, -0.25), rho = 0.8   coverage and power
#   SC2_null  d = (0.45, 0),     rho = 0     size of the tests d_2 = 0, S_12 = 0
#
# Output (results/boot_SC2/, or results/TEST_boot_SC2/ when test_mode is TRUE):
#   <scenario>_n<n>.rds   list of `nrep` replications, one file per block
#   _config.rds           settings and reference values, read when the results are summarised
#   _timing_runs.csv      wall-clock time per block

library(DCCA)
source("00 - utils_parallel.R")

test_mode <- FALSE # TRUE: tiny run to check the script; FALSE: full study
run_null <- TRUE
B <- 200
alpha <- 0.05
overlap <- TRUE
seed0 <- 20260926
mtrunc <- 50000
outdir <- "results/boot_SC2"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCENARIO 2 (SC2)
#         X1 ~ ARFIMA(1,d1,1)
#         X2 ~ ARFIMA(1,d2,1)
#         (e1, e2) ~ N(0,Sigma), Sigma = [[1, rho], [rho, 1]]
# with
#     phi1(z) = phi2(z)
#     theta1(z) = theta2(z)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCENARIO 2 NULL (SC2_null)
#         same as SC2 with d2 = 0 and rho = 0, so that
#         d2 = 0 and S12 = 0 hold (used for the size of the tests)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Sample sizes and number of Monte Carlo replications: the same for every
# scenario, so they are set once instead of being repeated inside `scen`.
ns <- c(500, 1000, 2000, 5000)
nrep <- 200

scen <- list(
  SC2 = list(d = c(0.45, -0.25), ar = -0.4, ma = 0.2, rho = 0.8),
  SC2_null = list(d = c(0.45, 0), ar = -0.4, ma = 0.2, rho = 0)
)
if (!run_null) scen$SC2_null <- NULL

# Window sizes: same grid as in the first set of simulations, fixed for every n
ms <- seq(10, 300, by = 20)

if (test_mode) {
  ns <- 500
  nrep <- 3
  ms <- seq(10, 60, by = 10)
  B <- 20
  mtrunc <- 2000
  outdir <- "results/TEST_boot_SC2"
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Reference (population) values for one scenario, used later to compute bias,
# RMSE and coverage: `par` for the estimated parameters (d_k, S_k, S12), `lim`
# and `ex` for the F2dfa/Fdcca/rho curves against m.
reference_values <- function(s, ms) {
  # Asymptotic (m -> infinity) limits of F2dfa1(m), F2dfa2(m), Fdcca(m) and,
  # jointly with n -> infinity, of rho(m); exact analytical expressions, not
  # an approximation evaluated at a large but finite m.
  lim <- rhoLim(
    m = ms,
    d1 = s$d[1], ar1 = s$ar, ma1 = s$ma, var.innov1 = 1,
    d2 = s$d[2], ar2 = s$ar, ma2 = s$ma, var.innov2 = 1,
    cov.innov = s$rho
  )

  # S_k = Psi_kk (both series share ar, ma, so Psi_11 = Psi_22 = Psi_12 = psi)
  # and S12 = Psi_12 * rho, as in the "Data generating process" subsection.
  psi <- ((1 + s$ma) / (1 - s$ar))^2
  s12 <- psi * s$rho

  # Exact (finite m) expectations of F2dfa1(m), F2dfa2(m), Fdcca(m), rho(m),
  # from the autocovariance (G1, G2) and cross-covariance (G12) matrices of
  # the two ARFIMA processes; these converge to `lim` as m, n -> infinity.
  mm <- max(ms)
  G1 <- toeplitz(arfima.acf(lag.max = mm, d = s$d[1], ar = s$ar, ma = s$ma)$acf)
  G2 <- toeplitz(arfima.acf(lag.max = mm, d = s$d[2], ar = s$ar, ma = s$ma)$acf)
  ccf <- arfima.ccf(
    lag.max = mm, d1 = s$d[1], ar1 = s$ar, ma1 = s$ma,
    d2 = s$d[2], ar2 = s$ar, ma2 = s$ma, cov.innov = s$rho
  )$ccf
  # G12 is not symmetric when d1 != d2, so it is rebuilt from the cross-lag
  # vector `ccf` (lags -mm, ..., mm) into a full (mm+1) x (mm+1) matrix.
  idx <- outer(1:(mm + 1), 1:(mm + 1), function(i, j) j - i + (mm + 1))
  ex <- rhoE(m = ms, G1 = G1, G2 = G2, G12 = matrix(ccf[idx], nrow = mm + 1))

  list(
    par = c(d1 = s$d[1], d2 = s$d[2], S1 = psi, S2 = psi, S12 = s12),
    # lim$rho is a single number (the rho limit requires m, n -> infinity
    # jointly), recycled by data.frame() to match the per-m columns F2dfa1,
    # F2dfa2 and Fdcca, which are proper curves.
    lim = data.frame(m = ms, F2dfa1 = lim$F2dfa1, F2dfa2 = lim$F2dfa2,
                     Fdcca = lim$Fdcca, rho = lim$rhodcca),
    ex = data.frame(m = ms, F2dfa1 = ex$F2dfa1, F2dfa2 = ex$F2dfa2,
                    Fdcca = ex$Fdcca, rho = ex$rhodcca)
  )
}

# One Monte Carlo replication (one scenario, sample size and seed); returns
# the result, which the caller collects together with the other replications
# of the same (scenario, n) block before writing them to a single file.
run_one <- function(name, n, r, s, ms, B, alpha, overlap, trunc, seed0) {
  # Seed depends on scenario, n and r, so results do not change if replications
  # are re-run in a different order or split across a different number of cores.
  set.seed(seed0 + 100000L * match(name, c("SC2", "SC2_null")) + 100L * n + r)
  t0 <- proc.time()

  # Simulate the bivariate ARFIMA(1,d,1) sample with the package's own
  # bootstrap routine (Cholesky-correlated innovations), as in the first set of simulations.
  psi <- list(
    X1 = arfima.expansion(d = s$d[1], ar = s$ar, ma = s$ma, trunc = trunc),
    X2 = arfima.expansion(d = s$d[2], ar = s$ar, ma = s$ma, trunc = trunc)
  )
  Sigma <- matrix(c(1, s$rho, s$rho, 1), 2,
                  dimnames = list(c("X1", "X2"), c("X1", "X2")))
  x <- DCCA:::.boot.simulate(n, trunc, psi, Sigma)

  # Point estimates: `point` are the raw DFA/DCCA statistics (used later to
  # assess bias and RMSE of the curves), `fit` are the estimated parameters
  # (d_k, S_k, S12) that feed the bootstrap below.
  point <- rhodcca(x[, 1], x[, 2], m = ms, overlap = overlap)
  fit <- estimate(x, type = "joint", m = ms, nu = 0, overlap = overlap)
  est <- c(
    d1 = unname(fit$d[1]), d2 = unname(fit$d[2]),
    S1 = unname(fit$Seps[1]), S2 = unname(fit$Seps[2]),
    S12 = unname(fit$S12eps[1])
  )
  t1 <- proc.time()
  out <- list(
    name = name, n = n, r = r, ms = ms, est = est, ok = FALSE,
    point = list(F2dfa1 = point$F2dfa1, F2dfa2 = point$F2dfa2,
                 Fdcca = point$Fdcca, rho = point$rhodcca)
  )

  # The bootstrap simulates ARFIMA(0, d_k, 0) processes with innovation
  # covariance matrix built from (S1, S2, S12); it requires that matrix to be
  # positive semidefinite (Cauchy-Schwarz: |S12| <= sqrt(S1*S2)), which the
  # estimated (S1, S2, S12) need not satisfy. Replications that fail this
  # check are left with ok = FALSE and no confidence intervals.
  if (abs(est["S12"]) <= sqrt(est["S1"] * est["S2"])) {
    bc <- boot.control(
      m = ms, nu = 0, overlap = overlap, trunc = trunc, n_boot = B,
      parallel = FALSE, report = FALSE, seed = sample.int(1e6, 1)
    )
    # Re-estimates all five parameters (and, with curves = TRUE, the four
    # DFA/DCCA/rho curves) in each of the B bootstrap samples.
    b <- boot.dS12(
      n = n, d = fit$d, Seps = fit$Seps, S12eps = fit$S12eps,
      boot_control = bc, curves = TRUE
    )
    # Percentile confidence interval: the (alpha/2, 1 - alpha/2) quantiles of
    # the bootstrap distribution, one row/column per parameter or window size.
    ci <- function(M) t(apply(as.matrix(M), 2, quantile,
                              probs = c(alpha / 2, 1 - alpha / 2), names = FALSE))
    out$ci_par <- rbind(
      d1 = ci(b$d[, 1]), d2 = ci(b$d[, 2]),
      S1 = ci(b$Seps[, 1]), S2 = ci(b$Seps[, 2]), S12 = ci(b$S12eps)
    )
    colnames(out$ci_par) <- c("lo", "hi")
    out$ci_curve <- lapply(
      list(F2dfa1 = b$F2dfa1, F2dfa2 = b$F2dfa2, Fdcca = b$Fdcca, rho = b$rho), ci
    )
    out$ok <- TRUE
  }
  # Time split between (data generation + point/joint estimation) and the
  # bootstrap itself, reported later in the computational cost table.
  out$time_gen_est <- unname((t1 - t0)["elapsed"])
  out$time_boot <- unname((proc.time() - t1)["elapsed"])
  out
}

# One block per (scenario, n); `file` collects all `nrep` replications of
# that block into a single .rds, and its existence is what makes the script
# resumable (a finished block is skipped; an unfinished one restarts from its
# first replication, not from where it stopped).
blocks <- do.call(rbind, lapply(names(scen), function(nm) data.frame(name = nm, n = ns)))
blocks$file <- file.path(outdir, sprintf("%s_n%d.rds", blocks$name, blocks$n))
blocks <- blocks[order(-blocks$n), ] # largest n (most expensive blocks) first
todo <- blocks[!file.exists(blocks$file), ]
cat(sprintf("%d blocks in total, %d pending (%d replications each, B = %d)\n",
            nrow(blocks), nrow(todo), nrep, B))

# Population values for every (scenario, n) combination, computed once (they
# do not depend on the replication r) and read later when the results are summarised.
reference <- list()
for (nm in names(scen)) {
  for (n in ns) reference[[paste0(nm, "_n", n)]] <- reference_values(scen[[nm]], ms)
}
saveRDS(
  list(reference = reference, scen = scen, ns = ns, nrep = nrep, B = B,
       alpha = alpha, mtrunc = mtrunc),
  file.path(outdir, "_config.rds")
)

n_cores <- if (test_mode) 3 else max(1, parallel::detectCores() - 1)
f_time <- file.path(outdir, "_timing_runs.csv")
for (i in seq_len(nrow(todo))) {
  bl <- todo[i, ]
  cat(sprintf("\n%s, n = %d (block %d/%d)\n", bl$name, bl$n, i, nrow(todo)))
  t_start <- proc.time()
  # The nrep replications of this block are spread across the cores; the
  # whole block is written to disk only once all of them have finished.
  res <- parallel_map(
    seq_len(nrep),
    function(r) run_one(bl$name, bl$n, r, scen[[bl$name]], ms, B, alpha, overlap, mtrunc, seed0),
    n_cores,
    export = c("run_one", "bl", "scen", "ms", "B", "alpha", "overlap", "mtrunc", "seed0")
  )
  saveRDS(res, bl$file)
  elapsed <- (proc.time() - t_start)["elapsed"]
  errors <- Filter(function(z) inherits(z, "error"), res) # replications that raised an error
  cat(sprintf("Elapsed: %.1f min (%d cores), %d errors\n", elapsed / 60, n_cores, length(errors)))
  write.table(
    data.frame(
      when = format(Sys.time()), scenario = bl$name, n = bl$n, nrep = nrep, B = B, n_cores = n_cores,
      elapsed_min = elapsed / 60, R = R.version.string, DCCA = as.character(packageVersion("DCCA"))
    ),
    f_time,
    sep = ";", row.names = FALSE, col.names = !file.exists(f_time), append = file.exists(f_time)
  )
  if (length(errors) > 0) print(errors[[1]]) # first error, for a quick look
}
