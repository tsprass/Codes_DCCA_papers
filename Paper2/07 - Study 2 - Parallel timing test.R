# Study II: single replication (one scenario, one n, one seed) of the
# estimation + parametric bootstrap pipeline used by sim2_bootstrap_SC2.R,
# but with the bootstrap step itself run in parallel (boot.control(parallel
# = TRUE)) instead of serially. sim2_bootstrap_SC2.R uses parallel = FALSE
# for the bootstrap because it already parallelizes the OUTER loop over the
# nrep Monte Carlo replications (nesting two parallel backends would not
# help and complicates cluster management); this script is for running, or
# timing, a single replication with the INNER B bootstrap draws themselves
# parallelized, e.g. to decide whether that is worthwhile when nrep is small
# or when only a handful of replications are needed.
#
# boot.control(parallel = TRUE) is already implemented by the package
# (PSOCK cluster via doSNOW + doRNG, see boot.control()'s "parallel"/
# "n_cores" arguments); this script does not add any new parallel
# infrastructure, it only exercises that existing option and times it.
#
# Edit the CONFIGURATION block below to select which replication to run.
#
# Output: printed point estimates, bootstrap percentile CIs, and the
#         elapsed time of the bootstrap step (parallel, and serial for
#         comparison when compare_serial = TRUE).

library(DCCA)

# --- CONFIGURATION ----------------------------------------------------------
scenario   <- "SC2"  # "SC2" (d = (0.45, -0.25), rho = 0.8) or "SC2_null"
n          <- 2000   # sample size
r          <- 1L     # Monte Carlo replication index (only affects the seed)
seed0      <- 20260926
mtrunc     <- 50000
overlap    <- TRUE
ms         <- seq(10, 300, by = 20)
alpha      <- 0.05
B          <- 1000  # number of bootstrap replications
n_cores    <- max(1, parallel::detectCores() - 1)
compare_serial <- FALSE # also time the same bootstrap with parallel = FALSE
# -----------------------------------------------------------------------------

scen <- list(
  SC2      = list(d = c(0.45, -0.25), ar = -0.4, ma = 0.2, rho = 0.8),
  SC2_null = list(d = c(0.45, 0),     ar = -0.4, ma = 0.2, rho = 0)
)[[scenario]]
if (is.null(scen)) stop("Unknown scenario: ", scenario)

# Same seeding rule as sim2_bootstrap_SC2.R, so this replication reproduces
# exactly the one that script would generate for the same (scenario, n, r).
set.seed(seed0 + 100000L * match(scenario, c("SC2", "SC2_null")) + 100L * n + r)

# ---------------------------------------------------------------------------
# Simulate the bivariate ARFIMA(1,d,1) sample and compute the point estimates
# that feed the bootstrap, exactly as in sim2_bootstrap_SC2.R's run_one().
# ---------------------------------------------------------------------------
psi <- list(
  X1 = arfima.expansion(d = scen$d[1], ar = scen$ar, ma = scen$ma, trunc = mtrunc),
  X2 = arfima.expansion(d = scen$d[2], ar = scen$ar, ma = scen$ma, trunc = mtrunc)
)
Sigma <- matrix(c(1, scen$rho, scen$rho, 1), 2,
                dimnames = list(c("X1", "X2"), c("X1", "X2")))
x <- DCCA:::.boot.simulate(n, mtrunc, psi, Sigma)

fit <- estimate(x, type = "joint", m = ms, nu = 0, overlap = overlap)
cat("Point estimates:\n")
print(fit$d)
print(fit$Seps)
print(fit$S12eps)

if (abs(fit$S12eps) > sqrt(fit$Seps[1] * fit$Seps[2])) {
  stop("Cauchy-Schwarz violated for this replication (S12 not admissible); ",
       "pick a different r.")
}

# ---------------------------------------------------------------------------
# Bootstrap with the inner B replications run in parallel.
# ---------------------------------------------------------------------------
bc_parallel <- boot.control(
  m = ms, nu = 0, overlap = overlap, trunc = mtrunc, n_boot = B,
  parallel = TRUE, n_cores = n_cores, report = TRUE, seed = sample.int(1e6, 1)
)

cat(sprintf("\nRunning bootstrap with parallel = TRUE (%d cores, B = %d)...\n",
            n_cores, B))
t_par <- system.time(
  boot_parallel <- boot.dS12(
    n = n, d = fit$d, Seps = fit$Seps, S12eps = fit$S12eps,
    boot_control = bc_parallel, curves = TRUE
  )
)
cat(sprintf("Elapsed (parallel): %.1f s\n", t_par["elapsed"]))

ci <- function(M) t(apply(as.matrix(M), 2, quantile,
                          probs = c(alpha / 2, 1 - alpha / 2), names = FALSE))
ci_par <- rbind(
  d1 = ci(boot_parallel$d[, 1]), d2 = ci(boot_parallel$d[, 2]),
  S1 = ci(boot_parallel$Seps[, 1]), S2 = ci(boot_parallel$Seps[, 2]),
  S12 = ci(boot_parallel$S12eps)
)
colnames(ci_par) <- c("lo", "hi")
cat("\n95% percentile CI (parallel bootstrap):\n")
print(ci_par)

# ---------------------------------------------------------------------------
# Optional: same bootstrap run serially, for a direct timing comparison on
# this machine before committing to parallel = TRUE for the full study.
# ---------------------------------------------------------------------------
if (compare_serial) {
  bc_serial <- boot.control(
    m = ms, nu = 0, overlap = overlap, trunc = mtrunc, n_boot = B,
    parallel = FALSE, report = TRUE, seed = bc_parallel$seed
  )
  cat(sprintf("\nRunning bootstrap with parallel = FALSE (B = %d) for comparison...\n", B))
  t_ser <- system.time(
    boot_serial <- boot.dS12(
      n = n, d = fit$d, Seps = fit$Seps, S12eps = fit$S12eps,
      boot_control = bc_serial, curves = TRUE
    )
  )
  cat(sprintf("Elapsed (serial):   %.1f s\n", t_ser["elapsed"]))
  cat(sprintf(
    "\nSpeedup (serial / parallel): %.2fx (%d cores)\n",
    t_ser["elapsed"] / t_par["elapsed"], n_cores
  ))
}
