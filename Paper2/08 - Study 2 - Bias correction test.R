# Study II: test whether bootstrap bias correction improves the point
# estimates, using the REAL parametric bootstrap already run by
# sim2_bootstrap_SC2.R (results/boot_SC2/SC2_n<n>.rds) -- not an
# approximation, the actual boot.dS12() output for scenario SC2.
#
# Bootstrap bias correction: theta_BC = 2*theta_hat - mean(theta_star), where
# theta_star are the bootstrap replicates. sim2_bootstrap_SC2.R only saved
# the 95% percentile interval (lo, hi) for each replication, not the B raw
# bootstrap draws, so mean(theta_star) is approximated here as (lo+hi)/2 --
# exact for a symmetric bootstrap distribution, an approximation otherwise.
# This still uses only data already produced; no new simulation is run.
#
# For each (scenario, n) block and each parameter (d1, d2, S1, S2, S12),
# this compares, across the nrep replications with ok = TRUE:
#   - bias/RMSE of the original point estimate theta_hat
#   - bias/RMSE of the corrected estimate theta_BC
#
# Input : results/boot_SC2/_config.rds, results/boot_SC2/SC2_n<n>.rds
# Output: results/study2_bias_correction.csv (per-parameter, per-n summary)
#         printed summary table

cfg <- readRDS("results/boot_SC2/_config.rds")
ns  <- cfg$ns
params <- c("d1", "d2", "S1", "S2", "S12")

summarise_block <- function(name, n) {
  res <- readRDS(sprintf("results/boot_SC2/%s_n%d.rds", name, n))
  truth <- cfg$reference[[paste0(name, "_n", n)]]$par

  ok <- vapply(res, function(z) isTRUE(z$ok), logical(1))
  res_ok <- res[ok]

  # ci_par's rownames are not reliably set to the parameter names (only
  # S12's row is labelled, as "S12eps"; the others are blank strings), but
  # its rows are always in the same order as names(z$est) = d1,d2,S1,S2,S12.
  # Index ci_par by that position instead of by name.
  out <- lapply(params, function(p) {
    pos        <- match(p, names(res_ok[[1]]$est))
    theta_hat  <- vapply(res_ok, function(z) z$est[p], numeric(1))
    lo         <- vapply(res_ok, function(z) z$ci_par[pos, "lo"], numeric(1))
    hi         <- vapply(res_ok, function(z) z$ci_par[pos, "hi"], numeric(1))
    theta_star_mean <- (lo + hi) / 2
    theta_BC   <- 2 * theta_hat - theta_star_mean

    true_val <- unname(truth[p])
    data.frame(
      scenario = name, n = n, parameter = p, n_ok = length(theta_hat),
      true_val = true_val,
      bias_before = mean(theta_hat - true_val),
      bias_after  = mean(theta_BC - true_val),
      rmse_before = sqrt(mean((theta_hat - true_val)^2)),
      rmse_after  = sqrt(mean((theta_BC - true_val)^2))
    )
  })
  do.call(rbind, out)
}

all_summary <- do.call(rbind, lapply(ns, function(n) summarise_block("SC2", n)))
write.csv(all_summary, "results/study2_bias_correction.csv", row.names = FALSE)

cat("Bootstrap bias correction test (SC2, using the real boot.dS12 output,\n")
cat("theta_star mean approximated as (lo+hi)/2 from the saved 95% CI):\n\n")
print(all_summary, row.names = FALSE, digits = 4)
