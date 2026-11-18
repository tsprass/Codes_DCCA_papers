# Parametric bootstrap experiments: all tables.
#
# This script only reads the replications saved by the simulation stage
# (results/boot_SC2/*.rds) -- no simulation is run here. It has three parts:
#
#   Part 1 (summaries):  bias, RMSE, coverage and width of the parameter
#     estimates, coverage of the F2dfa/Fdcca/rho curves, size and power of
#     the tests d_k = 0 and S_12 = 0, and computing times.
#   Part 2 (tables):     the summaries written to tables/.
#   Part 3 (appendix):   reduced coverage table for the five parameters.
#
# Output (tables/, or tables/TEST_* when test_mode is TRUE):
#   table_boot_params.{csv,tex}            parameters: bias, RMSE, coverage, width
#   table_boot_tests.{csv,tex}             size and power of the zero restrictions
#   table_boot_time.{csv,tex}              computing times
#   table_boot_curves.csv                  curve-level coverage, bias and RMSE
#   table_boot_rho.tex, table_boot_Fdcca.tex   curve coverage by n and window size (SC2)
#   table_coverage_params_SC2.{csv,tex}    reduced parameter coverage table
# tables/table_boot_params.csv and tables/table_boot_curves.csv are also the
# input of the figures.

suppressPackageStartupMessages(library(xtable))

test_mode <- FALSE # must match test_mode of the simulation stage
dir_in <- if (test_mode) "results/TEST_boot_SC2" else "results/boot_SC2"
prefix <- if (test_mode) "tables/TEST_" else "tables/"
dir.create("tables", showWarnings = FALSE)

# =============================================================================
# Part 1: summaries by (scenario, n) block
# =============================================================================

cfg <- readRDS(file.path(dir_in, "_config.rds"))
reference <- cfg$reference

# Each file is one (scenario, n) block: a list of `nrep` replications, saved
# together by the simulation stage.
files <- list.files(dir_in, pattern = "^SC2.*\\.rds$", full.names = TRUE)
groups <- lapply(files, readRDS)
names(groups) <- vapply(groups, function(g) paste0(g[[1]]$name, "_n", g[[1]]$n), "")
order_groups <- order(vapply(groups, function(g) g[[1]]$name, ""), vapply(groups, function(g) g[[1]]$n, 0))
groups <- groups[order_groups]

par_names <- c("d1", "d2", "S1", "S2", "S12")
curve_names <- c("F2dfa1", "F2dfa2", "Fdcca", "rho")
par_rows <- list()
curve_rows <- list()
test_rows <- list()
time_rows <- list()

for (key in names(groups)) {
  g_all <- groups[[key]]
  # Replications that raised an error inside run_one (e.g. a numerical
  # failure in estimate() or boot.dS12()) are stored as raw error/condition
  # objects by the simulation stage and lack every other field; they must be
  # removed before g is used, not just filtered by $ok (which is NULL, not
  # FALSE, on these objects and would misalign Filter()).
  crashed <- Filter(function(z) inherits(z, "error"), g_all)
  g <- Filter(function(z) !inherits(z, "error"), g_all)
  n_crash <- length(crashed)
  
  ref <- reference[[key]]
  name <- g[[1]]$name
  n <- g[[1]]$n
  ms <- g[[1]]$ms
  ok <- Filter(function(z) isTRUE(z$ok), g)
  n_rep <- length(g_all)
  n_ok <- length(ok)

  # Parameters: bias, RMSE, coverage, width
  est <- do.call(rbind, lapply(ok, function(z) z$est))
  for (j in seq_along(par_names)) {
    true_value <- ref$par[j]
    lo <- vapply(ok, function(z) z$ci_par[j, "lo"], 0)
    hi <- vapply(ok, function(z) z$ci_par[j, "hi"], 0)
    par_rows[[length(par_rows) + 1]] <- data.frame(
      scenario = name, n = n, param = par_names[j], true = unname(true_value),
      mean_est = mean(est[, j]), bias = mean(est[, j]) - true_value,
      rmse = sqrt(mean((est[, j] - true_value)^2)),
      coverage = mean(lo <= true_value & true_value <= hi), width = mean(hi - lo),
      n_ok = n_ok, n_rep = n_rep
    )
  }

  # Tests: H0 is rejected when zero lies outside the interval
  reject <- function(j) {
    mean(vapply(ok, function(z) !(z$ci_par[j, "lo"] <= 0 && 0 <= z$ci_par[j, "hi"]), TRUE))
  }
  test_rows[[length(test_rows) + 1]] <- data.frame(
    scenario = name, n = n, rej_d1 = reject(1), rej_d2 = reject(2), rej_S12 = reject(5),
    n_ok = n_ok, n_rep = n_rep, n_crash = n_crash, frac_not_psd = 1 - n_ok / n_rep
  )

  # Curves: coverage of the exact expectation and of the asymptotic limit
  for (q in curve_names) {
    point <- do.call(rbind, lapply(ok, function(z) z$point[[q]]))
    lo <- do.call(rbind, lapply(ok, function(z) z$ci_curve[[q]][, 1]))
    hi <- do.call(rbind, lapply(ok, function(z) z$ci_curve[[q]][, 2]))
    exact <- ref$ex[[q]]
    limit <- ref$lim[[q]]
    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      scenario = name, n = n, quantity = q, m = ms,
      true_exact = exact, true_limit = limit,
      cov_exact = colMeans(sweep(lo, 2, exact, "<=") & sweep(hi, 2, exact, ">=")),
      cov_limit = colMeans(sweep(lo, 2, limit, "<=") & sweep(hi, 2, limit, ">=")),
      width = colMeans(hi - lo),
      bias_point = colMeans(point) - exact,
      rmse_point = sqrt(colMeans(sweep(point, 2, exact)^2))
    )
  }

  # Time per replication (seconds, one core)
  t_gen <- vapply(g, function(z) z$time_gen_est, 0)
  t_boot <- vapply(g, function(z) z$time_boot, 0)
  time_rows[[length(time_rows) + 1]] <- data.frame(
    scenario = name, n = n, B = cfg$B, n_rep = n_rep,
    sec_gen_est = mean(t_gen), sec_boot = mean(t_boot),
    sec_total = mean(t_gen + t_boot), cpu_hours = sum(t_gen + t_boot) / 3600
  )
}

params <- do.call(rbind, par_rows)
curves <- do.call(rbind, curve_rows)
tests <- do.call(rbind, test_rows)
times <- do.call(rbind, time_rows)

# =============================================================================
# Part 2: write the tables
# =============================================================================

write_table <- function(x, file, digits) {
  write.table(x, paste0(prefix, file, ".csv"), sep = ";", row.names = FALSE)
  print(xtable(x, digits = digits),
    file = paste0(prefix, file, ".tex"), include.rownames = FALSE, floating = FALSE
  )
}
write_table(params, "table_boot_params", 4)
write_table(tests, "table_boot_tests", 3)
write_table(times, "table_boot_time", 2)
write.table(curves, paste0(prefix, "table_boot_curves.csv"), sep = ";", row.names = FALSE)

# Tables for SC2: every window size in the grid, for each n
for (q in c("rho", "Fdcca")) {
  x <- curves[curves$quantity == q & curves$scenario == "SC2",
              c("n", "m", "true_exact", "true_limit", "cov_exact", "cov_limit", "width")]
  print(xtable(x, digits = c(0, 0, 0, 4, 4, 3, 3, 4)),
    file = paste0(prefix, "table_boot_", q, ".tex"), include.rownames = FALSE, floating = FALSE
  )
}

# =============================================================================
# Part 3: reduced coverage table for the appendix. "Coverage" is the fraction
# of the nrep replications whose bootstrap 95% percentile interval contains
# the true (population) value of the parameter. Scenario, n, parameter, true
# value, coverage and mean width; bias and RMSE are left out because they
# are already in table_boot_params.
# =============================================================================
param_labels <- c(d1 = "d_1", d2 = "d_2", S1 = "S_1^epsilon",
                   S2 = "S_2^epsilon", S12 = "S_{12}^epsilon")

tab <- params[, c("scenario", "n", "param", "true", "coverage", "width", "n_ok", "n_rep")]
tab <- tab[order(tab$scenario, tab$n, match(tab$param, names(param_labels))), ]

write.table(tab, paste0(prefix, "table_coverage_params_SC2.csv"), sep = ";", row.names = FALSE)
print(xtable(tab, digits = c(0, 0, 0, 0, 4, 3, 4, 0, 0)),
  file = paste0(prefix, "table_coverage_params_SC2.tex"), include.rownames = FALSE
)
