# Study I: validation of the asymptotic theory (scenarios SC1 and SC2).
#
# This is the simulation stage of the Study I pipeline. It has two parts, run
# sequentially:
#
#   Part 1 (simulate): bivariate ARFIMA samples are generated with the
#     routine used internally by the parametric bootstrap of the DCCA
#     package (innovations correlated through the Cholesky factor of
#     Sigma), and the raw DFA/DCCA curves are computed for every
#     replication.
#   Part 2 (re-estimate): the point estimates (d1, d2, S1, S2, S12) are
#     computed from those raw curves, using the package's own internal
#     estimation routines, and written to a single tidy CSV.
#
# Note: when test_mode is TRUE, Part 1 writes to "results/TEST_*.csv"
# (see `prefix` below), but Part 2 always reads the full-run filenames
# ("results/result0d0.csv", "results/result1d1.csv"). 
#
# Output (results/):
#   result0d0.csv               SC1, no header, ";"-separated
#                               (rep, rho, d1, d2, m, dfa1, dfa2, dcca12, dccacor)
#   result1d1.csv               SC2, with header (rep, m, dfa, dfa2, dcca12, dccacor)
#   timing_theory.csv           wall-clock time per scenario
#   study1_full_estimates.csv   one row per replication/pair: re-estimated
#                                (d1, d2, S1, S2, S12), used to
#                                build the figures and tables

library(DCCA)
source("00 - utils_parallel.R")

# =============================================================================
# Part 1: simulate SC1 and SC2, save the raw DFA/DCCA curves
# =============================================================================

test_mode <- FALSE # TRUE: tiny run to check the script; FALSE: full study

n <- 2000
nrep <- 1000
mtrunc <- 50000
ms <- seq(10, 300, by = 20)
ds <- c(-0.45, -0.25, 0, 0.25, 0.45)
rhos <- c(-0.5, 0, 0.8)
overlap <- TRUE
# SC2 (ar = -0.4 means phi(z) = 1 + 0.4z, ma = 0.2 means theta(z) = 1 + 0.2z)
sc2 <- list(d = c(0.45, -0.25), ar = -0.4, ma = 0.2, rho = 0.8)
seed0 <- 1234

if (test_mode) {
  n <- 400
  nrep <- 4
  mtrunc <- 2000
  ms <- seq(10, 60, by = 10)
  ds <- c(-0.25, 0, 0.25)
  rhos <- c(0, 0.8)
}

dir.create("results", showWarnings = FALSE)
prefix <- if (test_mode) "results/TEST_" else "results/"
f_sc1 <- paste0(prefix, "result0d0.csv")
f_sc2 <- paste0(prefix, "result1d1.csv")
f_time <- paste0(prefix, "timing_theory.csv")
unlink(c(f_sc1, f_sc2, f_time))

n_cores <- if (test_mode) 2 else max(1, parallel::detectCores() - 1)

sigma_matrix <- function(rho) {
  matrix(c(1, rho, rho, 1), 2, dimnames = list(c("X1", "X2"), c("X1", "X2")))
}

# Row-bind the successful replications and report the failed ones
collect_results <- function(res) {
  failed <- vapply(res, inherits, logical(1), what = "error")
  if (any(failed)) warning(sum(failed), " replication(s) failed and were dropped")
  do.call(rbind, res[!failed])
}

log_time <- function(scenario, rho, elapsed) {
  write.table(
    data.frame(
      scenario = scenario, 
      rho = rho, 
      n = n, 
      nrep = nrep, 
      trunc = mtrunc,
      n_cores = n_cores, 
      elapsed_min = elapsed / 60,
      sec_per_rep = elapsed / nrep, 
      R = R.version.string,
      DCCA = as.character(packageVersion("DCCA"))
    ),
    f_time,
    sep = ";", row.names = FALSE,
    col.names = !file.exists(f_time), append = file.exists(f_time)
  )
}

#********************************************************************
# SCENARIO 1
#    X1 ~ ARFIMA(0,d1,0)
#    X2 ~ ARFIMA(0,d2,0)
#    (e1, e2) ~ N(0,Sigma), Sigma = [[1, rho], [rho, 1]]
# for every rho in `rhos` and every pair (d1, d2) with d1 <= d2
#********************************************************************
psi0 <- lapply(ds, function(d) arfima.expansion(d = d, trunc = mtrunc))
pairs <- which(outer(seq_along(ds), seq_along(ds), "<="), arr.ind = TRUE)

for (k in seq_along(rhos)) {
  rho <- rhos[k]
  cat("SC1: rho =", rho, "\n")
  t0 <- proc.time()
  res <- parallel_map(
    seq_len(nrep),
    function(r) {
      set.seed(seed0 + 1000L * k + r)
      out <- vector("list", nrow(pairs))
      for (p in seq_len(nrow(pairs))) {
        i <- pairs[p, 1]
        j <- pairs[p, 2]
        x <- DCCA:::.boot.simulate(
          n, mtrunc, list(X1 = psi0[[i]], X2 = psi0[[j]]), sigma_matrix(rho)
          )
        tmp <- rhodcca(x[, 1], x[, 2], m = ms, overlap = overlap)
        out[[p]] <- data.frame(
          rep = r, rho = rho, d1 = ds[i], d2 = ds[j], m = ms,
          dfa = tmp$F2dfa1, dfa2 = tmp$F2dfa2,
          dcca12 = tmp$Fdcca, dccacor = tmp$rhodcca
        )
      }
      do.call(rbind, out)
    },
    n_cores,
    export = c("n", "mtrunc", "psi0", "pairs", "ds", "ms", 
               "overlap", "rho", "k", "seed0", "sigma_matrix")
  )
  res <- collect_results(res)
  write.table(res, f_sc1, sep = ";", row.names = FALSE, 
              col.names = FALSE, append = file.exists(f_sc1))
  log_time("SC1", rho, (proc.time() - t0)["elapsed"])
}

#********************************************************************
# SCENARIO 2
#    X1 ~ ARFIMA(1,d1,1)
#    X2 ~ ARFIMA(1,d2,1)
#    (e1, e2) ~ N(0,Sigma), Sigma = [[1, rho], [rho, 1]]
# with
#    phi1(z) = phi2(z)
#    theta1(z) = theta2(z)
#********************************************************************
{
  cat("SC2\n")
  t0 <- proc.time()
  psi2 <- lapply(
    sc2$d, 
    function(d) arfima.expansion(d = d, ar = sc2$ar, ma = sc2$ma, trunc = mtrunc)
    )
  names(psi2) <- c("X1", "X2")
  res <- parallel_map(
    seq_len(nrep),
    function(r) {
      set.seed(seed0 + 900000L + r)
      x <- DCCA:::.boot.simulate(n, mtrunc, psi2, sigma_matrix(sc2$rho))
      tmp <- rhodcca(x[, 1], x[, 2], m = ms, overlap = overlap)
      data.frame(
        rep = r, m = ms, dfa = tmp$F2dfa1, dfa2 = tmp$F2dfa2,
        dcca12 = tmp$Fdcca, dccacor = tmp$rhodcca
      )
    },
    n_cores,
    export = c("n", "mtrunc", "psi2", "ms", 
               "overlap", "sc2", "seed0", "sigma_matrix")
  )
  res <- collect_results(res)
  write.table(res, f_sc2, sep = ";", row.names = FALSE)
  log_time("SC2", sc2$rho, (proc.time() - t0)["elapsed"])
}

# =============================================================================
# Part 2: re-estimate (d1, d2, S1, S2, S12) from the raw curves saved above
# =============================================================================
# Reads results/result0d0.csv and results/result1d1.csv, i.e. the
# full-run output of Part 1 above (test_mode = FALSE).

estimate.dS  <- get(".estimate.dS",  envir = asNamespace("DCCA"))
estimate.S12 <- get(".estimate.S12", envir = asNamespace("DCCA"))

# Re-estimates (d1, d2, S1, S2, S12) for one replication/pair from its saved
# F2dfa1, F2dfa2 and Fdcca curves, using the package's own internal routines.
compute_group <- function(g) {
  object_F <- list(
    F2dfa = cbind(X1 = g$dfa, X2 = g$dfa2),
    Fdcca = matrix(g$dcca12, ncol = 1, dimnames = list(NULL, "X1_X2")),
    m = g$m
  )
  attr(object_F, "series_names") <- c("X1", "X2")
  dS  <- estimate.dS(object_F)
  S12 <- estimate.S12(dS, object_F)
  data.frame(
    d1_hat = unname(dS$d["X1"]), 
    d2_hat = unname(dS$d["X2"]),
    S1_hat = unname(dS$Seps["X1"]), 
    S2_hat = unname(dS$Seps["X2"]),
    S12_hat = unname(S12["X1_X2"])
  )
}

# ---------------------------------------------------------------------------
# SC1: X1 ~ ARFIMA(0, d1, 0), X2 ~ ARFIMA(0, d2, 0); true S1 = S2 = 1,
# true S12 = rho (no ARMA short-memory component in either series).
# ---------------------------------------------------------------------------
sc1 <- read.table("results/result0d0.csv", sep = ";",
                   col.names = c("rep", "rho", "d1", "d2", "m",
                                 "dfa", "dfa2", "dcca12", "dccacor"))

groups1 <- split(sc1, list(sc1$rep, sc1$rho, sc1$d1, sc1$d2), drop = TRUE)
cat("SC1: re-estimating", length(groups1), "groups...\n")
pb1 <- txtProgressBar(min = 0, max = length(groups1), style = 3)
out1_list <- vector("list", length(groups1))
for (i in seq_along(groups1)) {
  g <- groups1[[i]]
  if (nrow(g) >= 2) {
    out1_list[[i]] <- cbind(compute_group(g), rho = g$rho[1],
                            d1_true = g$d1[1], d2_true = g$d2[1])
  }
  setTxtProgressBar(pb1, i)
}
close(pb1)
out1 <- do.call(rbind, out1_list)
out1$S1_true  <- 1
out1$S2_true  <- 1
out1$S12_true <- out1$rho
out1$scenario <- "SC1: ARFIMA(0,d,0)"

# ---------------------------------------------------------------------------
# SC2: X1 ~ ARFIMA(1, 0.45, 1), X2 ~ ARFIMA(1, -0.25, 1), ar = -0.4,
# ma = 0.2, rho = 0.8 (same ARMA short-memory filter on both series).
# ---------------------------------------------------------------------------
sc2_est <- read.table("results/result1d1.csv", sep = ";", header = TRUE)

psi_true <- ((1 + 0.2) / (1 - (-0.4)))^2  # Psi_kk = {theta(1)/phi(1)}^2

groups2 <- split(sc2_est, sc2_est$rep)
cat("SC2: re-estimating", length(groups2), "groups...\n")
pb2 <- txtProgressBar(min = 0, max = length(groups2), style = 3)
out2_list <- vector("list", length(groups2))
for (i in seq_along(groups2)) {
  g <- groups2[[i]]
  if (nrow(g) >= 2) {
    out2_list[[i]] <- cbind(compute_group(g), rho = 0.8, 
                            d1_true = 0.45, d2_true = -0.25)
  }
  setTxtProgressBar(pb2, i)
}
close(pb2)
out2 <- do.call(rbind, out2_list)
out2$S1_true  <- psi_true
out2$S2_true  <- psi_true
out2$S12_true <- psi_true * 0.8
out2$scenario <- "SC2: ARFIMA(1,d,1)"

all_est <- rbind(out1, out2)
write.csv(all_est, "results/study1_full_estimates.csv", row.names = FALSE)
cat("Re-estimation done. Total rows:", nrow(all_est), "\n")

# ---------------------------------------------------------------------------
# Numeric summary (mean bias by true d, pooled over scenario/rho)
# ---------------------------------------------------------------------------
long_dS <- rbind(
  data.frame(scenario = all_est$scenario,
             d_true = all_est$d1_true, d_hat = all_est$d1_hat,
             S_true = all_est$S1_true, S_hat = all_est$S1_hat),
  data.frame(scenario = all_est$scenario,
             d_true = all_est$d2_true, d_hat = all_est$d2_hat,
             S_true = all_est$S2_true, S_hat = all_est$S2_hat)
)
long_dS$d_bias <- long_dS$d_hat - long_dS$d_true
long_dS$S_bias <- long_dS$S_hat - long_dS$S_true
all_est$S12_bias <- all_est$S12_hat - all_est$S12_true

cat("\n--- mean bias of d-hat and S-hat by true d and scenario ---\n")
print(aggregate(cbind(d_bias, S_bias) ~ d_true + scenario, long_dS, mean))

cat("\n--- mean bias of S12-hat by (d1, d2, rho) ---\n")
print(aggregate(S12_bias ~ d1_true + d2_true + rho + scenario, all_est, mean))
