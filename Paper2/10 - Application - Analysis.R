library(DCCA)
source("09 - Application - Functions.R")

ptm <- proc.time()

############################################
#
#   Application: BTC + 4 stock indices
#   ------------------------------------
#   Rewritten to use the package's estimate()/bootstrap()/update() API
#   (this code predates that API; the previous version hand-rolled its
#   own bootstrap loops, referenced helper functions that were never
#   defined -- .C_func, estimate_parameters -- and its plotting code
#   referenced a field, m_grid, that the bootstrap never returned).
#
#   With 5 series (10 pairs), estimate(..., type = "joint") gives d_k,
#   S_k^eps for every series and S_{ij}^eps for every pair; bootstrap()
#   dispatches to the two-stage procedure (MB1: marginal bootstrap of
#   each series once; MB2: pairwise bootstrap of S_{ij}^eps with the
#   marginals held fixed), automatically excluding any pair that
#   violates the Cauchy-Schwarz inequality (with a warning) instead of
#   the manual check the old code did by hand; update() applies the
#   zero-in-CI restriction from Prass and Pumi (2026) automatically.
#
#   Output:
#     results/application_result.rds -- the full estimate()/bootstrap()/
#       update() result, so the (expensive) bootstrap does not need to
#       be rerun just to change a plot or a table.
#     tables/application_parameters.csv/.tex -- d_k, S_k^eps per series.
#     tables/application_pairs.csv/.tex -- S_{ij}^eps per pair.
#     figures/application_dfa_panel.pdf -- marginal DFA plots.
#     figures/application_pair_matrix.pdf -- F_DCCA(m)/rho_DCCA(m) matrix.
#
############################################

# ----- Quick test mode (small n_boot, to validate the script) -----
test_mode <- TRUE

files <- list(
  BTC  = "data/BTC-USD.csv",
  BVSP = "data/BVSP.csv",
  DJI  = "data/DJI.csv",
  GSPC = "data/GSPC.csv",
  IXIC = "data/IXIC.csv"
)
series_names <- names(files)

# ----- Load and merge the data, compute log-returns -----
data <- load_data(files)
rets <- data.frame(apply(data[, -1], 2, function(x) diff(log(x)) * 100))
head(rets)
cat(sprintf("\nSample size (n) after merging and differencing: %d\n", nrow(rets)))

if (!dir.exists("results")) dir.create("results")
if (!dir.exists("tables")) dir.create("tables")
if (!dir.exists("figures")) dir.create("figures")

# ----- Estimation grid (identical to the Simulation section) -----
ms <- seq(100, 300, by = 10)

n_boot <- if (test_mode) 20 else 1000
if (test_mode) cat("*** test_mode = TRUE: running a reduced bootstrap just to validate the script ***\n")

# ----- 1. Point estimates: d_k, S_k^eps (every series) and S_{ij}^eps (every pair) -----
est <- DCCA::estimate(as.matrix(rets), type = "joint", m = ms, nu = 0, overlap = TRUE)
print(est)

# ----- 2. Bootstrap (MB1 marginal + MB2 pairs; boot_control's m/nu/overlap -----
# are always aligned to the estimates by bootstrap(), so they never need
# to be set here -- see .boot.control.estimates() in the package)
boot_control <- DCCA::boot.control(
  n_boot = n_boot, parallel = TRUE,
  n_cores = max(1, parallel::detectCores() - 1)
)
result <- DCCA::bootstrap(est, boot_control = boot_control)
print(result)

# ----- 3. Apply the zero-in-CI restriction (Prass and Pumi, 2026), if any -----
# (update() is stats::update's generic; DCCA only registers methods for
# it, so it is never called as DCCA::update())
result <- update(result, rerun_boot = TRUE)
str(result$zero_in_ci)

saveRDS(result, "results/application_result.rds")

# ----- 4. Summaries -----
smry <- summary(result)

write.table(smry$marginal$parameters, "tables/application_parameters.csv", sep = ";", row.names = FALSE)
print(
  xtable::xtable(smry$marginal$parameters, digits = 4),
  file = "tables/application_parameters.tex", include.rownames = FALSE
)

pairs_table <- do.call(rbind, lapply(names(smry$pairs), function(pair) {
  cbind(pair = pair, smry$pairs[[pair]]$S12eps)
}))
write.table(pairs_table, "tables/application_pairs.csv", sep = ";", row.names = FALSE)
print(
  xtable::xtable(pairs_table, digits = 4),
  file = "tables/application_pairs.tex", include.rownames = FALSE
)

excluded_pairs <- attr(result$bootstrap$pairs, "excluded_pairs")
if (length(excluded_pairs)) {
  cat("\nPairs excluded (Cauchy-Schwarz violated):", paste(excluded_pairs, collapse = ", "), "\n")
}

# ----- 5. Figures -----
pdf("figures/application_dfa_panel.pdf", width = 20, height = 5)
grid.draw(plot_marginal_panel(smry$marginal$dfa, series_names))
dev.off()

pdf("figures/application_pair_matrix.pdf", width = 20, height = 20)
grid.draw(plot_pair_matrix(smry$pairs, series_names))
dev.off()

if (test_mode) {
  cat("\n*** test_mode = TRUE: set test_mode <- FALSE for the full bootstrap (n_boot = 1000). ***\n")
}

cat(sprintf("\nTotal script time: %.1f min\n", (proc.time() - ptm)["elapsed"] / 60))
