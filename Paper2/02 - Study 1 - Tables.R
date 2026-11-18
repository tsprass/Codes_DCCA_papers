# Study I: all tables (mean bias/SD summary and percentile-coverage cells).
#
# This script only reads results/study1_full_estimates.csv (produced by
# the simulation stage) -- no new simulation or re-estimation is run here. It
# has two independent parts:
#
#   Part 1 (bias/SD summary):     tables/table_bias_dS_summary.csv
#     Mean bias and SD of d-hat and S-hat, pooled over rho and series, by
#     true d and scenario -- this is exactly Table 1 in the paper.
#   Part 2 (percentile coverage): tables/table_percentile_coverage_cells.{csv,tex}
#     One row per (scenario, d1_true, d2_true, rho, quantity), fully
#     disaggregated -- the reference table behind the coverage heatmap.

suppressPackageStartupMessages(library(xtable))

dir.create("tables", showWarnings = FALSE)

all_est <- read.csv("results/study1_full_estimates.csv")
cat("Rows read:", nrow(all_est), "\n")

# =============================================================================
# Part 1: mean bias and SD of d-hat and S-hat by true d and scenario
# (pooled over rho and series -- the content of Table 1 in the paper)
# =============================================================================

long_dS <- rbind(
  data.frame(scenario = all_est$scenario, series = "X1", rho = all_est$rho,
             d_true = all_est$d1_true, d_hat = all_est$d1_hat,
             S_true = all_est$S1_true, S_hat = all_est$S1_hat),
  data.frame(scenario = all_est$scenario, series = "X2", rho = all_est$rho,
             d_true = all_est$d2_true, d_hat = all_est$d2_hat,
             S_true = all_est$S2_true, S_hat = all_est$S2_hat)
)
long_dS$d_bias <- long_dS$d_hat - long_dS$d_true
long_dS$S_bias <- long_dS$S_hat - long_dS$S_true

all_est$S12_bias <- all_est$S12_hat - all_est$S12_true

cat("\n--- mean bias and SD of d-hat and S-hat by true d and scenario ---\n")
bias_dS <- merge(
  merge(aggregate(d_bias ~ d_true + scenario, long_dS, mean),
        aggregate(d_hat  ~ d_true + scenario, long_dS, sd)),
  merge(aggregate(S_bias ~ d_true + scenario, long_dS, mean),
        aggregate(S_hat  ~ d_true + scenario, long_dS, sd))
)
names(bias_dS)[names(bias_dS) == "d_hat"] <- "sd_d"
names(bias_dS)[names(bias_dS) == "S_hat"] <- "sd_S"
names(bias_dS)[names(bias_dS) == "d_bias"] <- "bias_d"
names(bias_dS)[names(bias_dS) == "S_bias"] <- "bias_S"
bias_dS <- bias_dS[order(bias_dS$scenario, bias_dS$d_true),
                    c("scenario", "d_true", "bias_d", "sd_d", "bias_S", "sd_S")]
write.csv(bias_dS, "tables/table_bias_dS_summary.csv", row.names = FALSE)
print(bias_dS, digits = 3, row.names = FALSE)
cat("Wrote tables/table_bias_dS_summary.csv (source of Table 1 in the paper)\n")

cat("\n--- mean bias of S12-hat by (d1, d2, rho) ---\n")
print(aggregate(S12_bias ~ d1_true + d2_true + rho + scenario, all_est, mean))

# =============================================================================
# Part 2: percentile-coverage table, fully disaggregated by
# (scenario, d1_true, d2_true, rho, quantity) -- no pooling across cells or
# across quantities. Each cell's own replications are treated as if they
# were a bootstrap distribution centered on the true parameter, and the
# 2.5%-97.5% percentile interval is checked for coverage of the truth.
# =============================================================================

pct_interval <- function(x, probs = c(0.025, 0.975)) {
  q <- quantile(x, probs = probs, names = FALSE)
  c(lo = q[1], hi = q[2])
}

quantities <- list(
  d1  = c(hat = "d1_hat",  true = "d1_true"),
  d2  = c(hat = "d2_hat",  true = "d2_true"),
  S1  = c(hat = "S1_hat",  true = "S1_true"),
  S2  = c(hat = "S2_hat",  true = "S2_true"),
  S12 = c(hat = "S12_hat", true = "S12_true")
)

groups <- split(all_est, list(all_est$scenario, all_est$d1_true, all_est$d2_true, all_est$rho), drop = TRUE)

cells <- do.call(rbind, lapply(groups, function(g) {
  rows <- lapply(names(quantities), function(qname) {
    cols <- quantities[[qname]]
    q <- pct_interval(g[[cols["hat"]]])
    true_val <- g[[cols["true"]]][1]
    data.frame(
      scenario = g$scenario[1], d1_true = g$d1_true[1], d2_true = g$d2_true[1],
      rho = g$rho[1], quantity = qname, n_rep = nrow(g),
      mean_hat = mean(g[[cols["hat"]]]), lo = q["lo"], hi = q["hi"],
      true_val = true_val,
      covered = (q["lo"] <= true_val) & (true_val <= q["hi"])
    )
  })
  do.call(rbind, rows)
}))
cells <- cells[order(cells$scenario, cells$rho, cells$d1_true, cells$d2_true,
                     match(cells$quantity, names(quantities))), ]
rownames(cells) <- NULL

write.table(cells, "tables/table_percentile_coverage_cells.csv", sep = ";", row.names = FALSE)
print(xtable(cells, digits = c(0, 0, 2, 2, 2, 0, 0, 4, 4, 4, 4, 0)),
  file = "tables/table_percentile_coverage_cells.tex", include.rownames = FALSE
)
cat("Wrote tables/table_percentile_coverage_cells.{csv,tex} --", nrow(cells), "rows",
    "(", length(groups), "cells x", length(quantities), "quantities )\n")

cat("\n--- coverage rate by (scenario, quantity), over every disaggregated cell ---\n")
print(aggregate(covered ~ scenario + quantity, cells, mean))

cat("\n--- SC2 (single cell: d1=0.45, d2=-0.25, rho=0.8) ---\n")
print(cells[cells$scenario == "SC2: ARFIMA(1,d,1)", ], row.names = FALSE, digits = 4)
