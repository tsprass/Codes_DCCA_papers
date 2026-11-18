# Parametric bootstrap experiments: all figures.
#
# This script only reads the summaries written to tables/ -- no
# simulation is run here, so it can be run on a different machine, without the
# DCCA package, as long as those CSVs are available. It has two independent
# parts:
#
#   Part 1 (curve coverage):     tables/table_boot_curves.csv
#     -> figures/boot_coverage_SC2.pdf
#   Part 2 (parameter coverage): tables/table_boot_params.csv
#     -> figures/boot_coverage_params_SC2.pdf
#
# Both figures use the Okabe-Ito palette (colorblind-safe, distinguishable in
# grayscale print).

suppressPackageStartupMessages(library(ggplot2))

test_mode <- FALSE # must match test_mode of the simulation and table stages
prefix_in <- if (test_mode) "tables/TEST_" else "tables/"
prefix_out <- if (test_mode) "figures/TEST_" else "figures/"
dir.create("figures", showWarnings = FALSE)

nominal <- 0.95 # 1 - alpha, with alpha = 0.05 as set in the simulation stage

# =============================================================================
# Part 1: coverage of the exact expectation versus window size, SC2
# One panel per curve, one line per n.
# =============================================================================

curves <- read.csv(paste0(prefix_in, "table_boot_curves.csv"), sep = ";")

curve_names <- c("F2dfa1", "F2dfa2", "Fdcca", "rho")
# plotmath expressions (not LaTeX) for the facet strip labels
curve_expr <- c('F["1,DFA"]^2*(m)', 'F["2,DFA"]^2*(m)', 'F["DCCA"]*(m)', 'rho["DCCA"]*(m)')

sc2 <- curves[curves$scenario == "SC2", ]
sc2$quantity <- factor(sc2$quantity, levels = curve_names, labels = curve_expr)
sc2$n <- factor(sc2$n)
ms_grid <- sort(unique(sc2$m)) # the simulated window sizes, for the x-axis ticks

# Line type and point shape also vary by n, so the lines are still told apart
# in black and white.
p_curves <- ggplot(sc2, aes(x = m, y = cov_exact, color = n, linetype = n, shape = n)) +
  geom_hline(yintercept = nominal, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  facet_wrap(~quantity, nrow = 2, ncol = 2, labeller = label_parsed) +
  scale_x_continuous(breaks = ms_grid) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  labs(x = expression(m), y = "empirical coverage",
       color = expression(n), linetype = expression(n), shape = expression(n)) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave(paste0(prefix_out, "boot_coverage_SC2.pdf"), p_curves, width = 7, height = 6.5)
cat("Part 1 done:", paste0(prefix_out, "boot_coverage_SC2.pdf"), "\n")

# =============================================================================
# Part 2: coverage of the five parameters versus n
# One panel per parameter, one line per scenario, dashed line at the nominal
# level.
# =============================================================================

params <- read.csv(paste0(prefix_in, "table_boot_params.csv"), sep = ";")

param_names <- c("d1", "d2", "S1", "S2", "S12")
param_expr <- c("hat(d)[1]", "hat(d)[2]", "hat(S)[1]^epsilon",
                "hat(S)[2]^epsilon", "hat(S)[12]^epsilon")
params$param <- factor(params$param, levels = param_names, labels = param_expr)
params$n <- factor(params$n)

scenario_labels <- c(SC2 = expression(SC2), SC2_null = expression(SC2[0]))

p_params <- ggplot(params, aes(x = n, y = coverage, color = scenario,
                               shape = scenario, group = scenario)) +
  geom_hline(yintercept = nominal, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~param, nrow = 2, ncol = 3, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c(SC2 = "#0072B2", SC2_null = "#D55E00"),
                     labels = scenario_labels) +
  scale_shape_discrete(labels = scenario_labels) +
  labs(x = expression(n), y = "empirical coverage", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave(paste0(prefix_out, "boot_coverage_params_SC2.pdf"), p_params, width = 8, height = 5.5)
cat("Part 2 done:", paste0(prefix_out, "boot_coverage_params_SC2.pdf"), "\n")
