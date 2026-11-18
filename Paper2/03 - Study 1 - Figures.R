# Study I: all figures (SC1, SC2 and the derived bias/coverage figures).
#
# This script only reads the CSVs already produced by the simulation stage --
# it does not simulate or re-estimate anything, so it can be run on a
# different machine, without the DCCA package, as long as those CSVs are
# available. It has four parts, each independent of the others:
#
#   Part 1 (SC1 raw curves):    results/result0d0.csv
#     -> figures/rho{-05,0,08}DCCA.pdf, figures/dfa1.pdf
#   Part 2 (SC2 raw curves):    results/result1d1.csv
#     -> figures/SC2.pdf
#   Part 3 (bias boxplots):     results/study1_full_estimates.csv
#     -> figures/study1_bias_{d,S,S12}.pdf
#   Part 4 (coverage heatmap):  results/study1_full_estimates.csv
#     -> figures/study1_coverage_heatmap.pdf
#
# Parts 1-2 need the DCCA package (for rhoE/rhoLim, used to draw the
# theoretical expected value and asymptotic limit curves); Parts 3-4 do
# not.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(tidyr)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper functions ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gk.arfima <- function(m, d = 0, ar = NULL, ma = NULL, tau = 1, trunc = 50000) {
  toeplitz(DCCA::arfima.acf(
    lag.max = m, d = d, ar = ar, ma = ma, var.innov = tau, trunc = trunc
  )$acf)
}

G12.arfima <- function(m, d1 = 0, ar1 = NULL, ma1 = NULL,
                       d2 = 0, ar2 = NULL, ma2 = NULL, tau12 = 0, trunc = 50000) {
  # Create index matrix: lag = j - i, then convert to CCF index
  # CCF indices: 1 = lag -m, m+1 = lag 0, 2m+1 = lag +m
  # So index = (j - i) + (m + 1)
  gamma.12 <- DCCA::arfima.ccf(
    lag.max = m, d1 = d1, ar1 = ar1, ma1 = ma1,
    d2 = d2, ar2 = ar2, ma2 = ma2, cov.innov = tau12, trunc = trunc
  )$ccf
  index_matrix <- outer(1:(m + 1), 1:(m + 1), function(i, j) j - i + (m + 1))
  matrix(gamma.12[index_matrix], nrow = m + 1, ncol = m + 1)
}

main.plot <- function(plot_data, limit_data, fs = 12) {
  ggplot() +
    # Empirical boxplots
    geom_boxplot(data = plot_data, 
                 aes(x = factor(m), y = empirical, fill = fill), 
                 alpha = 0.7, width = 0.5, outlier.size = 0.7) +
    scale_fill_manual(values = c("dfa1" = "#9B59B6", "dfa2" = "#E74C3C",
                                 "dcca" = "#27AE60", "rho" = "#F1C40F"), 
                      guide = "none") +
    # Theoretical expected values
    geom_line(data = limit_data, aes(x = factor(m), y = expected, group = 1),
              color = "blue", linewidth = 0.8, linetype = "solid") +
    # Asymptotic limits
    geom_line(data = limit_data, aes(x = factor(m), y = limit, group = 1),
              color = "red", linewidth = 0.8, linetype = "dashed") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5, size = fs),
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = fs),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) + labs(x = "", y = "")
}

# Get empirical data for a specific (d1, d2) pair
get_data <- function(data, rho_val, d1_val, d2_val, ms) {
  # Always search for d1 <= d2 in the data
  dd1 <- ifelse(d1_val <= d2_val, d1_val, d2_val)
  dd2 <- ifelse(d1_val <= d2_val, d2_val, d1_val)
  data <- data %>% filter(rho == rho_val & d1 == dd1 & d2 == dd2 & m %in% ms)
  plot_data <- NULL
  
  # DCCA for upper triangle (including diagonal)
  if (d1_val >= d2_val) {
    plot_data <- rbind(
      plot_data,
      data %>%
        mutate(empirical = dcca12, type = "F[DCCA]", fill = "dcca") %>%
        select(rep, m, empirical, type, fill))
  }
  
  # rhoDCCA for lower triangle (including diagonal)
  if (d1_val <= d2_val) {
    plot_data <- rbind(
      plot_data,
      data %>%
        mutate(empirical = dccacor, type = "rho[DCCA]", fill = "rho") %>%
        select(rep, m, empirical, type, fill))
  }
  
  return(plot_data)
}

# Get theoretical limits for a specific (d1, d2) pair
get_limits <- function(data, rho_val, d1_val, d2_val, ms) {
  # Always search for d1 <= d2 in the data
  dd1 <- ifelse(d1_val <= d2_val, d1_val, d2_val)
  dd2 <- ifelse(d1_val <= d2_val, d2_val, d1_val)
  data <- data %>% filter(rho == rho_val & d1 == dd1 & d2 == dd2 & m %in% ms)
  
  plot_data <- NULL
  
  # DCCA limits for upper triangle
  if (d1_val >= d2_val) {
    plot_data <- rbind(
      plot_data,
      data %>%
        mutate(expected = Es12, limit = Lim12, type = "F[DCCA]") %>%
        select(m, expected, limit, type))
  }
  
  # rhoDCCA limits for lower triangle
  if (d1_val <= d2_val) {
    plot_data <- rbind(
      plot_data,
      data %>%
        mutate(expected = Es, limit = Lim, type = "rho[DCCA]") %>%
        select(m, expected, limit, type))
  }
  
  return(plot_data)
}

main.plot.dcca <- function(plot_data, limit_data, d1, d2, row_idx,
                           col_idx, n, fs = 12, yb = 4) {
  # create base plot
  p <- main.plot(plot_data, limit_data, fs = fs)
  
  # X-axis label only for bottom row
  if (row_idx != n) p <- p + theme(axis.text.x = element_blank())
  
  # Y-axis label only for first column (d2 values)
  if (col_idx == 1) {
    p <- p + labs(y = bquote(d[2] == .(d2))) +
      theme(axis.title.y = element_text(hjust = 0.5, face = "bold", size = fs))
  } else {
    p <- p + labs(y = "")
  }
  
  # Title with d1 value only for top row
  if (row_idx == 1) {
    p <- p + labs(title = bquote(d[1] == .(d1))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = fs))
  }
  p + scale_y_continuous(n.breaks = yb)
}

# Create final matrix plot with proper labels
create_final_matrix_plot <- function(
    rho_val, wd = 20, ht = 20, fs = 12, l1 = 6, l2 = length(ms), yb = 4, id = NULL) {
  n <- length(ds)
  plot_list <- list()
  
  # Create all individual plots
  for (row_idx in 1:n) { # rows: d2
    for (col_idx in 1:n) { # columns: d1
      d2 <- ds[row_idx]
      d1 <- ds[col_idx]
      
      lout <- ifelse(d1 == d2, l1, l2)
      my_ms <- ms[seq(1, length(ms), length.out = lout)]
      
      # Get data for this cell
      plot_data <- get_data(rresult, rho_val, d1, d2, my_ms) %>%
        mutate(type = factor(type, levels = c("rho[DCCA]", "F[DCCA]")))
      limit_data <- get_limits(limits, rho_val, d1, d2, my_ms) %>%
        mutate(type = factor(type, levels = c("rho[DCCA]", "F[DCCA]")))
      
      # Create base plot
      p <- main.plot.dcca(plot_data, limit_data, d1, d2, row_idx, col_idx, n, fs, yb)
      # Add faceting for diagonal cells (d1 = d2)
      if (d1 == d2) {
        p <- p + facet_wrap(. ~ type, scales = "free", labeller = label_parsed)
      }
      plot_list[[paste(row_idx, col_idx, sep = "_")]] <- p
    }
  }
  
  # Convert list to matrix format for grid.arrange
  grid_plots <- list()
  for (i in 1:(n * n)) {
    row <- ceiling(i / n)
    col <- i - (row - 1) * n
    grid_plots[[i]] <- plot_list[[paste(row, col, sep = "_")]]
  }
  
  # Create layout matrix
  layout_mat <- matrix(1:(n * n), nrow = n, ncol = n, byrow = TRUE)
  
  
  # Arrange all components
  final_grob <- arrangeGrob(
    grobs = grid_plots,
    nrow = n,
    ncol = n,
    layout_matrix = layout_mat
  )
  
  # Save to PDF
  filename <- paste0("figures/rho", gsub("\\.", "", rho_val), "DCCA.pdf")
  pdf(filename, width = wd, height = ht)
  grid.draw(final_grob)
  dev.off()
  
  return(invisible(final_grob))
}

pct_interval <- function(x, probs = c(0.025, 0.975)) {
  q <- quantile(x, probs = probs, names = FALSE)
  c(lo = q[1], hi = q[2])
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 1: SC1 raw curves (results/result0d0.csv) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters
m_min <- 10
m_max <- 300
ms <- seq(m_min, m_max, 20)
bks <- seq(1, length(ms), 4)
ds <- c(-0.45, -0.25, 0, 0.25, 0.45)
rhos <- c(-0.5, 0, 0.8)

# Read data
rresult <- read.table("results/result0d0.csv", sep = ";", header = TRUE)
colnames(rresult) <- c("rep", "rho", "d1", "d2", "m",
                       "dfa1", "dfa2", "dcca12", "dccacor")
head(rresult)

# Expected values and limits
# Create limits data - only for d1 <= d2
limits <- NULL
for (rho in rhos) {
  for (i in 1:length(ds)) {
    for (j in i:length(ds)) { # j starts at i to ensure d1 <= d2
      d1 <- ds[i]
      d2 <- ds[j]
      Gamma1 <- Gk.arfima(m_max, d = d1, tau = 1)
      Gamma2 <- Gk.arfima(m_max, d = d2, tau = 1)
      Gamma12 <- G12.arfima(m_max, d1 = d1, d2 = d2, tau12 = rho)
      E <- DCCA::rhoE(ms, G1 = Gamma1, G2 = Gamma2, G12 = Gamma12)
      lim <- DCCA::rhoLim(m = ms, d1 = d1, d2 = d2,
                          var.innov1 = 1, var.innov2 = 1, cov.innov = rho)
      limits <- rbind(
        limits,
        data.frame(rho = rho, d1 = d1, d2 = d2, m = ms,
                   Es1 = E$F2dfa1, Lim1 = lim$F2dfa1,
                   Es12 = E$Fdcca, Lim12 = lim$Fdcca,
                   Es = E$rhodcca, Lim = lim$rhodcca))
    }
  }
}


# Create plots for all rho values
for (rho in rhos) {
  create_final_matrix_plot(rho, wd = 15, ht = 10, fs = 10, l1 = 6, l2 = 12, yb = 4)
}

p <- main.plot(
  plot_data = rresult %>% filter(rho == 0 & d2 == 0.45) %>% mutate(empirical = dfa1, fill = "dfa1"),
  limit_data = limits %>% filter(rho == 0 & d2 == 0.45) %>% mutate(expected = Es1, limit = Lim1)
  ) + 
  facet_wrap(~d1, ncol = 5, labeller = label_bquote(d[1] == .(d1)), scales = "free") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 14),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none",
    legend.box = "horizontal"
  )

# Save the plot
pdf("figures/dfa1.pdf", width = 25, height = 5)
print(p)
dev.off()

cat("Part 1 done: figures/rho*DCCA.pdf, figures/dfa1.pdf\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 2: SC2 raw curves (results/result1d1.csv) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Scenario 2 parameters
d1_sc2 <- 0.45
d2_sc2 <- -0.25
rho_sc2 <- 0.8
phi_sc2 <- -0.4
theta_sc2 <- 0.2
ms_sc2 <- seq(10, 300, 20)
m_max_sc2 <- max(ms_sc2)

# Read simulation results
result_sc2 <- read.table("results/result1d1.csv", sep = ";", header = TRUE)
head(result_sc2)

# Expected values and limits
Gamma1_sc2 <- Gk.arfima(m = m_max_sc2, d = d1_sc2, ar = phi_sc2, ma = theta_sc2, tau = 1)
Gamma2_sc2 <- Gk.arfima(m = m_max_sc2, d = d2_sc2, ar = phi_sc2, ma = theta_sc2, tau = 1)
Gamma12_sc2 <- G12.arfima(
  m = m_max_sc2, d1 = d1_sc2, ar1 = phi_sc2, ma1 = theta_sc2,
  d2 = d2_sc2, ar2 = phi_sc2, ma2 = theta_sc2, tau12 = rho_sc2
)
E_sc2 <- DCCA::rhoE(m = ms_sc2, G1 = Gamma1_sc2, G2 = Gamma2_sc2, G12 = Gamma12_sc2)
lim_sc2 <- DCCA::rhoLim(
  m = ms_sc2, d1 = d1_sc2, ar1 = phi_sc2, ma1 = theta_sc2, var.innov1 = 1,
  d2 = d2_sc2, ar2 = phi_sc2, ma2 = theta_sc2, var.innov2 = 1, cov.innov = rho_sc2
)

limits_sc2 <- data.frame(
  m = ms_sc2,
  Es1 = E_sc2$F2dfa1, Lim1 = lim_sc2$F2dfa1,
  Es2 = E_sc2$F2dfa2, Lim2 = lim_sc2$F2dfa2,
  Es12 = E_sc2$Fdcca, Lim12 = lim_sc2$Fdcca,
  Es = E_sc2$rhodcca, Lim = lim_sc2$rhodcca
)

p1 <- main.plot(
  plot_data = result_sc2 %>% mutate(empirical = dfa, fill = "dfa1"),
  limit_data = limits_sc2 %>% mutate(expected = Es1, limit = Lim1)) + 
  ggtitle(~ {F["1,DFA"]^2}(m)) + 
  theme(
  axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(face = "bold", size = 20))

p2 <- main.plot(
  plot_data = result_sc2 %>% mutate(empirical = dfa2, fill = "dfa2"),
  limit_data = limits_sc2 %>% mutate(expected = Es2, limit = Lim2)) + 
  ggtitle(~ {F["2,DFA"]^2}(m)) + 
  theme(
  axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(face = "bold", size = 20))

p3 <- main.plot(
  plot_data = result_sc2 %>% mutate(empirical = dcca12, fill = "dcca"),
  limit_data = limits_sc2 %>% mutate(expected = Es12, limit = Lim12)) + 
  ggtitle(~ {F["DCCA"]}(m)) + 
  theme(
  axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(face = "bold", size = 20))

p4 <- main.plot(
  plot_data = result_sc2 %>% mutate(empirical = dccacor, fill = "rho"),
  limit_data = limits_sc2 %>% mutate(expected = Es, limit = Lim)) + 
  ggtitle(~ {rho["DCCA"]}(m)) + 
  theme(
  axis.text.x = element_text(angle = 90, hjust = 0.5, size = 10),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(face = "bold", size = 20))

p_final_sc2 <- (p1 | p2 | p3 | p4)

pdf("figures/SC2.pdf", width = 25, height = 5)
print(p_final_sc2)
dev.off()

cat("Part 2 done: figures/SC2.pdf\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 3: bias boxplots (results/study1_full_estimates.csv) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_est <- read.csv("results/study1_full_estimates.csv")
cat("Rows read:", nrow(all_est), "\n")

# Reshape to long format: one row per (series, replication) for d and S,
# so both series and both scenarios can be shown side by side.
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
long_dS$d_true_f <- factor(long_dS$d_true)

all_est$S12_bias <- all_est$S12_hat - all_est$S12_true

cat("Building bias figures...\n")
theme_set(theme_bw(base_size = 11))

# Legend labels for the scenarios in mathematical notation.
scenario_math <- function(x) {
  parse(text = sub("^(SC[0-9]+): ARFIMA\\((.*)\\)$", "\\1*\":\"~ARFIMA(\\2)", x))
}

# Aggregated over rho (checked separately that the bias of the marginal
# estimates d-hat and S-hat does not depend on rho, as expected, since a
# marginal series' DFA curve does not involve its partner series).
p_d <- ggplot(long_dS, aes(x = d_true_f, y = d_bias, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.3,
               position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = scenario_math) +
  labs(x = expression(paste("true ", d)), y = expression(hat(d) - d), fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/study1_bias_d.pdf", p_d, width = 7, height = 4.5)

p_S <- ggplot(long_dS, aes(x = d_true_f, y = S_bias, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.3,
               position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = scenario_math) +
  labs(x = expression(paste("true ", d)), y = expression(hat(S)^epsilon - S^epsilon), fill = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/study1_bias_S.pdf", p_S, width = 7, height = 4.5)

# Grid: rho in rows, d1 in columns, d2 on the x-axis within each panel
# (boxplot keeps the within-cell variability a heatmap of means would
# hide). scales/space = "free_x" lets each column show only the d2 values
# that actually occur for that d1 (SC1 only has d1 <= d2, so larger d1
# columns have fewer boxes) instead of empty panels.
d1_levels <- sort(unique(all_est$d1_true))
all_est$d1 <- factor(all_est$d1_true, levels = d1_levels,
                     labels = paste0("d[1] == ", d1_levels))
all_est$d2 <- factor(all_est$d2_true, levels = sort(unique(all_est$d2_true)))
rho_levels <- sort(unique(all_est$rho))
all_est$rho_f <- factor(all_est$rho, levels = rho_levels,
                        labels = paste0("rho == ", rho_levels))

p_S12 <- ggplot(all_est, aes(x = d2, y = S12_bias, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3,
               position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = scenario_math) +
  facet_grid(rows = vars(rho_f), cols = vars(d1),
             labeller = label_parsed,
             scales = "free_x", space = "free_x") +
  labs(x = expression(d[2]), y = expression(hat(S)[12]^epsilon - S[12]^epsilon), fill = NULL) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave("figures/study1_bias_S12.pdf", p_S12, width = 11, height = 6)
cat("Part 3 done: figures/study1_bias_{d,S,S12}.pdf\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 4: percentile-coverage heatmap (results/study1_full_estimates.csv) ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# For every SC1 cell (d1_true, d2_true, rho) individually -- no pooling --
# treats that cell's replications as if they were a bootstrap distribution
# centered on the true parameter, and checks whether the 2.5%-97.5%
# percentile interval covers the true value of d1, d2, S1, S2 and S12. The underlying
# numeric table is tables/table_percentile_coverage_cells.csv.

quantities <- list(
  d1  = c(hat = "d1_hat",  true = "d1_true"),
  d2  = c(hat = "d2_hat",  true = "d2_true"),
  S1  = c(hat = "S1_hat",  true = "S1_true"),
  S2  = c(hat = "S2_hat",  true = "S2_true"),
  S12 = c(hat = "S12_hat", true = "S12_true")
)

groups_cov <- split(all_est, list(all_est$scenario, all_est$d1_true, 
                                  all_est$d2_true, all_est$rho), drop = TRUE)

cells <- do.call(rbind, lapply(groups_cov, function(g) {
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

sc1_cov <- cells[cells$scenario == "SC1: ARFIMA(0,d,0)", ]
sc1_cov$pair <- sprintf("(%.2f, %.2f)", sc1_cov$d1_true, sc1_cov$d2_true)
pair_order <- unique(sc1_cov$pair[order(sc1_cov$d1_true, sc1_cov$d2_true)])
sc1_cov$pair <- factor(sc1_cov$pair, levels = pair_order)
sc1_cov$quantity <- factor(sc1_cov$quantity, levels = rev(names(quantities)))
rho_levels_sc1 <- sort(unique(sc1_cov$rho))
sc1_cov$rho_f <- factor(sc1_cov$rho, levels = rho_levels_sc1,
                        labels = paste0("rho == ", rho_levels_sc1))

p_heat <- ggplot(sc1_cov, aes(x = pair, y = quantity, fill = covered)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~rho_f, ncol = 1, labeller = label_parsed) +
  scale_y_discrete(labels = c(
    d1 = expression(d[1]), d2 = expression(d[2]),
    S1 = expression(S[1]^epsilon), S2 = expression(S[2]^epsilon),
    S12 = expression(S[12]^epsilon))) +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"),
                    labels = c(`TRUE` = "covers truth", `FALSE` = "misses truth")) +
  labs(x = expression((d[1]*", "*d[2])), y = NULL, fill = NULL,
       title = "SC1: does the 95% percentile interval of the replications cover the true value?") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "bottom",
    panel.grid = element_blank()
  )
ggsave("figures/study1_coverage_heatmap.pdf", p_heat, width = 10, height = 6.5)
cat("Part 4 done: figures/study1_coverage_heatmap.pdf\n")
