# =============================================================================
# SCRIPT 06: SIMULATION AND FIGURE - OVERLAPPING VS. NON-OVERLAPPING BOXES
# =============================================================================
# Paper: Prass, T. S. and Pumi, G. (2021). On the behavior of the DFA and DCCA
#        in trend-stationary processes. Journal of Multivariate Analysis,
#        182, 104703.
#
# Compares the finite-sample variability of rho_DCCA(m) computed with
# overlapping versus non-overlapping boxes (n = 100, m in {3, ..., 25}, 1,000
# replications) under Scenario 3a (MA(20) cross-correlation structure).
#
# Output:
#   - results/overlap/dcca1_ComparaOverlap.txt
#   - results/overlap/dcca1_ComparaNoOverlap.txt
#   - figures/comparison.pdf (Figure S5 in the Supplementary Material)
#
# Note: Pre-computed simulation files are already included in `results/overlap/`.
#       Set `run_simulation <- TRUE` below only if you want to regenerate them.
#       Updated for DCCA >= 1.0.0 (`rhoE()$rhodcca`).
# =============================================================================

library(DCCA)
library(ggplot2)

outdir_res <- "results/overlap"
outdir_fig <- "figures"
dir.create(outdir_res, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_fig, showWarnings = FALSE)

run_simulation <- FALSE

m_min <- 3
m_max <- 25
dm <- m_max - m_min + 1
n <- 100
n_rep <- 1000
m <- seq(m_min, m_max, 1)
thetas <- (20:1) / 10

f_over   <- file.path(outdir_res, "dcca1_ComparaOverlap.txt")
f_noover <- file.path(outdir_res, "dcca1_ComparaNoOverlap.txt")

# =============================================================================
# 1. Monte Carlo simulation (optional if pre-computed files are present)
# =============================================================================
if (run_simulation) {
  unlink(c(f_over, f_noover))
  for (r in 1:n_rep) {
    set.seed(100 + r)
    x <- rnorm(n + length(thetas))
    y <- numeric(n + length(thetas))
    for (i in (length(thetas) + 1):(n + length(thetas))) {
      y[i] <- x[i] + sum(thetas * x[(i - 1):(i - length(thetas))])
    }
    x <- x[-(1:length(thetas))]
    y <- y[-(1:length(thetas))]

    dccacor.o  <- DCCA::rhodcca(x, y, m, nu = 0, overlap = TRUE)$rhodcca
    dccacor.no <- DCCA::rhodcca(x, y, m, nu = 0, overlap = FALSE)$rhodcca

    write.table(t(dccacor.o),  file = f_over,   row.names = FALSE, col.names = FALSE, append = TRUE)
    write.table(t(dccacor.no), file = f_noover, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

# =============================================================================
# 2. Build comparison boxplot (Figure S5)
# =============================================================================
dcca1.o  <- read.table(f_over,   header = FALSE)
dcca1.no <- read.table(f_noover, header = FALSE)

colnames(dcca1.o)  <- m
colnames(dcca1.no) <- m

t1 <- data.frame(stack(dcca1.o),  Group = "Overlapping")
t2 <- data.frame(stack(dcca1.no), Group = "Non-overlapping")
dados <- rbind(t1, t2)
colnames(dados) <- c("rho", "m", "Group")

theta <- c(c(1, (20:1) / 10), rep(0, m_max - 20))
Gamma1  <- diag(m_max + 1)
Gamma2  <- matrix(0, ncol = m_max + 1, nrow = m_max + 1)
Gamma12 <- matrix(0, ncol = m_max + 1, nrow = m_max + 1)
for (t in 1:(m_max + 1)) {
  for (h in 0:(m_max + 1 - t)) {
    Gamma2[t, t + h]  <- sum(theta[1:(length(theta) - h)] * theta[(1 + h):length(theta)])
    Gamma2[t + h, t]  <- Gamma2[t, t + h]
    Gamma12[t, t + h] <- theta[h + 1]
  }
}

rho3 <- DCCA::rhoE(m = m_min:m_max, G1 = Gamma1, G2 = Gamma2, G12 = Gamma12)$rhodcca
rho3 <- data.frame(rho = rho3, m = m)

p <- ggplot(dados, aes(x = as.factor(m), y = rho)) +
  xlab("m") + ylab(expression(rho[DCCA](m))) +
  geom_boxplot(aes(fill = Group)) +
  geom_line(data = rho3, aes(x = as.factor(m), y = rho, group = 1), linewidth = 1.5, col = "blue") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )

pdf(file.path(outdir_fig, "comparison.pdf"), width = 15, height = 5)
print(p)
dev.off()

pdf(file.path(outdir_fig, "comparison2.pdf"), width = 18, height = 6)
print(p)
dev.off()
