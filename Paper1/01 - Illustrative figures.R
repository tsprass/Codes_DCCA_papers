# =============================================================================
# SCRIPT 01: ILLUSTRATIVE FIGURES OF DFA AND DCCA (MAIN PAPER FIGURES 1, 2 & 3)
# =============================================================================
# Paper: Prass, T. S. and Pumi, G. (2021). On the behavior of the DFA and DCCA
#        in trend-stationary processes. Journal of Multivariate Analysis,
#        182, 104703.
#
# Produces the conceptual and illustrative figures in Sections 2 and 3 of the
# main manuscript:
#   - figures/seriesXandR.pdf          (Figure 1: original series X_t and
#                                       integrated series R_t)
#   - figures/seriesXandRResiduals.pdf (Figure 2: overlapping vs.
#                                       non-overlapping boxes, local linear
#                                       fits, and detrended residuals)
#   - figures/DfaDcca.pdf              (Figure 3: simulated bivariate series,
#                                       integrated signals, and empirical vs.
#                                       theoretical DFA, DCCA, and rho_DCCA)
#   - figures/fDfaDccam10.pdf          (Box-level detrended variances and
#                                       covariances f_DFA(m, i) and f_DCCA(m, i)
#                                       for overlapping vs. non-overlapping boxes)
# =============================================================================

library(ggplot2)
library(ggpubr)
library(DCCA)

dir.create("figures", showWarnings = FALSE)

My_theme <- theme_bw() +
  theme(
    text = element_text(size = rel(4), face = "bold"),
    strip.text.x = element_text(size = rel(4)),
    strip.text.y = element_text(size = rel(4))
  )

# =============================================================================
# Helper functions for plotting
# =============================================================================

# Plots a time series (points + line) with optional panel tag (e.g., "(a)")
plot_series <- function(data, y_col, y_lab, tag = NULL, tag_x = NULL, tag_y = NULL) {
  p <- ggplot(data = data, aes(x = .data[["t"]], y = .data[[y_col]])) +
    geom_point() + geom_line() + My_theme +
    xlab("") + ylab(y_lab)
  if (!is.null(tag)) {
    p <- p + annotate("text", x = tag_x, y = tag_y, label = tag, size = 6)
  }
  p
}

# Adds a horizontal bracket with vertical end ticks to a ggplot object
add_bracket <- function(p, x1, x2, y1, y2, col = "black") {
  seg <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
  p +
    geom_segment(data = seg, aes(x = x1, xend = x1, y = y1, yend = y2), colour = col) +
    geom_segment(data = seg, aes(x = x2, xend = x2, y = y1, yend = y2), colour = col) +
    geom_segment(data = seg, aes(x = x1, xend = x2, y = y2, yend = y2), colour = col)
}

# Plots non-overlapping (blue) and overlapping (red) fluctuation curves against
# their asymptotic linear approximation (dark green dashed line)
plot_fluct_comp <- function(dfs, col_no, col_o, y_theo, y_lab, tag, tag_y) {
  ggplot() + My_theme +
    xlab("") + ylab(y_lab) +
    geom_point(data = dfs, aes(x = .data[["m"]], y = .data[[col_no]]), color = "blue") +
    geom_line(data = dfs,  aes(x = .data[["m"]], y = .data[[col_no]]), color = "blue") +
    geom_point(data = dfs, aes(x = .data[["m"]], y = .data[[col_o]]),  color = "red") +
    geom_line(data = dfs,  aes(x = .data[["m"]], y = .data[[col_o]]),  color = "red") +
    geom_line(aes(x = dfs$m, y = y_theo), color = "darkgreen", lwd = 1.05, lty = 2) +
    annotate("text", x = 3, y = tag_y, label = tag, size = 6)
}

# Plots box-level detrended statistics f(m, i) across boxes i with their mean
plot_box_level <- function(data, y_col, col, pch, hline_val, y_lab,
                           seg_x, seg_xend, seg_y, seg_yend, txt_x, txt_y, txt_label) {
  ggplot(data = data) + My_theme +
    geom_line(aes(x = .data[["i"]], y = .data[[y_col]]), col = col) +
    geom_point(aes(x = .data[["i"]], y = .data[[y_col]]), col = col, pch = pch) +
    geom_hline(yintercept = hline_val) +
    xlab("") + ylab(y_lab) +
    annotate("segment", x = seg_x, xend = seg_xend, y = seg_y, yend = seg_yend,
             arrow = arrow(), color = "black") +
    annotate("text", x = txt_x, y = txt_y, label = txt_label, parse = TRUE, size = 4.5)
}

# =============================================================================
# 1. Original time series X_t and integrated process R_t (Figure 1)
# =============================================================================
set.seed(1982)
x  <- rnorm(50)
df <- data.frame(x = x, r = cumsum(x), t = seq_along(x))

gx <- plot_series(df, "x", ~X[t])
gr <- plot_series(df, "r", ~R[t])

pdf("figures/seriesXandR.pdf", width = 20, height = 5)
print(ggarrange(gx, gr, ncol = 2))
dev.off()

# =============================================================================
# 2. Integrated process, local polynomial fits, and residuals (Figure 2)
# =============================================================================
cores <- c("blue", "red")
m <- 9

# --- Non-overlapping boxes ---
wno  <- seq(m + 1, length(x), m + 1)
grno <- ggplot(data = df[1:21, ], aes(x = t, y = r)) + My_theme + xlab("") + ylab(~R[t]) +
  geom_rect(aes(xmin = 0.7, xmax = wno[1] + 0.3, ymin = -2.3, ymax = 2),
            colour = "yellow", fill = "yellow", alpha = 0.01) +
  geom_rect(aes(xmin = wno[1] + 0.7, xmax = wno[2] + 0.3, ymin = -4, ymax = -1.4),
            colour = "darkgreen", fill = "darkgreen", alpha = 0.01)
grno <- add_bracket(grno, 0.7, wno[1] + 0.3, 2.1, 2.2)
grno <- add_bracket(grno, wno[1] + 0.7, wno[2] + 0.3, -1.3, -1.2)
grno <- grno +
  annotate("text", size = 4.5, x = 5, y = 2.45,
           label = "First~Box:~bold(R)[m+1]^{'<1>'}~'='~(R[1]~','~R[2]~', ... ,'~R[m+1])", parse = TRUE) +
  annotate("text", size = 4.5, x = 15.5, y = -0.9,
           label = "Second~Box:~bold(R)[2(m+1)]^{'<m+2>'}~'='~(R[m+2]~','~R[m+3]~', ... ,'~R[2(m+1)])", parse = TRUE) +
  annotate("segment", x = 12, xend = 9, y = 0.5, yend = -0.8, arrow = arrow(), color = "blue") +
  annotate("text", size = 4.5, x = 15.5, y = 0.7,
           label = "Fitted~Values:~tilde(bold(R))[1]~'='~(tilde(R)[1]~','~tilde(R)[2]~', ... ,'~tilde(R)[m+1])", parse = TRUE) +
  annotate("segment", x = 9, xend = 11, y = -3, yend = -3.1, arrow = arrow(), color = "red") +
  annotate("text", size = 4.5, x = 4.5, y = -2.8,
           label = "Fitted~Values:~tilde(bold(R))[2]~'='~(tilde(R)[m+2]~','~tilde(R)[m+3]~', ... ,'~tilde(R)[2(m+1)])", parse = TRUE) +
  geom_point() + geom_line()

for (i in 1:2) {
  grno <- grno + geom_smooth(data = df[(wno[i] - m):wno[i], ], formula = y ~ x,
                             se = FALSE, method = "lm", col = cores[i])
}

# --- Overlapping boxes ---
wo  <- seq(1, length(x) - m, 1)
gro <- ggplot(data = df[1:21, ], aes(x = t, y = r)) + My_theme + xlab("") + ylab(~R[t]) +
  geom_rect(aes(xmin = 1.7, xmax = (m + 2) + 0.3, ymin = -2.7, ymax = 2),
            colour = "darkgreen", fill = "darkgreen", alpha = 0.01) +
  geom_rect(aes(xmin = 0.7, xmax = (m + 1) + 0.3, ymin = -2.3, ymax = 2),
            colour = "yellow", fill = "yellow", alpha = 0.01)
gro <- add_bracket(gro, 0.7, (m + 1) + 0.3, 2.1, 2.2)
gro <- add_bracket(gro, 1.7, (m + 2) + 0.3, -2.8, -2.9)
gro <- gro +
  annotate("text", size = 4.5, x = 5, y = 2.45,
           label = "First~Box:~bold(R)[m+1]^{'<1>'}~'='~(R[1]~','~R[2]~', ... ,'~R[m+1])", parse = TRUE) +
  annotate("text", size = 4.5, x = 6, y = -3.2,
           label = "Second~Box:~bold(R)[m+2]^{'<2>'}~'='~(R[2]~','~R[3]~', ... ,'~R[m+2])", parse = TRUE) +
  annotate("segment", x = 12, xend = 9.5, y = 0.5, yend = -0.9, arrow = arrow(), color = "blue") +
  annotate("text", size = 4.5, x = 16.5, y = 0.7,
           label = "Fitted~Values:~tilde(bold(R))[1]~'='~(tilde(R)[1](1)~','~tilde(R)[2](1)~', ... ,'~tilde(R)[m+1](1))", parse = TRUE) +
  annotate("segment", x = 13, xend = 10.3, y = -0.5, yend = -1.5, arrow = arrow(), color = "red") +
  annotate("text", size = 4.5, x = 16.5, y = -0.3,
           label = "Fitted~Values:~tilde(bold(R))[2]~'='~(tilde(R)[2](2)~','~tilde(R)[3](2)~', ... ,'~tilde(R)[m+2](2))", parse = TRUE) +
  geom_point() + geom_line()

for (i in 1:2) {
  gro <- gro + geom_smooth(data = df[wo[i]:(wo[i] + m), ], formula = y ~ x,
                           se = FALSE, method = "lm", col = cores[i])
}

# --- Detrended residuals (non-overlapping and overlapping) ---
lm1 <- data.frame(eps = residuals(lm(r ~ t, data = df[1:10, ])), t = 1:10)
lm2 <- data.frame(eps = residuals(lm(r ~ t, data = df[11:20, ])), t = 11:20)

grnor <- ggplot(data = lm1, aes(x = t, y = eps)) + My_theme + xlab("") + ylab(~epsilon[t]) +
  geom_rect(aes(xmin = 10.7, xmax = 20.3, ymin = -1.1, ymax = 1.1),
            colour = "darkgreen", fill = "darkgreen", alpha = 0.025) +
  geom_rect(aes(xmin = 0.5, xmax = 10.3, ymin = -1.1, ymax = 1.3),
            colour = "yellow", fill = "yellow", alpha = 0.025) +
  geom_point(color = "blue") + geom_line(color = "blue") +
  geom_point(data = lm2, color = "red") +
  geom_line(data = lm2, color = "red")
grnor <- add_bracket(grnor, 0.5, 10.3, 1.4, 1.5)
grnor <- add_bracket(grnor, 10.7, 20.5, -1.3, -1.4)
grnor <- grnor +
  annotate("text", size = 4.5, x = 5.5, y = 1.7,
           label = "Residuals:~(epsilon[1]~','~epsilon[2]~', ... ,'~epsilon[m+1])", parse = TRUE) +
  annotate("text", size = 4.5, x = 15.5, y = -1.6,
           label = "Residuals:~(epsilon[m+2]~','~epsilon[m+3]~', ... ,'~epsilon[2(m+1)])", parse = TRUE)

lm1o <- data.frame(eps = residuals(lm(r ~ t, data = df[1:10, ])), t = 1:10)
lm2o <- data.frame(eps = residuals(lm(r ~ t, data = df[2:11, ])), t = 2:11)

gror <- ggplot(data = lm1o, aes(x = t, y = eps)) + My_theme + xlab("") + ylab(~epsilon[t]) +
  geom_rect(aes(xmin = 1.7, xmax = 11.3, ymin = -1.4, ymax = 1.05),
            colour = "darkgreen", fill = "darkgreen", alpha = 0.025) +
  geom_rect(aes(xmin = 0.5, xmax = 10.3, ymin = -1.1, ymax = 1.3),
            colour = "yellow", fill = "yellow", alpha = 0.025) +
  geom_point(color = "blue") + geom_line(color = "blue") +
  geom_point(data = lm2o, color = "red") +
  geom_line(data = lm2o, color = "red") + xlim(c(0, 20))
gror <- add_bracket(gror, 0.5, 10.3, 1.4, 1.5)
gror <- add_bracket(gror, 1.7, 11.5, -1.5, -1.6)
gror <- gror +
  annotate("text", size = 4.5, x = 5.5, y = 1.7,
           label = "Residuals:~(epsilon[1](1)~','~epsilon[2](1)~', ... ,'~epsilon[m+1](1))", parse = TRUE) +
  annotate("text", size = 4.5, x = 6.5, y = -1.8,
           label = "Residuals:~(epsilon[2](2)~','~epsilon[3](2)~', ... ,'~epsilon[m+2](2))", parse = TRUE)

pdf("figures/seriesXandRResiduals.pdf", width = 20, height = 10)
print(ggarrange(grno, gro, grnor, gror, ncol = 2, nrow = 2))
dev.off()

# =============================================================================
# 3. DFA, DCCA, and rho_DCCA illustration (Figure 3 and box-level plot)
# =============================================================================
set.seed(1982)
n <- 50
x <- rnorm(n)
alpha <- 2
a <- -2
b <- 2
y <- alpha * x + runif(length(x), a, b)
df <- data.frame(x1 = x, x2 = y, r1 = cumsum(x), r2 = cumsum(y), t = seq_along(x))

gx1 <- plot_series(df, "x1", ~X["1,t"], "(a)", 1, 2.5)
gx2 <- plot_series(df, "x2", ~X["2,t"], "(d)", 1, 4.5)
gr1 <- plot_series(df, "r1", ~R["1,t"], "(b)", 1, 4.0)
gr2 <- plot_series(df, "r2", ~R["2,t"], "(e)", 1, 10.6)

ms <- seq(3, 20, 1)
temp  <- DCCA::rhodcca(df$x1, df$x2, overlap = FALSE, m = ms)
tempo <- DCCA::rhodcca(df$x1, df$x2, overlap = TRUE, m = ms)

dfs <- data.frame(
  dfa1no = temp$F2dfa1,  dfa2no = temp$F2dfa2,
  dccano = temp$Fdcca,   rhono  = temp$rhodcca,
  dfa1o  = tempo$F2dfa1, dfa2o  = tempo$F2dfa2,
  dccao  = tempo$Fdcca,  rhoo   = tempo$rhodcca,
  m = ms
)

var_y <- alpha^2 + (b - a)^2 / 12
dfas1    <- plot_fluct_comp(dfs, "dfa1no", "dfa1o", ms / 15,         ~F["1,DFA"]^2, "(c)", 1.7)
dfas2    <- plot_fluct_comp(dfs, "dfa2no", "dfa2o", ms / 15 * var_y, ~F["2,DFA"]^2, "(f)", 7.1)
dcca     <- plot_fluct_comp(dfs, "dccano", "dccao", ms / 15 * alpha, ~F["DCCA"],    "(g)", 3.1)
rhosdcca <- plot_fluct_comp(dfs, "rhono",  "rhoo",  rep(alpha / sqrt(var_y), length(ms)), ~rho["DCCA"], "(h)", 0.9)

pdf("figures/DfaDcca.pdf", width = 15, height = 9)
print(ggarrange(
  ggarrange(gx1, gr1, dfas1, gx2, gr2, dfas2, ncol = 3, nrow = 2),
  ggarrange(dcca, rhosdcca, ncol = 2, nrow = 1),
  ncol = 1, nrow = 2, heights = c(2, 1)
))
dev.off()

# --- Box-level quantities f_DFA(m, i) and f_DCCA(m, i) at m = 9 ---
compute_box_series <- function(R1, R2, m, overlap = TRUE) {
  Q <- DCCA::Qm(m = m, nu = 0)
  starts <- if (overlap) 1:(length(R1) - m) else ((1:floor(length(R1) / (m + 1))) - 1) * (m + 1) + 1
  res <- vapply(starts, function(s) {
    idx  <- s:(s + m)
    Eps1 <- Q %*% R1[idx]
    Eps2 <- Q %*% R2[idx]
    c(f1 = sum(Eps1^2) / m, f2 = sum(Eps2^2) / m, f12 = sum(Eps1 * Eps2) / m)
  }, numeric(3))
  data.frame(i = seq_along(starts), t(res))
}

dfo   <- compute_box_series(df$r1, df$r2, m = ms[7], overlap = TRUE)
df_no <- compute_box_series(df$r1, df$r2, m = ms[7], overlap = FALSE)

g1o  <- plot_box_level(dfo,   "f1",  "blue", 19, mean(dfo$f1),   ~f["1,DFA"](m, i),  7,   9, 1.50, 0.60,  5.7, 1.6, "{F['1,DFA']^2}(m)")
g2o  <- plot_box_level(dfo,   "f2",  "blue", 17, mean(dfo$f2),   ~f["2,DFA"](m, i), 12,  10, 4.00, 3.10, 13.0, 4.4, "{F['2,DFA']^2}(m)")
g12o <- plot_box_level(dfo,   "f12", "blue", 15, mean(dfo$f12),  ~f[DCCA](m, i),    12,  10, 2.00, 1.10, 13.0, 2.3, "{F[DCCA]^2}(m)")

g1   <- plot_box_level(df_no, "f1",  "red",  19, mean(df_no$f1), ~f["1,DFA"](m, i),  2.9, 2, 0.95, 0.75,  3.5, 1.0, "{F['1,DFA']^2}(m)")
g2   <- plot_box_level(df_no, "f2",  "red",  17, mean(df_no$f2), ~f["2,DFA"](m, i),  3.0, 2, 5.00, 3.50,  3.5, 5.3, "{F['2,DFA']^2}(m)")
g12  <- plot_box_level(df_no, "f12", "red",  15, mean(dfo$f12),  ~f[DCCA](m, i),     3.0, 2, 1.60, 1.10,  3.5, 1.7, "{F[DCCA]^2}(m)")

pdf("figures/fDfaDccam10.pdf", width = 15, height = 9)
print(ggarrange(g1o, g1, g2o, g2, g12o, g12, ncol = 2, nrow = 3, widths = c(3, 1)))
dev.off()
