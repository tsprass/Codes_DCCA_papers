# =============================================================================
# SCRIPT 07: EMPIRICAL APPLICATION - BITCOIN AND STOCK MARKET INDICES
# =============================================================================
# Paper: Prass, T. S. and Pumi, G. (2021). On the behavior of the DFA and DCCA
#        in trend-stationary processes. Journal of Multivariate Analysis,
#        182, 104703 (Supplementary Material, Section 2).
#
# Analyzes the daily log-returns and absolute log-returns of Bitcoin (BTC) and
# four stock market indices (BVSP, DJI, GSPC, IXIC) from October 18, 2010 to
# December 30, 2019 (n = 2,217 common trading days), over window sizes
# m in {3, 5, ..., 299} with nu = 0.
#
# Input:  `data/BTC-USD.csv`, `data/BVSP.csv`, `data/DJI.csv`, `data/GSPC.csv`,
#         `data/IXIC.csv`
# Output:
#   - results/application/Returns_rev1.csv
#   - results/application/Rev1_2020_{dfas,dfasa,dccas,dccasa,rhodccas,rhodccasa}.csv
#   - figures/scatter.pdf  (Figure S6: time series and scatter plot matrix)
#   - figures/dfas.pdf     (Figure S7: F^2_DFA(m) and log F^2_DFA(m))
#   - figures/dccas.pdf    (Figure S8: F_DCCA(m) and log |F_DCCA(m)|)
#   - figures/rhodccas.pdf (Figure S9: rho_DCCA(m) and log |rho_DCCA(m)|)
# =============================================================================

library(lattice)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(reshape2)
library(DCCA)

datadir    <- "data"
outdir_res <- "results/application"
outdir_fig <- "figures"
dir.create(outdir_res, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_fig, showWarnings = FALSE)

My_theme <- theme_bw() +
  theme(
    text = element_text(size = rel(4), face = "bold"),
    strip.text.x = element_text(size = rel(4)),
    strip.text.y = element_text(size = rel(4)),
    plot.title = element_text(size = rel(4), face = "bold", hjust = 0.5)
  )

# =============================================================================
# 1. Load and align the five price series, compute log-returns
# =============================================================================
series_files <- c(
  BTC  = "BTC-USD.csv",
  BVSP = "BVSP.csv",
  DJI  = "DJI.csv",
  GSPC = "GSPC.csv",
  IXIC = "IXIC.csv"
)
series_names <- names(series_files)

read_close <- function(fname, sname) {
  df <- read.table(file.path(datadir, fname), sep = ",", header = TRUE, na.strings = "null")
  df <- df[!is.na(df$Adj.Close), c("Date", "Adj.Close")]
  colnames(df)[2] <- sname
  df
}

dados <- Reduce(
  function(x, y) merge(x, y, by = "Date"),
  Map(read_close, series_files, series_names)
)
dados$Date <- as.Date(dados$Date)

rets  <- data.frame(Date = dados[-1, "Date"],
                    apply(dados[, -1], 2, function(x) diff(log(x)) * 100))
arets <- abs(rets[, -1])
write.table(rets, file = file.path(outdir_res, "Returns_rev1.csv"), sep = ";", row.names = FALSE)

# =============================================================================
# 2. Time series and pairwise scatter plots (Figure S6: scatter.pdf)
# =============================================================================
teste <- data.frame(Date = rep(rets$Date, length(series_names)), melt(rets[, -1], id.vars = NULL))
colnames(teste) <- c("Date", "serie", "Ret")

gg0 <- ggplot(data = teste) +
  geom_line(aes(x = Date, y = Ret, group = serie), lwd = 0.1) +
  facet_grid(serie ~ ., scale = "free") + ylab("") + xlab("") +
  ggtitle("Log-returns") +
  My_theme

gg1 <- ggplot(rets[, -1], aes(x = .panel_x, y = .panel_y)) +
  geom_point(alpha = 0.2, shape = 16, size = 0.5) +
  geom_autodensity() +
  facet_matrix(vars(everything()), layer.diag = 2) + My_theme +
  theme(axis.text.x = element_text(size = rel(4), angle = 90, vjust = 0.3))

pdf(file = file.path(outdir_fig, "scatter.pdf"), width = 15, height = 8)
print(ggarrange(gg0, gg1, ncol = 2, nrow = 1, widths = c(1, 1.5)))
dev.off()

# =============================================================================
# 3. Compute DFA, DCCA, and rho_DCCA for log-returns and absolute log-returns
# =============================================================================
ms <- seq(3, 300, 2)

# Pair definitions matching the column naming convention of the saved CSVs
pair_specs <- list(
  BVSP_BTC  = c("BVSP", "BTC",  "BTC x BVSP"),
  BVSP_DJI  = c("BVSP", "DJI",  "BVSP x DJI"),
  BVSP_GSPC = c("BVSP", "GSPC", "BVSP x GSPC"),
  BVSP_IXIC = c("BVSP", "IXIC", "BVSP x IXIC"),
  BTC_DJI   = c("BTC",  "DJI",  "BTC x DJI"),
  BTC_GSPC  = c("BTC",  "GSPC", "BTC x GSPC"),
  BTC_IXIC  = c("BTC",  "IXIC", "BTC x IXIC"),
  DJI_GSPC  = c("DJI",  "GSPC", "DJI x GSPC"),
  DJI_IXIC  = c("DJI",  "IXIC", "DJI x IXIC"),
  GSPC_IXIC = c("GSPC", "IXIC", "GSPC x IXIC")
)

compute_all_curves <- function(mat) {
  pair_res <- lapply(pair_specs, function(sp) DCCA::rhodcca(mat[[sp[1]]], mat[[sp[2]]], m = ms))
  df_dfa <- data.frame(
    m    = ms,
    BTC  = pair_res$BVSP_BTC$F2dfa2,
    BVSP = pair_res$BVSP_BTC$F2dfa1,
    DJI  = pair_res$BVSP_DJI$F2dfa2,
    GSPC = pair_res$BVSP_GSPC$F2dfa2,
    IXIC = pair_res$BVSP_IXIC$F2dfa2
  )
  df_dcca <- data.frame(m = ms, lapply(pair_res, `[[`, "Fdcca"))
  df_rho  <- data.frame(m = ms, lapply(pair_res, `[[`, "rhodcca"))
  list(dfa = df_dfa, dcca = df_dcca, rho = df_rho)
}

res_ret  <- compute_all_curves(rets[, -1])
res_aret <- compute_all_curves(arets)

dfas      <- res_ret$dfa;   dfasa     <- res_aret$dfa
dccas     <- res_ret$dcca;  dccasa    <- res_aret$dcca
rhodccas  <- res_ret$rho;   rhodccasa <- res_aret$rho

write.table(dfas,      file.path(outdir_res, "Rev1_2020_dfas.csv"),      col.names = TRUE, sep = ";")
write.table(dfasa,     file.path(outdir_res, "Rev1_2020_dfasa.csv"),     col.names = TRUE, sep = ";")
write.table(dccas,     file.path(outdir_res, "Rev1_2020_dccas.csv"),     col.names = TRUE, sep = ";")
write.table(dccasa,    file.path(outdir_res, "Rev1_2020_dccasa.csv"),    col.names = TRUE, sep = ";")
write.table(rhodccas,  file.path(outdir_res, "Rev1_2020_rhodccas.csv"),  col.names = TRUE, sep = ";")
write.table(rhodccasa, file.path(outdir_res, "Rev1_2020_rhodccasa.csv"), col.names = TRUE, sep = ";")

# =============================================================================
# 4. Generate Figures S7 (dfas.pdf), S8 (dccas.pdf), and S9 (rhodccas.pdf)
# =============================================================================
grf <- function(df1, df2 = NULL, cl, ylb, yll = c(-0.15, 1.05),
                yl = FALSE, fit = FALSE, c1 = "darkblue", c2 = NULL) {
  if (is.null(c2)) c2 <- paste("dark", c1)
  f <- ggplot(data = df1) + My_theme +
    geom_point(aes(x = .data[["m"]], y = .data[[cl]]), col = c2, cex = 0.7) +
    xlab("") + ylab(ylb)
  if (!is.null(df2)) {
    f <- f + geom_point(data = df2, aes(x = .data[["m"]], y = .data[[cl]]), col = c1, cex = 0.7)
  }
  if (fit) {
    f <- f + geom_smooth(data = df1, formula = y ~ x, aes(x = .data[["m"]], y = .data[[cl]]),
                         method = lm, se = FALSE, col = "darkred")
    if (!is.null(df2)) {
      f <- f + geom_smooth(data = df2, formula = y ~ x, aes(x = .data[["m"]], y = .data[[cl]]),
                           method = lm, se = FALSE, col = "darkred")
    }
  }
  if (yl) f <- f + ylim(yll)
  f
}

# --- Figure S7: DFA curves and log-log DFA curves (dfas.pdf) ---
plot_dfa_col <- function(df_in, log_scale = FALSE) {
  lapply(series_names, function(s) {
    p <- grf(df1 = if (log_scale) log(df_in) else df_in, cl = s, ylb = s, c1 = "gray", c2 = "black")
    if (log_scale) {
      offset <- if (s == "BTC") log(var(abs(rets[, "BTC"])) / 15) else 0
      p <- p + geom_smooth(formula = y ~ x, aes(x = log(ms), y = log(ms) + offset),
                           method = "lm", col = "red")
    }
    p
  })
}

make_dfa_column <- function(plot_list, top_title, bottom_title) {
  annotate_figure(
    ggarrange(plotlist = plot_list, ncol = 1, nrow = length(series_names)),
    top    = text_grob(top_title,    face = "bold", size = 14),
    bottom = text_grob(bottom_title, face = "bold", size = 14)
  )
}

g1  <- make_dfa_column(plot_dfa_col(dfas,  FALSE), "Log-returns",          "m")
g2  <- make_dfa_column(plot_dfa_col(dfasa, FALSE), "Absolute log-returns", "m")
g1l <- make_dfa_column(plot_dfa_col(dfas,  TRUE),  "Log-returns",          "log(m)")
g2l <- make_dfa_column(plot_dfa_col(dfasa, TRUE),  "Absolute log-returns", "log(m)")

gga  <- annotate_figure(ggarrange(g1, g2),   top = text_grob(expression({F[DFA]^2}(m)),      face = "bold", size = 18))
ggal <- annotate_figure(ggarrange(g1l, g2l), top = text_grob(expression(log({F[DFA]^2}(m))), face = "bold", size = 18))

pdf(file.path(outdir_fig, "dfas.pdf"), width = 20, height = 15)
print(ggarrange(gga, ggal))
dev.off()

# --- Helper to assemble the 5x5 pairwise matrix for Figures S8 and S9 ---
build_5x5_matrix <- function(lower_plots, upper_plots) {
  diag_labels <- lapply(series_names, function(s) {
    ggplot() + theme_void() + annotate(geom = "text", x = 1, y = 1, label = s, size = 10)
  })
  names(diag_labels) <- series_names

  # Map each (row_i, col_j) cell of the 5x5 grid to the corresponding plot
  pair_key <- function(s1, s2) {
    k1 <- paste0(s1, "_", s2)
    if (k1 %in% names(pair_specs)) k1 else paste0(s2, "_", s1)
  }

  grid_list <- vector("list", 25)
  idx <- 1
  for (i in 1:5) {
    for (j in 1:5) {
      if (i == j) {
        grid_list[[idx]] <- diag_labels[[i]]
      } else if (i < j) {
        grid_list[[idx]] <- upper_plots[[pair_key(series_names[i], series_names[j])]]
      } else {
        grid_list[[idx]] <- lower_plots[[pair_key(series_names[i], series_names[j])]]
      }
      idx <- idx + 1
    }
  }
  ggarrange(plotlist = grid_list, ncol = 5, nrow = 5)
}

# --- Figure S8: DCCA curves and log-log DCCA curves (dccas.pdf) ---
dcca_lower <- lapply(names(pair_specs), function(pk) {
  grf(df1 = dccas, df2 = dccasa, cl = pk, ylb = pair_specs[[pk]][3], c1 = "blue", c2 = "darkgreen")
})
names(dcca_lower) <- names(pair_specs)

ld1 <- -log(abs(dccas));  ld1$m <- -ld1$m
ld2 <- -log(abs(dccasa)); ld2$m <- -ld2$m

dcca_upper <- lapply(names(pair_specs), function(pk) {
  lbl <- pair_specs[[pk]][3]
  if (pk == "BVSP_BTC") {
    grf(df1 = ld1, df2 = ld2, cl = pk, ylb = lbl, c1 = "blue", c2 = "darkgreen") +
      geom_smooth(formula = y ~ x, aes(x = log(ms), y = -log(ms)), method = "lm", col = "red")
  } else {
    grf(df1 = log(abs(dccas)), df2 = log(abs(dccasa)), cl = pk, ylb = lbl, c1 = "blue", c2 = "darkgreen") +
      geom_smooth(formula = y ~ x, aes(x = log(ms), y = log(ms)), method = "lm", col = "red")
  }
})
names(dcca_upper) <- names(pair_specs)

pdf(file.path(outdir_fig, "dccas.pdf"), width = 20, height = 15)
print(build_5x5_matrix(dcca_lower, dcca_upper))
dev.off()

# --- Figure S9: rho_DCCA curves and log-log rho_DCCA curves (rhodccas.pdf) ---
rho_lower <- lapply(names(pair_specs), function(pk) {
  grf(df1 = rhodccas, df2 = rhodccasa, cl = pk, ylb = pair_specs[[pk]][3],
      c1 = "blue", c2 = "darkgreen", yl = FALSE)
})
names(rho_lower) <- names(pair_specs)

ld1r <- -log(abs(rhodccas));  ld1r$m <- -ld1r$m
ld2r <- -log(abs(rhodccasa)); ld2r$m <- -ld2r$m

rho_upper <- lapply(names(pair_specs), function(pk) {
  lbl <- pair_specs[[pk]][3]
  if (pk == "BVSP_BTC") {
    grf(df1 = ld1r, df2 = ld2r, cl = pk, ylb = lbl, c1 = "blue", c2 = "darkgreen")
  } else {
    grf(df1 = log(abs(rhodccas)), df2 = log(abs(rhodccasa)), cl = pk, ylb = lbl,
        c1 = "blue", c2 = "darkgreen")
  }
})
names(rho_upper) <- names(pair_specs)

pdf(file.path(outdir_fig, "rhodccas.pdf"), width = 20, height = 15)
print(build_5x5_matrix(rho_lower, rho_upper))
dev.off()
