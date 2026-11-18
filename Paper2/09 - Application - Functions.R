library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Read and merge the raw price series, one CSV per series
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_data <- function(files) {
  series_list <- imap(files, ~ {
    df <- read.table(.x, sep = ",", header = TRUE, na.strings = "null")
    df <- df[!is.na(df$Adj.Close), ]
    df[, c("Date", "Adj.Close")]
  })

  data <- reduce(series_list, inner_join, by = "Date")
  names(data)[-1] <- names(files)
  data$Date <- as.Date(data$Date)
  data
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Plotting helpers
#
#  Both functions below read directly from the data frames returned by
#  summary() on a "joint_estimates_boot" object (see CodeApl.r): columns
#  `m`, `boot_mean`, `boot_ci_lower`, `boot_ci_upper`, `point_estimate`,
#  exactly as DCCA's own summary.dS_bootstrap()/summary.S12_bootstrap()
#  return them. No custom bootstrap or confidence-interval code is
#  needed here anymore -- estimate()/bootstrap()/update() already do
#  that, including the Cauchy-Schwarz check (boot.S12 excludes an
#  offending pair with a warning) and the zero-in-CI restriction
#  (update()).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

My_theme <- theme_bw() +
  theme(
    text = element_text(size = rel(4), face = "bold"),
    strip.text.x = element_text(size = rel(4)),
    strip.text.y = element_text(size = rel(4)),
    plot.title = element_text(size = rel(4), face = "bold", hjust = 0.5)
  )

# Bootstrap band (mean +/- percentile CI) with the point estimate overlaid,
# for one series/pair/quantity at a time.
curve_plot <- function(df, colore = "grey30", title = NULL) {
  ggplot(df, aes(x = m)) +
    geom_ribbon(aes(ymin = boot_ci_lower, ymax = boot_ci_upper), fill = "grey80", alpha = 0.5) +
    geom_line(aes(y = boot_mean), color = "grey20", linewidth = 1) +
    geom_line(aes(y = point_estimate), color = colore, linetype = "dashed", linewidth = 1) +
    labs(title = title, x = NULL, y = NULL) +
    My_theme
}

# Panel of marginal DFA plots, one per series (the marginal parameters
# d_k, S_k^eps -- and hence F^2_k,DFA(m) -- are bootstrapped once per
# series, not once per pair, so a single panel is all that is needed;
# unlike the pre-package code, there is no "DFA matrix by pair" to draw,
# since that would just repeat the same curve across a row/column).
plot_marginal_panel <- function(marginal_dfa, series_names, ncol = length(series_names)) {
  plots <- lapply(series_names, function(s) {
    curve_plot(marginal_dfa[marginal_dfa$series == s, ], colore = "#9B59B6", title = s)
  })
  arrangeGrob(grobs = plots, ncol = ncol)
}

# Matrix of pairwise plots: diagonal = series name, upper triangle =
# F_DCCA(m), lower triangle = rho_DCCA(m). `pairs_summary` is
# smry$pairs (a named list, one element per pair actually bootstrapped
# -- pairs excluded for violating Cauchy-Schwarz are simply absent and
# shown as a placeholder here, exactly as the pre-package code did).
plot_pair_matrix <- function(pairs_summary, series_names) {
  n <- length(series_names)
  plots <- vector("list", n * n)

  placeholder <- function(label, fontface = "italic", color = "grey50") {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = label, size = 5, color = color, fontface = fontface) +
      theme_void()
  }

  for (i in 1:n) {
    for (j in 1:n) {
      idx <- (i - 1) * n + j
      if (i == j) {
        plots[[idx]] <- placeholder(series_names[i], fontface = "bold", color = "black")
        next
      }
      # pair names follow the order series appear in the data matrix
      # (see estimate()/dfadcca.multi()): the lower-index series first
      pair_name <- paste0(series_names[min(i, j)], "_", series_names[max(i, j)])

      if (!(pair_name %in% names(pairs_summary))) {
        plots[[idx]] <- placeholder("CS violation")
        next
      }
      res <- pairs_summary[[pair_name]]
      plots[[idx]] <- if (i < j) {
        curve_plot(res$Fdcca, colore = "#27AE60")
      } else {
        curve_plot(res$rho, colore = "#F1C40F")
      }
    }
  }

  arrangeGrob(grobs = plots, nrow = n, ncol = n, layout_matrix = matrix(1:(n * n), nrow = n, byrow = TRUE))
}
