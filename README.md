# Codes_DCCA_papers

Reproducible R code, datasets, simulation outputs, and figures for the two theoretical and methodological papers on Detrended Fluctuation Analysis (DFA) and Detrended Cross-Correlation Analysis (DCCA):

1. **Paper 1 (`Paper1/`):** Prass, T. S. and Pumi, G. (2021). *On the behavior of the DFA and DCCA in trend-stationary processes*. Journal of Multivariate Analysis, 182, 104703.
2. **Paper 2 (`Paper2/`):** Prass, T. S. and Pumi, G. (2026). *A note on the behavior of the DFA and DCCA in trend-stationary ARFIMA processes*.

> **Status:** All simulation scripts, pre-computed Monte Carlo outputs, summary tables, and figures for Paper 1 (main text and Supplementary Material) and Paper 2 (Studies 1 and 2) have been reviewed and organized. All scripts are compatible with `DCCA >= 1.0.0`.

Run the scripts in each folder (`Paper1/` or `Paper2/`) with that folder set as the working directory.

---

## Paper 1: Trend-Stationary Short-Memory Processes (`Paper1/`)

Scripts used to produce all figures in the main manuscript and Supplementary Material of Prass and Pumi (2021). Pre-computed Monte Carlo outputs are already included in `Paper1/results/`, so the simulation scripts (`02`, `04`, and the simulation block of `06`) only need to be run if you wish to regenerate the raw outputs from scratch.

| File | Description |
|---|---|
| `01 - Illustrative figures.R` | Generates the conceptual and illustrative figures in Sections 2 and 3 of the main paper: `figures/seriesXandR.pdf` (Figure 1: original series $X_t$ and integrated series $R_t$), `figures/seriesXandRResiduals.pdf` (Figure 2: overlapping vs. non-overlapping boxes, local polynomial fits, and detrended residuals), `figures/DfaDcca.pdf` (Figure 3: empirical vs. theoretical $F^2_{k,\text{DFA}}(m)$, $F_{\text{DCCA}}(m)$, and $\rho_{\text{DCCA}}(m)$), and `figures/fDfaDccam10.pdf` (box-level detrended variances and covariances). |
| `02 - Simulation - Polynomial order and m fixed.R` | Runs two Monte Carlo experiments ($1{,}000$ replications each): (1) effect of polynomial detrending order $\nu \in \{0, 2, 5\}$ for large windows $m \in \{1000, 1010, \dots, 1100\}$ ($n = 2{,}000$) under independent $\text{AR}(1)$ and $\text{MA}(1)$ processes; and (2) consistency as $n$ increases from $500$ to $10{,}000$ for fixed $m = 100$ under two i.i.d. series and under independent $\text{AR}(1)$ and $\text{MA}(1)$ series. Saves raw outputs to `results/part1_nu_mfixed/`. |
| `03 - Figures - Polynomial order and m fixed.R` | Reads `results/part1_nu_mfixed/` and produces: (1) `figures/figNu025.pdf` (Figure 4 in the main paper: boxplots for $\nu \in \{0, 2, 5\}$ and pairwise differences between $\nu = 2$ and $\nu = 5$); and (2) `figures/c1mfixed.pdf` and `figures/c2mfixed.pdf` (Figure S2 in the Supplementary Material: $\rho_{\text{DCCA}}(100)$ vs. $n$). |
| `04 - Simulation - Scenarios 1 to 4.R` | Runs $1{,}000$ Monte Carlo replications ($n = 2{,}000$, $m \in \{3, 5, \dots, 101\}$, $\nu = 0$) across all 10 short-memory configurations of Scenarios 1 to 4 in Section 1 of the Supplementary Material (independent i.i.d. and $\text{AR}(1)/\text{MA}(1)$ processes, signal-plus-noise models, correlated bivariate white noise with $\rho \in \{0.5, 0.8\}$, $\text{MA}(20)$ and $\text{AR}(1)$ cross-correlation structures, and common $\text{AR}(1)$ innovations driving $\text{AR}(1)$ or $\text{MA}(1)$ filters). Saves outputs to `results/part2_scenarios/`. |
| `05 - Figures - Scenarios 1 to 4.R` | Reads `results/part2_scenarios/`, computes exact finite-$m$ theoretical expectations via `DCCA::rhoE()`, and generates all figures in Section 1 of the Supplementary Material: `figures/cn1a.pdf`, `cn1b.pdf` (Figure S1), `cn2cor5.pdf`, `cn2cor8.pdf`, `cn2a.pdf`, `cn2b.pdf` (Figure S3), and `cn3a.pdf`, `cn3b.pdf`, `cn4a.pdf`, `cn4b.pdf` (Figure S4). |
| `06 - Simulation and Figure - Overlap comparison.R` | Compares $\rho_{\text{DCCA}}(m)$ computed with overlapping versus non-overlapping boxes ($n = 100$, $m \in \{3, \dots, 25\}$, $1{,}000$ replications) under the $\text{MA}(20)$ cross-correlation scenario. Reads/writes `results/overlap/` and produces `figures/comparison.pdf` (Figure S5 in the Supplementary Material). |
| `07 - Empirical application.R` | Empirical analysis of Bitcoin (`BTC`) and four stock market indices (`BVSP`, `DJI`, `GSPC`, `IXIC`) from October 18, 2010 to December 30, 2019 ($n = 2{,}217$ common trading days), evaluating log-returns and absolute log-returns across $m \in \{3, 5, \dots, 299\}$. Saves processed series in `results/application/` and produces `figures/scatter.pdf` (Figure S6), `figures/dfas.pdf` (Figure S7), `figures/dccas.pdf` (Figure S8), and `figures/rhodccas.pdf` (Figure S9). |

### Key Helper Functions in `Paper1/`

| Function | Script | Description |
|---|---|---|
| `plot_series()`, `add_bracket()`, `plot_fluct_comp()`, `compute_box_series()`, `plot_box_level()` | `01 - Illustrative figures.R` | Reusable plotting and computation helpers for time series panels, box brackets, empirical vs. theoretical fluctuation curves, and box-level detrended statistics $f_{\text{DFA}}(m, i)$ and $f_{\text{DCCA}}(m, i)$. |
| `read_sim_stack()`, `plot_nu_box()`, `build_stat_row()`, `plot_mfixed()` | `03 - Figures - Polynomial order and m fixed.R` | Loads simulation outputs by polynomial order $\nu \in \{0, 2, 5\}$, builds each 4-panel row of Figure 4 (including the $\nu = 2$ minus $\nu = 5$ difference panel), and plots fixed-$m$ boxplots across sample sizes $n$. |
| `run_curve()`, `save_rep()` | `04 - Simulation - Scenarios 1 to 4.R` | Evaluates $F^2_{1,\text{DFA}}(m)$, $F^2_{2,\text{DFA}}(m)$, $F_{\text{DCCA}}(m)$, and $\rho_{\text{DCCA}}(m)$ across the window grid $m$ for a single Monte Carlo replication and appends the result to the scenario output files. |
| `polyinv.ar()`, `poly.prod()`, `cov_matrix_ma()`, `ccv_matrix_common_innov()`, `ccv_matrix_unidirectional()`, `plot_scenario()` | `05 - Figures - Scenarios 1 to 4.R` | Computes $\text{MA}(\infty)$ weights by polynomial inversion/multiplication, constructs the exact $(m_{\max}+1)\times(m_{\max}+1)$ auto- and cross-covariance matrices $\Gamma_1, \Gamma_2, \Gamma_{12}$, evaluates `DCCA::rhoE()`, and renders each 4-panel scenario figure. |
| `read_close()`, `compute_all_curves()`, `plot_dfa_col()`, `build_5x5_matrix()` | `07 - Empirical application.R` | Reads and aligns the raw price series, computes DFA, DCCA, and $\rho_{\text{DCCA}}$ curves across all 5 series and 10 pairs, and assembles the multi-panel and $5\times 5$ matrix figures. |

### Data and Results in `Paper1/`

| Path | Description |
|---|---|
| `Paper1/data/` | Raw Yahoo Finance daily price series (`BTC-USD.csv`, `BVSP.csv`, `DJI.csv`, `GSPC.csv`, `IXIC.csv`) and dataset description (`Description.txt`). |
| `Paper1/results/part1_nu_mfixed/` | Pre-computed simulation outputs from `02 - Simulation - Polynomial order and m fixed.R`. |
| `Paper1/results/part2_scenarios/` | Pre-computed simulation outputs from `04 - Simulation - Scenarios 1 to 4.R`. |
| `Paper1/results/overlap/` | Pre-computed overlapping vs. non-overlapping simulation outputs from `06 - Simulation and Figure - Overlap comparison.R`. |
| `Paper1/results/application/` | Aligned log-returns (`Returns_rev1.csv`) and empirical DFA, DCCA, and $\rho_{\text{DCCA}}$ curves (`Rev1_2020_*.csv`) generated by `07 - Empirical application.R`. |
| `Paper1/figures/` | All PDF figures for the main paper (Figures 1–4) and Supplementary Material (Figures S1–S9). |

---

## Paper 2: Trend-Stationary Bivariate ARFIMA Processes & Bootstrap Inference (`Paper2/`)

Scripts used to validate the asymptotic theory for bivariate $\text{ARFIMA}(p_k, d_k, q_k)$ processes with mixed memory parameters, evaluate the finite-sample performance of the log-log regression point estimators $(\hat{d}_1, \hat{d}_2, \hat{S}^\epsilon_{11}, \hat{S}^\epsilon_{22}, \hat{S}^\epsilon_{12})$ and parametric bootstrap inference, and run the financial empirical application.

| File | Description |
|---|---|
| `00 - utils_parallel.R` | Cross-platform parallel map helper (`parallel_map()`), using `mclapply` on Unix/macOS and PSOCK clusters (`parLapply`) on Windows. Not meant to be run on its own. |
| `01 - Study 1 - Simulation.R` | **Study 1 simulation and re-estimation stage.** Part 1 simulates $1{,}000$ replications ($n = 2{,}000$, truncation $T = 50{,}000$, $m \in \{10, 30, \dots, 290\}$) for Scenario SC1 ($\text{ARFIMA}(\mathbf{0}, \mathbf{d}, \mathbf{0})$ across all pairs $d_1 \le d_2 \in \{-0.45, -0.25, 0, 0.25, 0.45\}$ and $\rho \in \{-0.5, 0, 0.8\}$) and Scenario SC2 ($\text{ARFIMA}(1, d_k, 1)$ with $d = (0.45, -0.25)$, $\phi_k(z) = 1 + 0.4z$, $\theta_k(z) = 1 + 0.2z$, $\rho = 0.8$), saving `results/result0d0.csv`, `results/result1d1.csv`, and `results/timing_theory.csv`. Part 2 re-estimates $(\hat{d}_1, \hat{d}_2, \hat{S}^\epsilon_{11}, \hat{S}^\epsilon_{22}, \hat{S}^\epsilon_{12})$ from the saved curves and writes `results/study1_full_estimates.csv`. |
| `02 - Study 1 - Tables.R` | Reads `results/study1_full_estimates.csv` and produces: (1) `tables/table_bias_dS_summary.csv` (Table 1 in the paper: pooled mean bias and SD of $\hat{d}$ and $\hat{S}^\epsilon$ by true $d$ and scenario); and (2) `tables/table_percentile_coverage_cells.{csv,tex}` (disaggregated cell-by-cell $95\%$ percentile coverage reference table). |
| `03 - Study 1 - Figures.R` | Reads `results/result0d0.csv`, `results/result1d1.csv`, and `results/study1_full_estimates.csv` and produces: `figures/dfa1.pdf` (Figure 1), `figures/rho08DCCA.pdf` (Figure 2), `figures/SC2.pdf` (Figure 3), `figures/study1_bias_d.pdf` and `figures/study1_bias_S.pdf` (Figure 4), `figures/study1_bias_S12.pdf` (Figure 5), `figures/study1_coverage_heatmap.pdf` (Figure 6), and `figures/rho-05DCCA.pdf` & `figures/rho0DCCA.pdf` (Appendix Figure B.8). |
| `04 - Study 2 - Bootstrap simulation.R` | **Study 2 parametric bootstrap simulation stage.** Runs $R = 200$ Monte Carlo replications with $B = 200$ bootstrap draws ($T = 50{,}000$, $m \in \{10, 30, \dots, 290\}$) across $n \in \{500, 1000, 2000, 5000\}$ for Scenario SC2 ($d = (0.45, -0.25)$, $\rho = 0.8$) and the null Scenario $\text{SC2}_0$ ($d = (0.45, 0)$, $\rho = 0$). Supports block-level checkpointing/resuming. Saves per-block `.rds` files, `_config.rds`, and `_timing_runs.csv` in `results/boot_SC2/`. |
| `05 - Study 2 - Bootstrap tables.R` | Reads `results/boot_SC2/*.rds` and generates all Study 2 summary tables in `tables/`: `table_boot_params.{csv,tex}` (Table 2: parameter bias, RMSE, coverage, and width), `table_boot_rho.tex` and `table_boot_Fdcca.tex` (Table 3: curve-level coverage by $m$ and $n$), `table_boot_tests.{csv,tex}` (Table 4: empirical size and power of $H_0: d_k = 0$ and $H_0: S^\epsilon_{12} = 0$), `table_boot_time.{csv,tex}` (Table 5: computational cost), `table_boot_curves.csv`, and `table_coverage_params_SC2.{csv,tex}`. |
| `06 - Study 2 - Bootstrap figures.R` | Reads `tables/table_boot_curves.csv` and `tables/table_boot_params.csv` and produces `figures/boot_coverage_SC2.pdf` (Figure 7 in the paper: empirical coverage of the $95\%$ bootstrap intervals for $F^2_{1,\text{DFA}}(m)$, $F^2_{2,\text{DFA}}(m)$, $F_{\text{DCCA}}(m)$, and $\rho_{\text{DCCA}}(m)$) and `figures/boot_coverage_params_SC2.pdf`. |
| `07 - Study 2 - Parallel timing test.R` | Diagnostic benchmark timing a single replication of `boot.dS12()` with inner parallelization (`boot.control(parallel = TRUE)`) versus serial execution. |
| `08 - Study 2 - Bias correction test.R` | Exploratory script evaluating whether bootstrap bias correction ($\hat{\theta}_{\text{BC}} = 2\hat{\theta} - \bar{\theta}^*$) reduces the fixed-grid bias of the point estimators in `results/boot_SC2/`. |
| `09 - Application - Functions.R` | Helper functions for `10 - Application - Analysis.R` (`load_data()`, `curve_plot()`, `plot_marginal_panel()`, `plot_pair_matrix()`). Not meant to be run on its own. |
| `10 - Application - Analysis.R` | Joint estimation and two-stage parametric bootstrap inference (`DCCA::estimate()`, `DCCA::bootstrap()`, `update()`) for the 5 financial series (`BTC`, `BVSP`, `DJI`, `GSPC`, `IXIC`, 10 pairs). Saves the full fitted object to `results/application_result.rds`, parameter and pairwise summary tables to `tables/application_{parameters,pairs}.{csv,tex}`, and figures to `figures/application_{dfa_panel,pair_matrix}.pdf`. |

### Functions in `00 - utils_parallel.R` and `09 - Application - Functions.R` (`Paper2/`)

| Function | File | Description |
|---|---|---|
| `parallel_map()` | `00 - utils_parallel.R` | Runs a function over a vector/list in parallel using `mclapply` on Unix/macOS or a PSOCK cluster (`parLapply`, exporting required objects and loading `DCCA`) on Windows, with automatic `tryCatch` wrapping per worker. |
| `load_data()` | `09 - Application - Functions.R` | Reads each raw CSV file, drops missing adjusted close observations, and inner-joins all series by `Date`. |
| `curve_plot()` | `09 - Application - Functions.R` | Plots the bootstrap confidence band (`boot_ci_lower` to `boot_ci_upper`), bootstrap mean (`boot_mean`), and empirical point estimate (`point_estimate`) against window size $m$. |
| `plot_marginal_panel()` | `09 - Application - Functions.R` | Arranges the marginal $F^2_{k,\text{DFA}}(m)$ bootstrap plots across all series into a single horizontal panel. |
| `plot_pair_matrix()` | `09 - Application - Functions.R` | Builds the $p \times p$ pairwise matrix of plots with series names on the diagonal, $F_{\text{DCCA}}(m)$ on the upper triangle, and $\rho_{\text{DCCA}}(m)$ on the lower triangle (marking any pair excluded for Cauchy–Schwarz violation). |

### Data and Results in `Paper2/`

| Path | Description |
|---|---|
| `Paper2/data/` | Raw daily price series (`BTC-USD.csv`, `BVSP.csv`, `DJI.csv`, `GSPC.csv`, `IXIC.csv`) and `Description.txt`. |
| `Paper2/results/result0d0.csv`, `result1d1.csv` | Raw simulated DFA, DCCA, and $\rho_{\text{DCCA}}$ curves for Study 1 (SC1 and SC2, $1{,}000$ replications, $n = 2{,}000$). |
| `Paper2/results/study1_full_estimates.csv` | Re-estimated parameters $(\hat{d}_1, \hat{d}_2, \hat{S}^\epsilon_{11}, \hat{S}^\epsilon_{22}, \hat{S}^\epsilon_{12})$ for all $46{,}000$ series pairs in Study 1. |
| `Paper2/results/timing_theory.csv` | Wall-clock execution times for Study 1. |
| `Paper2/results/boot_SC2/` | Checkpointed `.rds` outputs for all 8 blocks of Study 2 (`SC2` and `SC2_null`, $n \in \{500, 1000, 2000, 5000\}$), reference configuration `_config.rds`, and `_timing_runs.csv`. |
| `Paper2/tables/` | CSV and LaTeX tables produced by `02 - Study 1 - Tables.R`, `05 - Study 2 - Bootstrap tables.R`, and `10 - Application - Analysis.R`. |
| `Paper2/figures/` | PDF figures produced by `03 - Study 1 - Figures.R`, `06 - Study 2 - Bootstrap figures.R`, and `10 - Application - Analysis.R`. |

---

## Requirements

- **R** ($\ge 4.2$)
- **DCCA** ($\ge 1.0.0$)
- **CRAN packages:** `ggplot2`, `ggpubr`, `ggforce`, `gridExtra`, `patchwork`, `dplyr`, `tidyr`, `purrr`, `reshape2`, `lattice`, `MASS`, `xtable`
