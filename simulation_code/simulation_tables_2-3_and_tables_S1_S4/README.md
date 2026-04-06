---
output:
  pdf_document: default
  html_document: default
---
# Simulation Code

This folder contains all R code needed to reproduce the simulation studies for Tables 2 and 3 of main text and Table S1 and Table S4 in Web Appendix. The simulation evaluates generalized win fraction regression model estimators for composite endpoints with a survival priority outcome (death, D) and a non-fatal event (hospitalization, H), under covariate-dependent censoring, and compares them to the proportional win-ratio regression estimator (pwreg) and win odds regression estimator.

---

## File Overview

### Scripts you run

| File | Purpose |
|---|---|
| `true_value.R` | **Step 1.** Compute Monte Carlo true parameter values for each scenario |
| `main_sim.R` | **Step 2.** Run the main simulation and fit all estimators |
| `summarize_all.R` | **Step 3.** Aggregate results into a single summary table |

### Settings and inputs

| File | Purpose |
|---|---|
| `setting.csv` | Full simulation design: 128 rows, one per simulation setting |

### Helper scripts in `code/` (sourced automatically — do not run directly)

| File | Purpose |
|---|---|
| `data_gen.R` | Data-generating process: latent event times, censoring, observed data, truth pairs |
| `CreateScoreFun.R` | Score function construction for PIM estimation |
| `censoring.weight.est.R` | Censoring weight estimation (IPCW) |
| `ipcw_weight.m1.R` | IPCW weights for win odds regression model |
| `ipcw.weight.m2.R` | IPCW weights for win fraction regression model |
| `ipcw.weight.mo.R` | IPCW weights for oracle variant |
| `merge.csd.R` | Merge composite survival data structures |
| `model.matrix.pim.csd.R` | Model matrix construction for PIM with composite survival data |
| `pim.csd.fit.R` | Core PIM fitting routine |
| `pim.csd.m1.R` | Win odds regression estimator (unadjusted / IPCW) |
| `pim.csd.m2.R` | Win fraction regression estimator (unadjusted / IPCW) |
| `pim.csd.mo.R` | PIM estimator — oracle variant |
| `win.strategy.R` | Win/loss/tie determination for composite endpoint pairs |

---

## Requirements

**R packages** (install once from CRAN before running):

```r
install.packages(c("survival", "MASS", "nleqslv", "pim", "WR"))
```

The `stats` package is part of base R and requires no installation.

---

## Directory Layout

The expected folder structure is:

```
your-folder/
├── setting.csv
├── true_value.R
├── main_sim.R
├── summarize_all.R
│
├── code/                        ← all helper scripts go here
│   ├── data_gen.R
│   ├── CreateScoreFun.R
│   ├── censoring.weight.est.R
│   ├── ipcw_weight.m1.R
│   ├── ipcw.weight.m2.R
│   ├── ipcw.weight.mo.R
│   ├── merge.csd.R
│   ├── model.matrix.pim.csd.R
│   ├── pim.csd.fit.R
│   ├── pim.csd.m1.R
│   ├── pim.csd.m2.R
│   ├── pim.csd.mo.R
│   └── win.strategy.R
│
├── true_value/                  ← created by true_value.R
│   ├── true_value_001.csv
│   ├── true_value_002.csv
│   └── ...
│
└── results/                     ← created by main_sim.R
    ├── sim_result_001.csv
    ├── sim_result_002.csv
    ├── ...
    └── summary_all.csv          ← created by summarize_all.R
```

The three main scripts (`true_value.R`, `main_sim.R`, `summarize_all.R`) automatically look for helper scripts in the `code/` subfolder. If the helper scripts are placed in the same directory as the main scripts instead, they will be found there as well — no changes needed.

---

## How to Run

The three scripts must be run **in order**. Open a terminal, navigate to the folder containing the scripts, and run:

### Step 1 — Compute true parameter values

```bash
Rscript true_value.R
```

This reads the 32 unique truth scenarios from `setting.csv` and estimates the true regression coefficients (β₁, β₂) for each via Monte Carlo integration (100,000 pairs × 100 iterations per scenario). Results are written to `true_value/true_value_XXX.csv`.

**Runtime note:** Each scenario runs 100 GLM fits on 100,000 pairs. This step is computationally intensive. Consider running overnight or on a subset of scenarios first (see "Running a Subset" below).

### Step 2 — Run the main simulation

```bash
Rscript main_sim.R
```

This loops over all 128 simulation settings in `setting.csv`. For each setting it generates 1,000 replicate datasets and fits three estimators:

- **m1** — Win odds regression (unadjusted or IPCW-adjusted)
- **m2** — Win fraction regression (unadjusted or IPCW-adjusted)
- **pwreg** — Proportional win-ratio regression (from the `WR` package)

Results are written to `results/sim_result_XXX.csv`, one file per `sim_id`.

**Runtime note:** 128 settings × 1,000 replications × 3 estimators is a large job. Expect several hours to days on a personal computer depending on hardware.

### Step 3 — Summarize results

```bash
Rscript summarize_all.R
```

This reads all files from `results/` and `true_value/`, matches each simulation result to its true parameter values, and computes per-setting summary statistics:

- Number of converged replications
- Mean estimate and Monte Carlo SD
- Average standard error (ASE)
- 95% confidence interval coverage probability (CP)
- Relative bias (%)
- For probit-link settings: all of the above after transforming estimates to the logit scale via the delta method

Output is written to `results/summary_all.csv`.

---

## Simulation Design (`setting.csv`)

Each row of `setting.csv` defines one simulation setting. The key columns are:

| Column | Description |
|---|---|
| `sim_id` | Unique ID for this simulation setting (1–128) |
| `truth_id` | Links to the corresponding true-value scenario |
| `beta_case` | Label for the coefficient scenario (e.g. B1, B2) |
| `n` | Sample size per replicate |
| `alpha` | Gumbel–Hougaard copula dependence parameter (α = 1 → independence) |
| `lambda_D`, `lambda_H`, `lambda_C` | Baseline hazard rates for death, hospitalization, censoring |
| `betaD1`, `betaD2` | True log-hazard-ratio coefficients for death |
| `betaH1`, `betaH2` | True log-hazard-ratio coefficients for hospitalization |
| `betaC1`, `betaC2` | Censoring hazard coefficients |
| `L` | Follow-up restriction time (`Inf` = no restriction) |
| `link` | Link function: `logit` or `probit` |
| `adjusted` | Censoring adjustment: `unadjusted` or `IPCW` |
| `n_rep` | Number of simulation replications (1,000) |
| `seed_base` | Base random seed; replicate r uses seed = seed_base + r |
| `truth_n_pair` | Pairs simulated per truth iteration (100,000) |
| `truth_n_iter` | Number of truth iterations averaged over (100) |
| `truth_seed` | Base seed for truth computation |

---

## Running a Subset of Settings

To test the code or run only specific scenarios, edit the loop index at the bottom of `main_sim.R` or `true_value.R`. For example, to run only `sim_id` 1 and 2, change:

```r
for (i in seq_len(nrow(settings))) {
```

to:

```r
for (i in 1:2) {
```

---

## Output Files

### `true_value/true_value_XXX.csv`

One row per truth scenario. Key columns:

| Column | Description |
|---|---|
| `truth_id` | Scenario identifier |
| `true_beta1`, `true_beta2` | Monte Carlo mean of estimated coefficients (the "true" values used to assess bias) |
| `true_beta1_mc_sd`, `true_beta2_mc_sd` | Monte Carlo SD across iterations (measures numerical precision) |
| `truth_ok_iter` | Number of iterations that converged successfully |

### `results/sim_result_XXX.csv`

One row per replicate × estimator combination. Key columns:

| Column | Description |
|---|---|
| `sim_id`, `rep` | Setting and replicate identifiers |
| `method` | Estimator: `pim_m1`, `pim_m2`, or `pwreg_rho0` |
| `ok` | 1 if the fit converged with finite estimates, 0 otherwise |
| `est_beta1`, `est_beta2` | Point estimates |
| `se_beta1`, `se_beta2` | Standard errors |
| `lcl_beta1`, `ucl_beta1` | 95% confidence interval bounds for β₁ |
| `lcl_beta2`, `ucl_beta2` | 95% confidence interval bounds for β₂ |
| `censor_rate_L` | Observed censoring rate before time L in that replicate |

### `results/summary_all.csv`

One row per simulation setting. For each of the three estimators (m1, m2, pwreg) the following columns are reported with the estimator name as a prefix (e.g. `m2_cp_beta1`):

| Suffix | Description |
|---|---|
| `n_ok` | Number of replications with a successful fit |
| `mean_beta1`, `mean_beta2` | Mean estimate across replications |
| `mcsd_beta1`, `mcsd_beta2` | Monte Carlo SD of estimates |
| `ase_beta1`, `ase_beta2` | Mean estimated standard error |
| `cp_beta1`, `cp_beta2` | 95% CI coverage probability |
| `rbpct_beta1`, `rbpct_beta2` | Relative bias (%) |

For probit-link settings, additional columns with prefix `m2_logit_equiv_` and `m1_logit_equiv_` report all of the above after transforming estimates to the logit scale using the delta method (transformation: `qlogis(pnorm(β))`).

---

## Changing Default Paths

Each script has a short block at the top labelled **"User-configurable paths"**. Edit those lines if your files are in different locations:

```r
# In true_value.R
SETTING_FILE <- "setting.csv"
OUT_DIR      <- "true_value"

# In main_sim.R
SETTING_FILE <- "setting.csv"
OUT_DIR      <- "results"
PWREG_RHO    <- 0              # rho parameter for WR::pwreg

# In summarize_all.R
SETTING_FILE <- "setting.csv"
RESULTS_DIR  <- "results"
TRUTH_DIR    <- "true_value"
OUT_FILE     <- file.path(RESULTS_DIR, "summary_all.csv")
```
