## =============================================================================
## true_value.R  --  Monte Carlo truth-value computation (local version)
##
## Place this file in the same directory as setting.csv and data_gen.R,
## then run:
##
##   Rscript true_value.R
##
## Results are written to true_value/true_value_XXX.csv, one file per unique
## truth_id.  All settings (truth_n_pair, truth_n_iter, truth_seed, L, link,
## ...) are read directly from setting.csv -- no command-line arguments needed.
##
## To change paths, edit the two lines under "User-configurable paths" below.
## =============================================================================

suppressPackageStartupMessages({
  library(stats)
})

# ---------------------------------------------------------------------------
# User-configurable paths  (edit these if your layout differs)
# ---------------------------------------------------------------------------
SETTING_FILE <- "setting.csv"   # simulation settings table
OUT_DIR      <- "true_value"    # output folder for true_value_XXX.csv

# ---------------------------------------------------------------------------
# Source data generation helper (assumed to be in the same directory)
# ---------------------------------------------------------------------------
source("data_gen.R")

# ---------------------------------------------------------------------------
# Extract unique truth settings from the full setting table
# ---------------------------------------------------------------------------
keep_cols <- c(
  "truth_id", "beta_case", "alpha", "tau",
  "lambda_D", "lambda_H", "lambda_C",
  "betaD1", "betaD2", "betaH1", "betaH2", "betaC1", "betaC2",
  "L", "link", "truth_n_pair", "truth_n_iter", "truth_seed"
)
settings       <- read.csv(SETTING_FILE, stringsAsFactors = FALSE, check.names = FALSE)
truth_settings <- unique(settings[, keep_cols, drop = FALSE])
truth_settings <- truth_settings[order(truth_settings$truth_id), , drop = FALSE]

# ---------------------------------------------------------------------------
# Fit the truth GLM for one iteration's simulated pairs
# ---------------------------------------------------------------------------
fit_truth_glm <- function(pair_df, link = c("logit", "probit")) {
  link <- match.arg(link)
  use  <- pair_df$resolved == 1L & is.finite(pair_df$z1) & is.finite(pair_df$z2)
  if (!any(use))
    return(list(ok = FALSE, beta = c(NA_real_, NA_real_),
                msg = "No resolved truth pairs."))
  
  fit <- tryCatch(
    glm.fit(
      x        = as.matrix(pair_df[use, c("z1", "z2")]),
      y        = pair_df$y[use],
      weights  = pair_df$resolved[use],
      family   = binomial(link = link),
      intercept = FALSE,
      control  = glm.control(epsilon = 1e-8, maxit = 100)
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error"))
    return(list(ok = FALSE, beta = c(NA_real_, NA_real_),
                msg = conditionMessage(fit)))
  
  beta_hat <- fit$coefficients
  ok <- isTRUE(fit$converged) && length(beta_hat) == 2L && all(is.finite(beta_hat))
  list(ok = ok, beta = as.numeric(beta_hat),
       msg = if (ok) "" else "Truth GLM did not converge cleanly.")
}

# ---------------------------------------------------------------------------
# Read settings and loop over all unique truth_ids sequentially
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Computing truth values for %d unique setting(s). Output -> '%s/'\n",
            nrow(truth_settings), OUT_DIR))

for (i in seq_len(nrow(truth_settings))) {
  row      <- truth_settings[i, , drop = FALSE]
  truth_id <- as.integer(row$truth_id)
  n_pair   <- as.integer(row$truth_n_pair)
  n_iter   <- as.integer(row$truth_n_iter)
  seed0    <- as.integer(row$truth_seed)
  beta_D   <- c(row$betaD1, row$betaD2)
  beta_H   <- c(row$betaH1, row$betaH2)
  
  cat(sprintf(
    "\n[%d/%d] truth_id=%d | link=%s | L=%s | n_pair=%d | n_iter=%d\n",
    i, nrow(truth_settings), truth_id, row$link,
    ifelse(is.infinite(row$L), "Inf", format(row$L)),
    n_pair, n_iter))
  flush.console()
  
  beta_mat <- matrix(NA_real_, nrow = n_iter, ncol = 2L)
  ok_iter  <- rep(FALSE, n_iter)
  
  for (b in seq_len(n_iter)) {
    pair_df <- generate_truth_pairs(
      n_pair   = n_pair,
      alpha    = row$alpha,
      lambda_D = row$lambda_D,
      lambda_H = row$lambda_H,
      beta_D   = beta_D,
      beta_H   = beta_H,
      L        = row$L,
      seed     = seed0 + b
    )
    
    fit          <- fit_truth_glm(pair_df, link = row$link)
    beta_mat[b,] <- fit$beta
    ok_iter[b]   <- isTRUE(fit$ok)
    
    if (b == 1L || b %% 5L == 0L || b == n_iter) {
      cat(sprintf("  iter=%d/%d | ok=%d\n", b, n_iter, as.integer(fit$ok)))
      flush.console()
    }
  }
  
  keep <- ok_iter & is.finite(beta_mat[, 1]) & is.finite(beta_mat[, 2])
  
  out <- data.frame(
    truth_id      = truth_id,     beta_case = row$beta_case,
    alpha         = row$alpha,    tau       = row$tau,
    lambda_D      = row$lambda_D, lambda_H  = row$lambda_H,
    lambda_C      = row$lambda_C,
    betaD1        = row$betaD1,   betaD2    = row$betaD2,
    betaH1        = row$betaH1,   betaH2    = row$betaH2,
    betaC1        = row$betaC1,   betaC2    = row$betaC2,
    L             = row$L,        link      = row$link,
    truth_n_pair  = n_pair,       truth_n_iter = n_iter,
    true_beta1       = if (any(keep)) mean(beta_mat[keep, 1]) else NA_real_,
    true_beta2       = if (any(keep)) mean(beta_mat[keep, 2]) else NA_real_,
    true_beta1_mc_sd = if (any(keep)) stats::sd(beta_mat[keep, 1]) else NA_real_,
    true_beta2_mc_sd = if (any(keep)) stats::sd(beta_mat[keep, 2]) else NA_real_,
    truth_ok_iter    = sum(keep),
    stringsAsFactors = FALSE
  )
  
  outfile <- file.path(OUT_DIR, sprintf("true_value_%03d.csv", truth_id))
  write.csv(out, outfile, row.names = FALSE)
  cat(sprintf("  Saved %s  (ok_iter=%d/%d)\n", outfile, sum(keep), n_iter))
  flush.console()
}

cat(sprintf("\nAll done. '%s/'\n", OUT_DIR))