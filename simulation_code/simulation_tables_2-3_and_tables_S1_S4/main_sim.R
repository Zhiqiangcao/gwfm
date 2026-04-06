## =============================================================================
## main_sim.R  --  Monte Carlo simulation (local version)
##
## Place this file in the same directory as setting.csv and all helper .R
## files, then run:
##
##   Rscript main_sim.R
##
## Results are written to results/sim_result_XXX.csv, one file per sim_id.
## All settings (n, alpha, L, link, n_rep, seed_base, ...) are read directly
## from setting.csv -- no command-line arguments are needed.
##
## To change paths or the pwreg rho parameter, edit the three lines under
## "User-configurable paths" below.
## =============================================================================

suppressPackageStartupMessages({
  library(stats)
  library(survival)
  library(MASS)
  library(nleqslv)
  library(pim)
  library(WR)
})

# ---------------------------------------------------------------------------
# User-configurable paths  (edit these if your layout differs)
# ---------------------------------------------------------------------------
SETTING_FILE <- "setting.csv"   # simulation settings table
OUT_DIR      <- "results"       # output folder for sim_result_XXX.csv
PWREG_RHO    <- 0               # rho parameter passed to WR::pwreg

# ---------------------------------------------------------------------------
# Source all helper scripts (assumed to be in the same directory)
# ---------------------------------------------------------------------------
source("code/CreateScoreFun.R")
source("code/censoring.weight.est.R")
source("code/ipcw.weight.m1.R")
source("code/ipcw.weight.m2.R")
source("code/ipcw.weight.mo.R")
source("code/merge.csd.R")
source("code/model.matrix.pim.csd.R")
source("code/pim.csd.fit.R")
source("code/pim.csd.m1.R")
source("code/pim.csd.m2.R")
source("code/pim.csd.mo.R")
source("code/win.strategy.R")
source("data_gen.R")

# ---------------------------------------------------------------------------
# Import pim internal functions required by the helper scripts
# ---------------------------------------------------------------------------
for (.nm in c("new.pim.env", "new.pim.formula", "poset",
              "has.specials", "has.intercept", "sandwich.vcov")) {
  if (!exists(.nm, mode = "function", inherits = TRUE))
    assign(.nm, getFromNamespace(.nm, "pim"), envir = .GlobalEnv)
}
rm(.nm)

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
safe_diag_se <- function(vcov_mat, p = 2L) {
  out <- rep(NA_real_, p)
  if (is.null(vcov_mat)) return(out)
  vv <- tryCatch(as.matrix(vcov_mat), error = function(e) NULL)
  if (is.null(vv) || nrow(vv) < p || ncol(vv) < p) return(out)
  dd <- diag(vv)[seq_len(p)]
  dd[!is.finite(dd) | dd < 0] <- NA_real_
  sqrt(dd)
}

calc_ci <- function(beta, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  list(
    lower = ifelse(is.finite(beta) & is.finite(se), beta - z * se, NA_real_),
    upper = ifelse(is.finite(beta) & is.finite(se), beta + z * se, NA_real_)
  )
}

# ---------------------------------------------------------------------------
# Data reshaping helpers
# ---------------------------------------------------------------------------
restrict_subject_data <- function(dat, L) {
  y_o     <- as.matrix(dat[, c("y1", "y2")])
  delta_o <- as.matrix(dat[, c("delta1", "delta2")])
  y_n     <- pmin(y_o, L)
  delta_n <- if (is.infinite(L)) delta_o else ifelse(y_n == L, 1L, delta_o)
  data.frame(
    id     = dat$id,
    y1     = y_n[, 1],     y2     = y_n[, 2],
    delta1 = as.integer(delta_n[, 1]),
    delta2 = as.integer(delta_n[, 2]),
    x1     = dat$x1,       x2     = dat$x2,
    stringsAsFactors = FALSE
  )
}

wide_to_pwreg <- function(dat, L) {
  dat_r <- restrict_subject_data(dat, L)
  out   <- vector("list", 2L * nrow(dat_r))
  k     <- 0L
  for (i in seq_len(nrow(dat_r))) {
    if (dat_r$delta2[i] == 1L) {
      k <- k + 1L
      out[[k]] <- data.frame(id = dat_r$id[i], time = dat_r$y2[i],
                             status = 2L,
                             x1 = dat_r$x1[i], x2 = dat_r$x2[i],
                             stringsAsFactors = FALSE)
    }
    k <- k + 1L
    out[[k]] <- data.frame(id = dat_r$id[i], time = dat_r$y1[i],
                           status = as.integer(dat_r$delta1[i]),
                           x1 = dat_r$x1[i], x2 = dat_r$x2[i],
                           stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, out[seq_len(k)])
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Model fitting
# ---------------------------------------------------------------------------
fit_with_pim <- function(dat, row, which = c("m1", "m2")) {
  which   <- match.arg(which)
  pim_fun <- switch(which, m1 = pim.csd.m1, m2 = pim.csd.m2)
  mname   <- paste0("pim_", which)
  
  fit_obj <- tryCatch(
    pim_fun(
      formula      = y1 ~ x1 + x2 - 1,
      data         = dat[, c("id", "y1", "y2", "delta1", "delta2", "x1", "x2")],
      SurvTime     = "y1",
      Status       = "delta1",
      covariateC   = c("x1", "x2"),
      priority     = 1:2,
      L            = row$L,
      adjusted     = row$adjusted,
      apply.weight = TRUE,
      link         = row$link,
      model        = "difference",
      compare      = "unique",
      keep.data    = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(fit_obj, "error"))
    return(list(method = mname, fit_detail = row$adjusted,
                ok = 0L, converged = 0L, flag = NA_integer_,
                beta = c(NA_real_, NA_real_), se = c(NA_real_, NA_real_),
                msg = conditionMessage(fit_obj)))
  
  beta_hat  <- tryCatch(as.numeric(fit_obj$coef), error = function(e) rep(NA_real_, 2L))
  se_hat    <- safe_diag_se(fit_obj$vcov, p = 2L)
  flag      <- if (!is.null(fit_obj$flag)) as.integer(fit_obj$flag) else NA_integer_
  converged <- as.integer(isTRUE(flag == 1L))
  ok        <- as.integer(
    converged == 1L && length(beta_hat) == 2L && length(se_hat) == 2L &&
      all(is.finite(beta_hat)) && all(is.finite(se_hat))
  )
  list(method = mname, fit_detail = row$adjusted,
       ok = ok, converged = converged, flag = flag,
       beta = beta_hat, se = se_hat,
       msg = if (ok == 1L) "" else "PIM did not converge or returned non-finite output.")
}

fit_with_pwreg <- function(dat, row, rho = 0) {
  mname  <- paste0("pwreg_rho", format(rho, trim = TRUE, scientific = FALSE))
  pw_dat <- wide_to_pwreg(dat, row$L)
  
  fit_obj <- tryCatch(
    WR::pwreg(
      ID     = pw_dat$id,
      time   = pw_dat$time,
      status = pw_dat$status,
      Z      = as.matrix(pw_dat[, c("x1", "x2")]),
      rho    = rho
    ),
    error = function(e) e
  )
  
  if (inherits(fit_obj, "error"))
    return(list(method = mname, fit_detail = paste0("rho=", rho),
                ok = 0L, converged = NA_integer_, flag = NA_integer_,
                beta = c(NA_real_, NA_real_), se = c(NA_real_, NA_real_),
                msg = conditionMessage(fit_obj)))
  
  beta_hat  <- tryCatch(as.numeric(fit_obj$beta), error = function(e) rep(NA_real_, 2L))
  se_hat    <- safe_diag_se(fit_obj$Var, p = 2L)
  converged <- if (!is.null(fit_obj$converged)) as.integer(isTRUE(fit_obj$converged)) else NA_integer_
  ok        <- as.integer(
    length(beta_hat) == 2L && length(se_hat) == 2L &&
      all(is.finite(beta_hat)) && all(is.finite(se_hat))
  )
  list(method = mname, fit_detail = paste0("rho=", rho),
       ok = ok, converged = converged, flag = NA_integer_,
       beta = beta_hat, se = se_hat,
       msg = if (ok == 1L) "" else "pwreg returned non-finite output.")
}

make_result_row <- function(fit, row, rep_id, seed, censor_rate_L) {
  ci <- calc_ci(fit$beta, fit$se)
  data.frame(
    sim_id        = row$sim_id,     truth_id   = row$truth_id,
    method        = fit$method,     fit_detail = fit$fit_detail,
    beta_case     = row$beta_case,  n          = row$n,
    alpha         = row$alpha,      tau        = row$tau,
    lambda_D      = row$lambda_D,   lambda_H   = row$lambda_H,
    lambda_C      = row$lambda_C,
    betaD1        = row$betaD1,     betaD2     = row$betaD2,
    betaH1        = row$betaH1,     betaH2     = row$betaH2,
    betaC1        = row$betaC1,     betaC2     = row$betaC2,
    L             = row$L,          link       = row$link,
    adjusted      = row$adjusted,
    rep           = rep_id,         seed       = seed,
    ok            = fit$ok,         converged  = fit$converged,
    flag          = fit$flag,
    est_beta1     = fit$beta[1],    est_beta2  = fit$beta[2],
    se_beta1      = fit$se[1],      se_beta2   = fit$se[2],
    lcl_beta1     = ci$lower[1],    ucl_beta1  = ci$upper[1],
    lcl_beta2     = ci$lower[2],    ucl_beta2  = ci$upper[2],
    censor_rate_L = censor_rate_L,
    msg           = fit$msg,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Read settings and loop over all sim_ids sequentially
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

settings <- read.csv(SETTING_FILE, stringsAsFactors = FALSE, check.names = FALSE)
settings <- settings[order(settings$sim_id), , drop = FALSE]

cat(sprintf("Running %d simulation setting(s). Output -> '%s/'\n",
            nrow(settings), OUT_DIR))

for (i in seq_len(nrow(settings))) {
  row <- settings[i, , drop = FALSE]
  
  cat(sprintf(
    "\n[%d/%d] sim_id=%d | n=%d | alpha=%.3f | L=%s | link=%s | adjusted=%s\n",
    i, nrow(settings), row$sim_id, as.integer(row$n), row$alpha,
    ifelse(is.infinite(row$L), "Inf", format(row$L)),
    row$link, row$adjusted))
  flush.console()
  
  n_rep          <- as.integer(row$n_rep)
  seed_base      <- as.integer(row$seed_base)
  progress_every <- if ("progress_every" %in% names(row)) as.integer(row$progress_every) else 25L
  beta_D         <- c(row$betaD1, row$betaD2)
  beta_H         <- c(row$betaH1, row$betaH2)
  beta_C         <- c(row$betaC1, row$betaC2)
  
  out_list <- vector("list", 3L * n_rep)
  idx      <- 0L
  
  for (r in seq_len(n_rep)) {
    sim_seed <- seed_base + r
    sim_dat  <- generate_observed_subjects(
      n        = as.integer(row$n),
      alpha    = row$alpha,
      lambda_D = row$lambda_D,
      lambda_H = row$lambda_H,
      lambda_C = row$lambda_C,
      beta_D   = beta_D,
      beta_H   = beta_H,
      beta_C   = beta_C,
      seed     = sim_seed,
      return_latent = TRUE
    )
    
    dat_obs       <- sim_dat$observed[, c("id", "y1", "y2", "delta1", "delta2", "x1", "x2")]
    censor_rate_L <- subject_censor_rate_at_L(sim_dat$latent, row$L)
    
    fit_list <- list(
      fit_with_pim(dat_obs, row, which = "m1"),
      fit_with_pim(dat_obs, row, which = "m2"),
      fit_with_pwreg(dat_obs, row, rho = PWREG_RHO)
    )
    
    for (fit in fit_list) {
      idx <- idx + 1L
      out_list[[idx]] <- make_result_row(fit, row, r, sim_seed, censor_rate_L)
    }
    
    if (r == 1L || r %% progress_every == 0L || r == n_rep) {
      cat(sprintf(
        "  rep=%d/%d | censor=%.3f | m1_ok=%d (%.4f,%.4f) | m2_ok=%d (%.4f,%.4f) | pw_ok=%d (%.4f,%.4f)\n",
        r, n_rep, censor_rate_L,
        fit_list[[1]]$ok, fit_list[[1]]$beta[1], fit_list[[1]]$beta[2],
        fit_list[[2]]$ok, fit_list[[2]]$beta[1], fit_list[[2]]$beta[2],
        fit_list[[3]]$ok, fit_list[[3]]$beta[1], fit_list[[3]]$beta[2]
      ))
      flush.console()
    }
  }
  
  out     <- do.call(rbind, out_list[seq_len(idx)])
  outfile <- file.path(OUT_DIR, sprintf("sim_result_%03d.csv", row$sim_id))
  write.csv(out, outfile, row.names = FALSE)
  cat(sprintf("  Saved %s\n", outfile))
  flush.console()
}

cat(sprintf("\nAll done. '\n", OUT_DIR))