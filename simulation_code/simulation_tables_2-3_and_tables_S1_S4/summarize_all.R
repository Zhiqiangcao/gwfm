suppressPackageStartupMessages({
  library(stats)
})

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

safe_rel_bias_pct <- function(est, truth) {
  ifelse(is.finite(est) & is.finite(truth) & abs(truth) > 1e-8,
         100 * (est - truth) / truth,
         NA_real_)
}

clamp_prob <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

probit_to_logit_coef <- function(beta) {
  stats::qlogis(clamp_prob(stats::pnorm(beta)))
}

probit_to_logit_deriv <- function(beta) {
  p <- clamp_prob(stats::pnorm(beta))
  stats::dnorm(beta) / (p * (1 - p))
}

empty_block <- function() {
  list(
    n_ok = 0L,
    mean_beta1 = NA_real_,
    mean_beta2 = NA_real_,
    mcsd_beta1 = NA_real_,
    mcsd_beta2 = NA_real_,
    ase_beta1 = NA_real_,
    ase_beta2 = NA_real_,
    cp_beta1 = NA_real_,
    cp_beta2 = NA_real_,
    rbpct_beta1 = NA_real_,
    rbpct_beta2 = NA_real_
  )
}

prefix_block <- function(x, prefix) {
  names(x) <- paste0(prefix, names(x))
  x
}

normalize_method <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x0))
  
  out[x0 %in% c("m2", "pim_m2", "pim.csd.m2", "pim_csd_m2")] <- "m2"
  out[x0 %in% c("m1", "pim_m1", "pim.csd.m1", "pim_csd_m1")] <- "m1"
  out[grepl("^pwreg($|_)", x0) | x0 %in% c("pwreg", "mao", "wr", "wr_pwreg")] <- "pwreg"
  
  out
}

infer_method_from_path <- function(path) {
  x <- tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
  if (grepl("(^|/|_|-|\\.)m2($|/|_|-|\\.)|pim.?m2", x)) return("m2")
  if (grepl("(^|/|_|-|\\.)m1($|/|_|-|\\.)|pim.?m1", x)) return("m1")
  if (grepl("pwreg|mao|wr", x)) return("pwreg")
  NA_character_
}

is_truth_dir <- function(path) {
  dir.exists(path) && length(list.files(path, pattern = "^true_value_.*\\.csv$")) > 0L
}

is_result_dir <- function(path) {
  dir.exists(path) && length(list.files(path, pattern = "^sim_result_.*\\.csv$")) > 0L
}

standardize_sim_df <- function(df, method_label = NA_character_) {
  if (nrow(df) == 0L) return(df)
  
  if (!"rep" %in% names(df)) {
    if ("iter" %in% names(df)) {
      df$rep <- df$iter
    } else {
      df$rep <- seq_len(nrow(df))
    }
  }
  
  if (!"ok" %in% names(df)) {
    if ("converged" %in% names(df)) {
      df$ok <- as.integer(df$converged %in% c(1, TRUE, "1", "TRUE", "True"))
    } else {
      df$ok <- 1L
    }
  }
  
  needed_num <- c("est_beta1", "est_beta2", "se_beta1", "se_beta2", "censor_rate_L")
  for (nm in needed_num) {
    if (!nm %in% names(df)) df[[nm]] <- NA_real_
  }
  
  z975 <- stats::qnorm(0.975)
  
  if (!"lcl_beta1" %in% names(df)) df$lcl_beta1 <- ifelse(is.finite(df$est_beta1) & is.finite(df$se_beta1), df$est_beta1 - z975 * df$se_beta1, NA_real_)
  if (!"ucl_beta1" %in% names(df)) df$ucl_beta1 <- ifelse(is.finite(df$est_beta1) & is.finite(df$se_beta1), df$est_beta1 + z975 * df$se_beta1, NA_real_)
  if (!"lcl_beta2" %in% names(df)) df$lcl_beta2 <- ifelse(is.finite(df$est_beta2) & is.finite(df$se_beta2), df$est_beta2 - z975 * df$se_beta2, NA_real_)
  if (!"ucl_beta2" %in% names(df)) df$ucl_beta2 <- ifelse(is.finite(df$est_beta2) & is.finite(df$se_beta2), df$est_beta2 + z975 * df$se_beta2, NA_real_)
  
  if (!"method" %in% names(df)) {
    df$method <- method_label
  } else {
    bad <- is.na(df$method) | !nzchar(trimws(as.character(df$method)))
    if (!is.na(method_label) && any(bad)) df$method[bad] <- method_label
  }
  
  df$ok <- as.integer(df$ok %in% c(1, TRUE, "1", "TRUE", "True"))
  df
}

read_sim_dir <- function(result_dir, method_label = NA_character_) {
  if (is.null(result_dir) || is.na(result_dir) || !nzchar(result_dir) || !dir.exists(result_dir)) {
    return(data.frame())
  }
  files <- list.files(result_dir, pattern = "^sim_result_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0L) return(data.frame())
  out <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  standardize_sim_df(out, method_label = method_label)
}

read_truth_df <- function(truth_dir) {
  truth_files <- list.files(truth_dir, pattern = "^true_value_.*\\.csv$", full.names = TRUE)
  if (length(truth_files) == 0L) stop("No truth-value files found in ", truth_dir)
  do.call(rbind, lapply(truth_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
}

summarize_method_block <- function(sim_df, truth_beta) {
  out <- empty_block()
  if (is.null(sim_df) || nrow(sim_df) == 0L) return(out)
  
  keep <- sim_df$ok == 1L &
    is.finite(sim_df$est_beta1) & is.finite(sim_df$est_beta2) &
    is.finite(sim_df$se_beta1) & is.finite(sim_df$se_beta2) &
    is.finite(sim_df$lcl_beta1) & is.finite(sim_df$ucl_beta1) &
    is.finite(sim_df$lcl_beta2) & is.finite(sim_df$ucl_beta2)
  
  if (!any(keep)) return(out)
  
  est1 <- sim_df$est_beta1[keep]
  est2 <- sim_df$est_beta2[keep]
  se1 <- sim_df$se_beta1[keep]
  se2 <- sim_df$se_beta2[keep]
  lcl1 <- sim_df$lcl_beta1[keep]
  ucl1 <- sim_df$ucl_beta1[keep]
  lcl2 <- sim_df$lcl_beta2[keep]
  ucl2 <- sim_df$ucl_beta2[keep]
  
  t1 <- truth_beta[1]
  t2 <- truth_beta[2]
  
  out$n_ok <- sum(keep)
  out$mean_beta1 <- mean_or_na(est1)
  out$mean_beta2 <- mean_or_na(est2)
  out$mcsd_beta1 <- sd_or_na(est1)
  out$mcsd_beta2 <- sd_or_na(est2)
  out$ase_beta1 <- mean_or_na(se1)
  out$ase_beta2 <- mean_or_na(se2)
  out$cp_beta1 <- mean(lcl1 <= t1 & t1 <= ucl1, na.rm = TRUE)
  out$cp_beta2 <- mean(lcl2 <= t2 & t2 <= ucl2, na.rm = TRUE)
  out$rbpct_beta1 <- safe_rel_bias_pct(out$mean_beta1, t1)
  out$rbpct_beta2 <- safe_rel_bias_pct(out$mean_beta2, t2)
  out
}

summarize_logit_equiv_block <- function(sim_df, link_value, logit_truth_beta) {
  out <- empty_block()
  if (!identical(as.character(link_value), "probit")) return(out)
  if (is.null(sim_df) || nrow(sim_df) == 0L) return(out)
  
  keep <- sim_df$ok == 1L &
    is.finite(sim_df$est_beta1) & is.finite(sim_df$est_beta2) &
    is.finite(sim_df$se_beta1) & is.finite(sim_df$se_beta2) &
    is.finite(sim_df$lcl_beta1) & is.finite(sim_df$ucl_beta1) &
    is.finite(sim_df$lcl_beta2) & is.finite(sim_df$ucl_beta2)
  
  if (!any(keep)) return(out)
  
  est1 <- sim_df$est_beta1[keep]
  est2 <- sim_df$est_beta2[keep]
  se1 <- sim_df$se_beta1[keep]
  se2 <- sim_df$se_beta2[keep]
  lcl1 <- sim_df$lcl_beta1[keep]
  ucl1 <- sim_df$ucl_beta1[keep]
  lcl2 <- sim_df$lcl_beta2[keep]
  ucl2 <- sim_df$ucl_beta2[keep]
  
  tr_est1 <- probit_to_logit_coef(est1)
  tr_est2 <- probit_to_logit_coef(est2)
  tr_se1 <- abs(probit_to_logit_deriv(est1)) * se1
  tr_se2 <- abs(probit_to_logit_deriv(est2)) * se2
  
  tr_lcl1 <- probit_to_logit_coef(lcl1)
  tr_ucl1 <- probit_to_logit_coef(ucl1)
  tr_lcl2 <- probit_to_logit_coef(lcl2)
  tr_ucl2 <- probit_to_logit_coef(ucl2)
  
  low1 <- pmin(tr_lcl1, tr_ucl1)
  high1 <- pmax(tr_lcl1, tr_ucl1)
  low2 <- pmin(tr_lcl2, tr_ucl2)
  high2 <- pmax(tr_lcl2, tr_ucl2)
  
  t1 <- logit_truth_beta[1]
  t2 <- logit_truth_beta[2]
  
  out$n_ok <- sum(keep)
  out$mean_beta1 <- mean_or_na(tr_est1)
  out$mean_beta2 <- mean_or_na(tr_est2)
  out$mcsd_beta1 <- sd_or_na(tr_est1)
  out$mcsd_beta2 <- sd_or_na(tr_est2)
  out$ase_beta1 <- mean_or_na(tr_se1)
  out$ase_beta2 <- mean_or_na(tr_se2)
  out$cp_beta1 <- mean(low1 <= t1 & t1 <= high1, na.rm = TRUE)
  out$cp_beta2 <- mean(low2 <= t2 & t2 <= high2, na.rm = TRUE)
  out$rbpct_beta1 <- safe_rel_bias_pct(out$mean_beta1, t1)
  out$rbpct_beta2 <- safe_rel_bias_pct(out$mean_beta2, t2)
  out
}

mean_censor_rate_from_all_methods <- function(sim_df) {
  if (is.null(sim_df) || nrow(sim_df) == 0L || !all(c("rep", "censor_rate_L") %in% names(sim_df))) {
    return(NA_real_)
  }
  tmp <- sim_df[is.finite(sim_df$censor_rate_L), c("rep", "censor_rate_L"), drop = FALSE]
  if (nrow(tmp) == 0L) return(NA_real_)
  tmp <- tmp[!duplicated(tmp$rep), , drop = FALSE]
  mean_or_na(tmp$censor_rate_L)
}

match_truth_current <- function(setting_row, truth_df) {
  hit <- truth_df[
    truth_df$truth_id == setting_row$truth_id &
      truth_df$link == setting_row$link &
      truth_df$L == setting_row$L,
    , drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop("Could not uniquely match current-link truth row for sim_id=", setting_row$sim_id)
  }
  hit
}

match_truth_logit_target <- function(setting_row, truth_df) {
  key_cols <- c(
    "beta_case", "alpha", "tau", "lambda_D", "lambda_H", "lambda_C",
    "betaD1", "betaD2", "betaH1", "betaH2", "betaC1", "betaC2", "L"
  )
  hit <- truth_df[truth_df$link == "logit", , drop = FALSE]
  for (nm in key_cols) {
    hit <- hit[hit[[nm]] == setting_row[[nm]], , drop = FALSE]
  }
  if (nrow(hit) != 1L) {
    stop("Could not uniquely match logit-target truth row for sim_id=", setting_row$sim_id)
  }
  hit
}

load_results_flexible <- function(args) {
  if (length(args) >= 6L && is_truth_dir(args[2]) && is_result_dir(args[4])) {
    mode <- "separate"
    setting_file <- args[1]
    truth_dir <- args[2]
    outfile <- args[3]
    m2_dir <- args[4]
    m1_dir <- args[5]
    pwreg_dir <- args[6]
    
    sim_df <- rbind(
      read_sim_dir(m2_dir, method_label = "m2"),
      read_sim_dir(m1_dir, method_label = "m1"),
      read_sim_dir(pwreg_dir, method_label = "pwreg")
    )
  } else {
    mode <- "combined"
    setting_file <- if (length(args) >= 1L) args[1] else "setting.csv"
    result_dir <- if (length(args) >= 2L) args[2] else "results"
    truth_dir <- if (length(args) >= 3L) args[3] else "true_value"
    outfile <- if (length(args) >= 4L) args[4] else file.path(result_dir, "summary_all.csv")
    
    sim_df <- read_sim_dir(result_dir)
    if (nrow(sim_df) == 0L) {
      stop("No simulation-result files found in ", result_dir)
    }
    
    if (!"method" %in% names(sim_df) || all(is.na(sim_df$method) | !nzchar(trimws(as.character(sim_df$method))))) {
      inferred <- infer_method_from_path(result_dir)
      if (is.na(inferred)) {
        stop(
          "Result files do not contain a usable 'method' column. ",
          "Use either a combined results directory with method labels, or run as: ",
          "Rscript summarize_all_clean_v2.R setting.csv true_value summary.csv results_m2 results_m1 results_pwreg"
        )
      }
      sim_df$method <- inferred
    }
  }
  
  sim_df <- standardize_sim_df(sim_df)
  sim_df$method_group <- normalize_method(sim_df$method)
  
  unknown_methods <- unique(sim_df$method[is.na(sim_df$method_group)])
  if (length(unknown_methods) > 0L) {
    warning("Unrecognized method labels were ignored: ", paste(unknown_methods, collapse = ", "))
  }
  
  list(
    mode = mode,
    setting_file = setting_file,
    truth_dir = truth_dir,
    outfile = outfile,
    sim_df = sim_df
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- load_results_flexible(args)
  
  settings <- read.csv(cfg$setting_file, stringsAsFactors = FALSE, check.names = FALSE)
  settings <- settings[order(settings$sim_id), , drop = FALSE]
  truth_df <- read_truth_df(cfg$truth_dir)
  sim_df <- cfg$sim_df
  
  out_list <- vector("list", nrow(settings))
  
  for (i in seq_len(nrow(settings))) {
    srow <- settings[i, , drop = FALSE]
    sim_cur <- sim_df[sim_df$sim_id == srow$sim_id, , drop = FALSE]
    
    truth_cur <- match_truth_current(srow, truth_df)
    truth_logit <- match_truth_logit_target(srow, truth_df)
    
    sim_m2 <- sim_cur[sim_cur$method_group == "m2", , drop = FALSE]
    sim_m1 <- sim_cur[sim_cur$method_group == "m1", , drop = FALSE]
    sim_pw <- sim_cur[sim_cur$method_group == "pwreg", , drop = FALSE]
    
    base_row <- data.frame(
      sim_id = srow$sim_id,
      truth_id = srow$truth_id,
      beta_case = srow$beta_case,
      n = srow$n,
      alpha = srow$alpha,
      tau = if ("tau" %in% names(srow)) srow$tau else NA_real_,
      adjusted = srow$adjusted,
      link = srow$link,
      betaD1 = srow$betaD1,
      betaD2 = srow$betaD2,
      betaH1 = srow$betaH1,
      betaH2 = srow$betaH2,
      betaC1 = srow$betaC1,
      betaC2 = srow$betaC2,
      lambda_D = srow$lambda_D,
      lambda_H = srow$lambda_H,
      lambda_C = srow$lambda_C,
      L = srow$L,
      n_rep = srow$n_rep,
      mean_censor_rate_L = mean_censor_rate_from_all_methods(sim_cur),
      true_beta1 = truth_cur$true_beta1,
      true_beta2 = truth_cur$true_beta2,
      logit_target_true_beta1 = if (identical(as.character(srow$link), "probit")) truth_logit$true_beta1 else NA_real_,
      logit_target_true_beta2 = if (identical(as.character(srow$link), "probit")) truth_logit$true_beta2 else NA_real_,
      probit_to_logit_transform = if (identical(as.character(srow$link), "probit")) "qlogis(pnorm(beta))" else NA_character_,
      stringsAsFactors = FALSE
    )
    
    blk_m2 <- summarize_method_block(sim_m2, c(truth_cur$true_beta1, truth_cur$true_beta2))
    blk_m1 <- summarize_method_block(sim_m1, c(truth_cur$true_beta1, truth_cur$true_beta2))
    blk_pw <- summarize_method_block(sim_pw, c(truth_cur$true_beta1, truth_cur$true_beta2))
    
    blk_m2_logit_equiv <- summarize_logit_equiv_block(sim_m2, srow$link, c(truth_logit$true_beta1, truth_logit$true_beta2))
    blk_m1_logit_equiv <- summarize_logit_equiv_block(sim_m1, srow$link, c(truth_logit$true_beta1, truth_logit$true_beta2))
    
    out_list[[i]] <- as.data.frame(c(
      as.list(base_row[1, , drop = FALSE]),
      prefix_block(blk_m2, "m2_"),
      prefix_block(blk_m1, "m1_"),
      prefix_block(blk_pw, "pwreg_"),
      prefix_block(blk_m2_logit_equiv, "m2_logit_equiv_"),
      prefix_block(blk_m1_logit_equiv, "m1_logit_equiv_")
    ), stringsAsFactors = FALSE)
  }
  
  out <- do.call(rbind, out_list)
  
  desired_order <- c(
    "sim_id", "truth_id", "beta_case", "n", "alpha", "tau", "adjusted", "link",
    "betaD1", "betaD2", "betaH1", "betaH2", "betaC1", "betaC2",
    "lambda_D", "lambda_H", "lambda_C", "L", "n_rep", "mean_censor_rate_L",
    "true_beta1", "true_beta2", "logit_target_true_beta1", "logit_target_true_beta2", "probit_to_logit_transform",
    "m2_n_ok", "m2_mean_beta1", "m2_mean_beta2", "m2_mcsd_beta1", "m2_mcsd_beta2", "m2_ase_beta1", "m2_ase_beta2", "m2_cp_beta1", "m2_cp_beta2", "m2_rbpct_beta1", "m2_rbpct_beta2",
    "m1_n_ok", "m1_mean_beta1", "m1_mean_beta2", "m1_mcsd_beta1", "m1_mcsd_beta2", "m1_ase_beta1", "m1_ase_beta2", "m1_cp_beta1", "m1_cp_beta2", "m1_rbpct_beta1", "m1_rbpct_beta2",
    "pwreg_n_ok", "pwreg_mean_beta1", "pwreg_mean_beta2", "pwreg_mcsd_beta1", "pwreg_mcsd_beta2", "pwreg_ase_beta1", "pwreg_ase_beta2", "pwreg_cp_beta1", "pwreg_cp_beta2", "pwreg_rbpct_beta1", "pwreg_rbpct_beta2",
    "m2_logit_equiv_n_ok", "m2_logit_equiv_mean_beta1", "m2_logit_equiv_mean_beta2", "m2_logit_equiv_mcsd_beta1", "m2_logit_equiv_mcsd_beta2", "m2_logit_equiv_ase_beta1", "m2_logit_equiv_ase_beta2", "m2_logit_equiv_cp_beta1", "m2_logit_equiv_cp_beta2", "m2_logit_equiv_rbpct_beta1", "m2_logit_equiv_rbpct_beta2",
    "m1_logit_equiv_n_ok", "m1_logit_equiv_mean_beta1", "m1_logit_equiv_mean_beta2", "m1_logit_equiv_mcsd_beta1", "m1_logit_equiv_mcsd_beta2", "m1_logit_equiv_ase_beta1", "m1_logit_equiv_ase_beta2", "m1_logit_equiv_cp_beta1", "m1_logit_equiv_cp_beta2", "m1_logit_equiv_rbpct_beta1", "m1_logit_equiv_rbpct_beta2"
  )
  
  desired_order <- desired_order[desired_order %in% names(out)]
  out <- out[, desired_order, drop = FALSE]
  
  write.csv(out, cfg$outfile, row.names = FALSE)
  cat(sprintf("Saved %s\n", cfg$outfile))
}

main()
