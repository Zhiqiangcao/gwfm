# Probability index model for composite survival data
# Revised old-weight implementation. Shared alignment and K_i(t) bugs are fixed.

pim.csd.mo = function(formula, data, SurvTime, Status, covariateC = NULL, priority = 1:2,
                      L = Inf, adjusted = c("unadjusted", "IPCW"), apply.weight = TRUE,
                      start = NULL, control = list(), link = c("logit", "probit", "identity"),
                      model = c("difference", "marginal", "regular", "customized"),
                      compare = if (model == "marginal") "all" else "unique",
                      na.action = getOption("na.action"), keep.data = FALSE, ...){

  model <- match.arg(model)
  link <- match.arg(link)
  adjusted <- match.arg(adjusted)
  if (is.null(na.action)) na.action <- "na.fail"
  if (!is.character(na.action)) na.action <- deparse(substitute(na.action))
  nafun <- match.fun(na.action)

  if (!is.data.frame(data)) data <- as.data.frame(data)

  index <- order(data[, SurvTime])
  data <- data[index, , drop = FALSE]
  rownames(data) <- NULL

  n_ep <- length(priority)
  if (!(model %in% c("difference", "marginal"))) {
    stop("The revised pim.csd.mo() currently supports model='difference' or 'marginal'.")
  }

  .build_pair_design_csd <- function(formula, data, trt_con, model, nafun) {
    mf <- model.frame(formula, data = data, na.action = nafun, drop.unused.levels = TRUE)
    tt <- attr(mf, "terms")
    mm <- model.matrix(tt, data = mf)
    if (model == "difference") {
      x <- mm[trt_con$pid_trt, , drop = FALSE] - mm[trt_con$pid_con, , drop = FALSE]
    } else {
      x <- mm[trt_con$pid_trt, , drop = FALSE]
    }
    list(x = as.matrix(x), penv = list(L = trt_con$pid_trt, R = trt_con$pid_con))
  }

  row_id <- seq_len(nrow(data))
  y_o <- as.matrix(data[, 2:(n_ep + 1), drop = FALSE])
  delta_o <- as.matrix(data[, (n_ep + 2):(1 + 2 * n_ep), drop = FALSE])
  y_n <- pmin(y_o, L)
  delta_n <- ifelse(y_n == L, 1, delta_o)
  datap_n <- as.matrix(cbind(row_id, y_n, delta_n))

  trt_con <- merge.csd(datap = datap_n, model = model, n_ep = n_ep)
  design_obj <- .build_pair_design_csd(formula, data, trt_con, model, nafun)
  x <- design_obj$x
  pair_poset <- design_obj$penv

  win_status <- win.strategy(trt_con = trt_con, priority = priority)

  determined_flag <- NULL
  for (q in seq(3, (2 + 2 * n_ep - 1), 2)) {
    win_qth_res <- win_status[, q:(q + 1), drop = FALSE]
    qth_flag <- rowSums(as.matrix(win_qth_res))
    determined_flag <- cbind(determined_flag, qth_flag)
  }
  ties_flag <- rowSums(determined_flag)
  y <- rowSums(as.matrix(win_status[, seq(3, (2 + 2 * n_ep - 1), 2), drop = FALSE]))

  x <- as.matrix(ties_flag * x)
  p <- ncol(x)
  if (is.null(start)) start <- rep(0, p)

  if (adjusted == "unadjusted") {
    weight_est <- weight_var <- NULL
  } else if (!isTRUE(apply.weight)) {
    warning("adjusted='IPCW' but apply.weight=FALSE, so IPCW weights are not applied.")
    weight_est <- weight_var <- NULL
  } else {
    if (is.null(covariateC) || length(covariateC) == 0L) {
      stop("covariateC must be supplied when adjusted='IPCW'.")
    }
    tdata <- data.frame(datap_n, data[, covariateC, drop = FALSE], check.names = FALSE)
    colnames(tdata) <- c(colnames(data[, c(1:(1 + 2 * n_ep)), drop = FALSE]), covariateC)
    censoring_weight <- censoring.weight.est(
      Data = tdata,
      SurvTime = SurvTime,
      Status = Status,
      covariateC = covariateC
    )
    weight_res <- ipcw.weight.mo(
      tdata = tdata,
      censoring_weight = censoring_weight,
      win_status = win_status,
      priority = priority,
      determined_flag = determined_flag,
      ties_flag = ties_flag,
      L = L,
      n_ep = n_ep
    )
    weight_est <- weight_res$weight_est
    weight_var <- weight_res$weight_var
  }

  est_res <- pim.csd.fit(
    x = x,
    y = y,
    link = link,
    start = start,
    control = control,
    weight_est = weight_est,
    weight_var = weight_var,
    penv = pair_poset,
    ...
  )

  names(est_res$coefficients) <- colnames(x)
  if (!keep.data) {
    x <- matrix(nrow = 0, ncol = 0)
    y <- numeric(0)
  }

  flag <- est_res$flag
  Estimate <- est_res$coefficients
  Std.Error <- sqrt(diag(est_res$vcov))
  if (any(is.na(Std.Error)) || any(is.infinite(Std.Error))) {
    flag <- 0
    Std.Error <- rep(1, p)
  }
  z.value <- Estimate / Std.Error
  p.value <- 2 * pnorm(-abs(z.value))

  coefficients <- data.frame(Estimate, Std.Error, z.value, p.value)
  final_res <- list(
    summary = coefficients,
    coef = est_res$coefficients,
    vcov = est_res$vcov,
    fitted = est_res$fitted,
    flag = flag,
    penv = pair_poset,
    link = link,
    estimators = est_res$estim,
    model.matrix = x,
    response = y,
    keep.data = keep.data,
    model = model
  )
  return(final_res)
}
