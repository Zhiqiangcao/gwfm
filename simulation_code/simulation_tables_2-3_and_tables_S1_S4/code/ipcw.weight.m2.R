ipcw.weight.m2 = function(tdata, trt_con, censoring_weight, win_status,
                          priority, determined_flag, L, n_ep){
  n <- nrow(tdata)
  y1 <- tdata[, 2]
  npim <- nrow(win_status)

  Kij_ext <- if (!is.null(censoring_weight$Kij_ext)) {
    censoring_weight$Kij_ext
  } else {
    cbind(rep(1, n), censoring_weight$Kij)
  }
  time_grid <- if (!is.null(censoring_weight$time_grid)) {
    censoring_weight$time_grid
  } else {
    c(0, y1)
  }

  Omega <- censoring_weight$Omega
  Omega_inv <- MASS::ginv(Omega)
  U <- censoring_weight$U
  mU <- colMeans(U)
  pc <- length(mU)

  if (!is.null(censoring_weight$Gij)) {
    Gij <- censoring_weight$Gij
  } else {
    G <- censoring_weight$G
    Gij <- array(0, dim = c(n, n + 1L, ncol(G)))
    for (i in seq_len(n)) {
      Gij[i, i + 1L, ] <- G[i, ]
    }
  }

  con_id <- win_status$pid_con
  trt_id <- win_status$pid_trt

  step_index <- function(tt) {
    findInterval(tt, time_grid, rightmost.closed = TRUE)
  }

  extract_surv <- function(id, idx) {
    Kij_ext[cbind(id, idx)]
  }

  extract_G <- function(id, idx) {
    out <- matrix(0, nrow = length(id), ncol = pc)
    for (k in seq_len(pc)) {
      out[, k] <- Gij[cbind(id, idx, k)]
    }
    out
  }

  IPCW_W <- epsilon_W <- matrix(0, npim, n_ep)
  Ieqif <- rep(1, npim)

  for (q in seq_len(n_ep)) {
    ind.prior <- priority[q]
    eventq <- tdata[, 1 + n_ep + ind.prior]

    if (q == 1L) {
      cen_timeq <- tdata[, 1 + ind.prior]
    } else {
      trt_centimq_1 <- trt_con[, 1 + n_ep + ind.prior - 1L]
      con_centimq_1 <- trt_con[, 2 + 3 * n_ep + ind.prior - 1L]
      if (is.infinite(L)) {
        cen_timeq <- tdata[, 1 + ind.prior]
        Ieqqth <- as.numeric(trt_centimq_1 == con_centimq_1)
      } else {
        cen_timeq <- rep(L, n)
        Ieqqth <- as.numeric(trt_centimq_1 == L & con_centimq_1 == L)
      }
      Ieqif <- Ieqif * Ieqqth
    }

    idx_all <- step_index(cen_timeq)
    trt_index <- idx_all[trt_id]
    con_index <- idx_all[con_id]

    surv_trt <- extract_surv(trt_id, trt_index)
    surv_con <- extract_surv(con_id, con_index)

    deltai <- eventq[trt_id]
    deltaj <- eventq[con_id]

    double_weight <- surv_trt * surv_con
    double_weight[double_weight < 0.01] <- 0.01
    IPCW_W[, q] <- Ieqif * (deltai * deltaj) / double_weight

    G.i <- extract_G(trt_id, trt_index)
    G.j <- extract_G(con_id, con_index)
    G.ij <- G.i + G.j
    epsilon_W[, q] <- as.numeric(G.ij %*% Omega_inv %*% matrix(mU, ncol = 1L))
  }

  weight_est <- rowSums(IPCW_W * determined_flag)
  weight_var0 <- rowSums(epsilon_W * determined_flag)
  weight_var <- weight_est * (1 + weight_var0)

  res <- list(weight_est = weight_est, weight_var = weight_var)
  return(res)
}
