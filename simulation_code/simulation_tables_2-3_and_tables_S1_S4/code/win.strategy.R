# Hierarchical observed strict win function for composite survival outcomes.
# This revised implementation is exact for the death + one nonfatal endpoint
# setting used in the current simulations. It also enforces sequential
# comparison so that lower-priority endpoints are only considered when all
# higher-priority endpoints are tied on the observed restricted times.
win.strategy <- function(trt_con, priority){
  n_ep <- length(priority)
  npair <- nrow(trt_con)

  get_col <- function(prefix, ep, arm) {
    paste0(prefix, "_", ep, "_", arm)
  }

  win_mat <- matrix(0, nrow = npair, ncol = 2L * n_ep)
  eligible <- rep(TRUE, npair)

  for (q in seq_len(n_ep)) {
    ep <- priority[q]
    delta_trt <- trt_con[[get_col("Delta", ep, "trt")]]
    delta_con <- trt_con[[get_col("Delta", ep, "con")]]
    y_trt <- trt_con[[get_col("Y", ep, "trt")]]
    y_con <- trt_con[[get_col("Y", ep, "con")]]

    win_trt <- eligible & (delta_con == 1) & (y_trt > y_con)
    win_con <- eligible & (delta_trt == 1) & (y_con > y_trt)

    win_mat[, 2L * q - 1L] <- as.integer(win_trt)
    win_mat[, 2L * q] <- as.integer(win_con)

    # Only exact ties on the observed restricted endpoint allow the next
    # lower-priority endpoint to be considered.
    eligible <- eligible & !(win_trt | win_con) & (y_trt == y_con)
  }

  colnames(win_mat) <- paste0(rep(c("Trt", "Con"), n_ep), "_Endpoint", rep(priority, each = 2L))
  win_status <- as.data.frame(win_mat)
  win_status <- data.frame(
    pid_trt = trt_con[, 1],
    pid_con = trt_con[, 2 + 2 * n_ep],
    win_status,
    stringsAsFactors = FALSE
  )
  return(win_status)
}
