suppressPackageStartupMessages({
  library(stats)
})

# -----------------------------
# Covariate generation
# -----------------------------

r_x1 <- function(n) {
  lower <- pnorm(-1)
  upper <- pnorm(1)
  qnorm(runif(n, min = lower, max = upper))
}

r_x2 <- function(n) {
  sample(c(-1, 1), size = n, replace = TRUE)
}

sample_covariates <- function(n) {
  data.frame(
    x1 = r_x1(n),
    x2 = r_x2(n),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Gumbel-Hougaard generator
# -----------------------------

r_positive_stable <- function(n, gamma) {
  if (abs(gamma - 1) < 1e-12) return(rep(1, n))
  if (gamma <= 0 || gamma > 1) stop("gamma must lie in (0, 1].")

  U <- runif(n, 0, pi)
  W <- rexp(n, 1)
  part1 <- sin(gamma * U) / (sin(U))^(1 / gamma)
  part2 <- (sin((1 - gamma) * U) / W)^((1 - gamma) / gamma)
  part1 * part2
}

generate_latent_subjects <- function(n,
                                     alpha,
                                     lambda_D,
                                     lambda_H,
                                     beta_D,
                                     beta_H,
                                     seed = NULL,
                                     x = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta_D) != 2L || length(beta_H) != 2L) {
    stop("beta_D and beta_H must both have length 2.")
  }

  if (is.null(x)) {
    x <- sample_covariates(n)
  }
  if (nrow(x) != n) stop("nrow(x) must equal n.")

  xmat <- as.matrix(x[, c("x1", "x2")])
  rate_D <- lambda_D * exp(-drop(xmat %*% beta_D))
  rate_H <- lambda_H * exp(-drop(xmat %*% beta_H))

  if (abs(alpha - 1) < 1e-12) {
    D <- rexp(n, rate = rate_D)
    H <- rexp(n, rate = rate_H)
  } else {
    gamma <- 1 / alpha
    V <- r_positive_stable(n, gamma)
    E1 <- rexp(n, 1)
    E2 <- rexp(n, 1)
    U1 <- exp(- (E1 / V)^gamma)
    U2 <- exp(- (E2 / V)^gamma)
    U1 <- pmin(pmax(U1, 1e-12), 1 - 1e-12)
    U2 <- pmin(pmax(U2, 1e-12), 1 - 1e-12)
    D <- -log(U1) / rate_D
    H <- -log(U2) / rate_H
  }

  data.frame(
    x1 = x$x1,
    x2 = x$x2,
    D = D,
    H = H,
    stringsAsFactors = FALSE
  )
}

generate_censoring_times <- function(x,
                                     lambda_C,
                                     beta_C,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta_C) != 2L) stop("beta_C must have length 2.")

  xmat <- as.matrix(x[, c("x1", "x2")])
  rate_C <- lambda_C * exp(drop(xmat %*% beta_C))
  C <- rexp(nrow(x), rate = rate_C)

  data.frame(C = C, rate_C = rate_C, stringsAsFactors = FALSE)
}

latent_to_observed <- function(latent_df, censor_df = NULL) {
  if (is.null(censor_df)) {
    censor_df <- data.frame(C = rep(Inf, nrow(latent_df)), stringsAsFactors = FALSE)
  }
  if (nrow(latent_df) != nrow(censor_df)) {
    stop("latent_df and censor_df must have the same number of rows.")
  }

  D <- latent_df$D
  H <- latent_df$H
  C <- censor_df$C

  y1 <- pmin(D, C)
  delta1 <- as.integer(D <= C)

  y2 <- pmin(H, D, C)
  delta2 <- as.integer(H <= pmin(D, C))

  data.frame(
    id = seq_len(nrow(latent_df)),
    y1 = y1,
    y2 = y2,
    delta1 = delta1,
    delta2 = delta2,
    x1 = latent_df$x1,
    x2 = latent_df$x2,
    stringsAsFactors = FALSE
  )
}

generate_observed_subjects <- function(n,
                                       alpha,
                                       lambda_D,
                                       lambda_H,
                                       lambda_C,
                                       beta_D,
                                       beta_H,
                                       beta_C,
                                       seed = NULL,
                                       return_latent = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  x <- sample_covariates(n)
  latent <- generate_latent_subjects(
    n = n,
    alpha = alpha,
    lambda_D = lambda_D,
    lambda_H = lambda_H,
    beta_D = beta_D,
    beta_H = beta_H,
    x = x
  )
  censor_df <- generate_censoring_times(
    x = x,
    lambda_C = lambda_C,
    beta_C = beta_C
  )
  observed <- latent_to_observed(latent, censor_df)

  if (isTRUE(return_latent)) {
    return(list(observed = observed, latent = cbind(latent, censor_df)))
  }
  observed
}

# -----------------------------
# Restriction and pair outcomes
# -----------------------------

restrict_observed_data <- function(dat, L) {
  xi1 <- pmin(dat$y1, L)
  d1 <- ifelse(is.finite(L) & xi1 == L, 1L, dat$delta1)

  xi2 <- pmin(dat$y2, L)
  d2 <- ifelse(is.finite(L) & xi2 == L, 1L, dat$delta2)

  data.frame(
    id = dat$id,
    x1 = dat$x1,
    x2 = dat$x2,
    y1 = dat$y1,
    y2 = dat$y2,
    delta1 = dat$delta1,
    delta2 = dat$delta2,
    xi1 = xi1,
    xi2 = xi2,
    d1 = d1,
    d2 = d2,
    stringsAsFactors = FALSE
  )
}

compute_pair_outcomes <- function(rest_df) {
  n <- nrow(rest_df)
  if (n < 2L) stop("Need at least two subjects.")

  idx <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]

  xi1_i <- rest_df$xi1[i]
  xi1_j <- rest_df$xi1[j]
  d1_i <- rest_df$d1[i]
  d1_j <- rest_df$d1[j]

  xi2_i <- rest_df$xi2[i]
  xi2_j <- rest_df$xi2[j]
  d2_i <- rest_df$d2[i]
  d2_j <- rest_df$d2[j]

  win1 <- (xi1_i > xi1_j) & (d1_j == 1L)
  lose1 <- (xi1_j > xi1_i) & (d1_i == 1L)
  tie1 <- (xi1_i == xi1_j)

  win2 <- tie1 & (xi2_i > xi2_j) & (d2_j == 1L)
  lose2 <- tie1 & (xi2_j > xi2_i) & (d2_i == 1L)

  resolved_q1 <- win1 | lose1
  resolved_q2 <- (!resolved_q1) & (win2 | lose2)
  resolved <- resolved_q1 | resolved_q2
  y_win <- as.integer(win1 | win2)

  data.frame(
    i = i,
    j = j,
    z1 = rest_df$x1[i] - rest_df$x1[j],
    z2 = rest_df$x2[i] - rest_df$x2[j],
    y = y_win,
    resolved = as.integer(resolved),
    resolved_q1 = as.integer(resolved_q1),
    resolved_q2 = as.integer(resolved_q2),
    xi1_i = xi1_i,
    xi1_j = xi1_j,
    xi2_i = xi2_i,
    xi2_j = xi2_j,
    d1_i = d1_i,
    d1_j = d1_j,
    d2_i = d2_i,
    d2_j = d2_j,
    stringsAsFactors = FALSE
  )
}

generate_truth_pairs <- function(n_pair,
                                 alpha,
                                 lambda_D,
                                 lambda_H,
                                 beta_D,
                                 beta_H,
                                 L,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  x_i <- sample_covariates(n_pair)
  x_j <- sample_covariates(n_pair)

  lat_i <- generate_latent_subjects(
    n = n_pair,
    alpha = alpha,
    lambda_D = lambda_D,
    lambda_H = lambda_H,
    beta_D = beta_D,
    beta_H = beta_H,
    x = x_i
  )
  lat_j <- generate_latent_subjects(
    n = n_pair,
    alpha = alpha,
    lambda_D = lambda_D,
    lambda_H = lambda_H,
    beta_D = beta_D,
    beta_H = beta_H,
    x = x_j
  )

  obs_i <- latent_to_observed(lat_i)
  obs_j <- latent_to_observed(lat_j)

  ri <- restrict_observed_data(obs_i, L)
  rj <- restrict_observed_data(obs_j, L)

  win1 <- (ri$xi1 > rj$xi1) & (rj$d1 == 1L)
  lose1 <- (rj$xi1 > ri$xi1) & (ri$d1 == 1L)
  tie1 <- (ri$xi1 == rj$xi1)

  win2 <- tie1 & (ri$xi2 > rj$xi2) & (rj$d2 == 1L)
  lose2 <- tie1 & (rj$xi2 > ri$xi2) & (ri$d2 == 1L)

  resolved_q1 <- win1 | lose1
  resolved_q2 <- (!resolved_q1) & (win2 | lose2)
  resolved <- resolved_q1 | resolved_q2
  y_win <- as.integer(win1 | win2)

  data.frame(
    z1 = x_i$x1 - x_j$x1,
    z2 = x_i$x2 - x_j$x2,
    y = y_win,
    resolved = as.integer(resolved),
    stringsAsFactors = FALSE
  )
}

subject_censor_rate_at_L <- function(latent_df, L) {
  if (!all(c("D", "C") %in% names(latent_df))) {
    stop("latent_df must contain D and C.")
  }
  if (is.infinite(L)) {
    mean(latent_df$C < latent_df$D)
  } else {
    mean(latent_df$C < pmin(latent_df$D, L))
  }
}
