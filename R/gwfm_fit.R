#' Fitter function for a generalized win fraction regression model for composite outcomes
#' 
#' This is the basic computing engine called by \code{\link{gwfm}} to get the estimates for the coefficients and the variance-
#' covariance matrices. This function currently only spits out these components using the sandwich estimators.
#'
#' @param x a model matrix with as many rows as \code{y}.
#' @param y a vector with the pseudo-responses
#' @param link a character vector with a link function
#' @param start A numeric vector with the starting values for the fitting algorithm, if required.
#' @param control conditions for convergence of \code{\link[nleqslv]{nleqslv}}
#' @param method "Newton" (default) and "Broyden" can be chosen in \code{\link[nleqslv]{nleqslv}}
#' @param vcov.estim a function to determine the variance-covariance matrix.
#' Possibilities are \code{\link{sandwich.vcov}} and \code{link{score.vcov}}.
#' Defaults to \code{sandwich.vcov}
#' @param weight_est weights used in parameter estimation
#' @param weight_var weights used in variance estimation
#' @param penv An environment
#' @param ... Further arguments that need to be passed to the estimation function. The most
#' relevant is construct, allowing you to write your own score function for the
#' numerical optimization. See also estimators
#'
#' @return A list with the following elements
#' \describe{
#'  \item{coefficients}{a numeric vector with the coefficients}
#'  \item{fitted}{a numeric vector with the raw fitted values}
#'  \item{flag}{a flag indicating nleqslv is converged or not, its value is equal to termcd of nleqslv}
#'  \item{estim}{a list with two components named \code{coef} and \code{vcov}
#'  containing information on the used estimators for both.}
#'  \item{vcov}{a numeric matrix with the variance-covariance matrix for
#'  the coefficients}
#' }
#' 
#' @section NOTE: This function is mainly quoted from \code{\link{pim}} package by Joris Meys et al.(2020) 
#'          (https://cran.r-project.org/web/packages/pim/index.html)
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' @importFrom nleqslv nleqslv
#' @importFrom pim poset
#' 
#' 


gwfm_fit <- function(x, y, link = "logit", start = rep(0, ncol(x)), 
                       control = list(), method = "Newton", vcov.estim = "sandwich.vcov", 
                       weight_est = NULL, weight_var = NULL, penv,...){
  fn <- CreateScoreFun(Z=x, Y=y, link=link, W=weight_est)
  res <- nleqslv(start, fn, method = method, control = control)
  fits <- x %*% res$x
  flag <- res$termcd    #1-converge; others-not fully converged
  vcov.estimF <- match.fun(vcov.estim)
  dim(fits) <- NULL
  if (!is.list(penv)) penv <- poset(penv, as.list = TRUE)
  vc <- vcov.estimF(fitted=fits, X=x, Y=y, W=weight_var, link=link, poset=penv)
  return(list(coefficients = res$x, vcov = vc, fitted = fits, flag = flag,
              estim = list(coef = "estimator.nleqslv", vcov = vcov.estim)))
}
