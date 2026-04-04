#' Generate a simulated data for generalized win fraction regression model with identity link
#' 
#' This function can generated a simulated data for generalized win fraction regression model with identity link,
#' and two composite survival outcomes are generated from exponential distributions with different shape
#' parameters, and then modified through linear predictors with different regression coefficients
#' 
#' @param n sample size in the simulated data
#' @param beta_d true regression coefficient in linear predictors for the first endpoint
#' @param beta_h true regression coefficient in linear predictors for the second endpoint
#' @param lambda_d shape parameter in exponential distribution for the first endpoint
#' @param lambda_h shape parameter in exponential distribution for the second endpoint
#' @param betac true regression coefficients in Cox model for censoring
#' 
#' @return A data frame including id, y1(death time), y2(hospitalization time), 
#'         delta1(censored indicator for y1), delta2(censored indicator for y2),
#'         x(covariate used in gwrm) and z(covariate used in censoring model)
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @importFrom stats rbinom 
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom copula normalCopula
#' @importFrom copula rMvdc
#' 
#' 
#' @examples
#' n = 200  
#' beta_d = 10 
#' beta_h = 8
#' lambda_d = 0.15
#' lambda_h = 0.3 
#' betac = -5.7;  
#' mydata = sim_gwfm_identity(n,beta_d,beta_h,lambda_d,lambda_h,betac)
#' 

sim_gwfm_identity <- function(n,beta_d,beta_h,lambda_d,lambda_h,betac){
  #covaraite, binary veriable
  x = rbinom(n,1,0.5)
  #event time
  dataz = rmvnorm(n, mean = c(0, 0), sigma = matrix(c(1, 0.5, 0.5, 1),2,2))
  u1 = pnorm(dataz[,1]) #
  u2 = pnorm(dataz[,2]) #
  D = -log(1-u1)/lambda_d
  H = -log(1-u2)/lambda_h
  D = D + x*beta_d
  H = H + x*beta_h
  z = pnorm(x)
  C = -log(runif(n))/exp(z*betac) #censoring time
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<=C & H<=Do] = 1
  status2[H>Do] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = data.frame(id=1:n, y1= Do, y2 = Ho, delta1 = status1,
                   delta2 = status2,x,z)
  #note: y1 can be treated as observed death time and y2 can be treated as hospitalization time
  return(data_set)
}
