#' Generate a simulated data for generalized win fraction regression model with probit link
#' 
#' This function can generated a simulated data for generalized win fraction regression model with probit link,
#' and two composite survival outcomes are generated from a conditional Gaussian copula model.
#' 
#' @param n sample size in the simulated data
#' @param alpha true regression coefficients in marginal model
#' @param betac true regression coefficients in Cox model for censoring
#' @param u intercept term in both marginal models
#' @param sig standard error used in Gaussian copula model 
#' 
#' @return A data frame including id, y1(death time), y2(hospitalization time), 
#'         delta1(censored indicator for y1), delta2(censored indicator for y2),
#'         x1 and x2(covariates)
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @importFrom stats rnorm 
#' @importFrom stats runif
#' @importFrom copula mvdc
#' @importFrom copula normalCopula
#' @importFrom copula rMvdc
#' 
#' 
#' @examples
#' n = 200  
#' alpha = 1
#' betac = -0.79
#' u = 6
#' sig = 1
#' set.seed(123456)
#' mydata = sim_gwfm_probit(n,alpha,betac,u,sig)
#' 

sim_gwfm_probit <- function(n,alpha,betac,u,sig){
  x = seq(3,u,length=n)
  y = alpha*x
  mv.NE = mvdc(normalCopula(0.75), c("norm", "norm"),
               list(list(mean = 0, sd =sig), list(mean = 0, sd = sig)))
  x.sample = rMvdc(n, mv.NE)  ##error term
  D = y + x.sample[,1]
  H = y + x.sample[,2]
  C = -log(runif(n))/exp(x*betac) #censoring time
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
                        delta2 = status2, x)
  #note: y1 can be treated as observed death time and y2 can be treated as hospitalization time
  return(data_set)
}
