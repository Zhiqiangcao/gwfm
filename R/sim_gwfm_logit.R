#' Generate a simulated data for generalized win fraction regression model with logit link
#' 
#' This function can generated a simulated data for generalized win fraction regression model with logit link,
#' and two composite survival outcomes are generated from a conditional Gumbel-Hougaard
#' copula model with Cox proportional hazards models as two marginal models
#' 
#' @param n sample size in the simulated data
#' @param beta true regression coefficients in gwfrm model
#' @param betac true regression coefficients in Cox model for censoring
#' @param lambda_d baseline hazard in the marginal Cox model for the death time 
#' @param lambda_h baseline hazard in the marginal Cox model for the hospitalization time
#' @param alpha Kendall’s tau between hospitalization time and death time 
#' 
#' @return A data frame including id, y1(death time), y2(hospitalization time), 
#'         delta1(censored indicator for y1), delta2(censored indicator for y2),
#'         x1 and x2(covariates)
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @importFrom stats rnorm 
#' @importFrom stats rbinom
#' @importFrom stats runif
#' @importFrom gumbel rgumbel
#' 
#' 
#' @examples
#' n = 200  
#' beta = c(-0.5,0.5) 
#' betac = c(-1,0.5) 
#' lambda_d = 0.2
#' lambda_h = 2
#' alpha = 2
#' mydata = sim_gwfm_logit(n,beta,betac,lambda_d,lambda_h,alpha)
#' 

sim_gwfm_logit <- function(n,beta,betac,lambda_d,lambda_h,alpha){
  x1 = rnorm(n)
  x1 = x1*(-1 <= x1 & x1 <= 1) - 1*(x1 < -1) + 1*(x1>1)
  x2 = 2*rbinom(n,1,0.5)-1
  outcome = rgumbel(n,alpha=alpha,dim=2)
  Lambda_D = as.numeric(lambda_d*exp(-cbind(x1,x2)%*%beta))
  Lambda_H = as.numeric(lambda_h*exp(-cbind(x1,x2)%*%beta))
  D = -log(outcome[,1])/Lambda_D #survival time
  H = -log(outcome[,2])/Lambda_H #nonfatal event time
  #censoring time
  C = -log(runif(n))/exp(cbind(x1,x2)%*%betac) #censoring time
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<C & H<Do] = 1
  status2[H>pmin(Do,C)] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = data.frame(id=1:n, y1= Do, y2 = Ho, delta1 = status1,
                   delta2 = status2, x1, x2)
  #note: y1 can be treated as observed death time and y2 can be treated as hospitalization time
  return(data_set)
}
