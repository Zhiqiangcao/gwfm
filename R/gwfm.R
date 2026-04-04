#' A generalized win fraction regression model for composite survival data
#' 
#' This function fits a generalized win fraction regression model (gwfm), also known as GWFM. 
#' It can be used to fit many different flavours of models that can be reformulated as a gwfm. The most
#' general models are implemented, but the flexible formula interface allows you to specify a wide
#' variety of different models, like the \code{\link{pim}} package
#' 
#' @param formula An object of class formula (or one that can be coerced to that class): A symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data a data frame including id, multiple endpoints (e.g., y1, y2), corresponding censoring status (e.g., delta1, delta2) and covariates
#' @param SurvTime name of survival time for fatal endpoint (e.g., death)
#' @param Status name of status for fatal endpoint 
#' @param covariateC covariates used in Cox regression model for censoring when IPCW method is chosen.
#'        If IPCW = FALSE, then this setting is useless (default NULL)
#' @param priority composite endpoints ranked in descending order of importance (e.g., priority = c(1,2) means y1 is more important than y2)
#' @param L pre-specified time horizon to construct restriction event endpoints (default Inf)
#' @param t_val truncation value used in IPCW weight (default 0.01)
#' @param IPCW using IPCW adjustment or not (default TRUE)
#' @param start start value for fitting a gwfm model for composite survival data (default NULL)
#' @param control control settings for \code{\link[nleqslv]{nleqslv}} when fitting gwfm model for composite survival data (default list())
#' @param link a character vector with a single value that determines the used link function. Possible values are "logit", "probit" and "identity". The default is "logit". 
#' @param model a single character value with possible values "difference" (the default), "marginal", or "customized". If the formula indicates a customized model (by the
#'        use of L() or R() in \code{\link{pim}} package), this parameter is set automatically to "customized".
#' @param compare a character vector with a single value that describes how the model compares observations. It can take the values "unique" or "all". Alternatively you can pass
#'        a matrix with two columns. Each row represents the row numbers in the original data frame that should be compared to each other. See Details.      
#' @param na.action the name of a function which indicates what should happen when the data contains NAs. The default is set by the na.action setting of options, and is na.fail when unset.
#' @param keep.data a logical value indicating whether the model matrix should be saved in the object. Defaults to FALSE. See Details
#' @param ... extra parameters sent to \code{\link{gwfm_fit}}
#' 
#' @details
#' GWFMs are based on a set of pseudo-observations constructed from the comparison between a range
#' of possible combinations of 2 observations with composite survival endpoints. We call the set of pseudo observations poset in the
#' context of this package, like \code{\link{pim}} package. 
#' 
#' By default, this poset takes every unique combination of 2 observations (compare = "unique"). You
#' can either use a character value, or use a matrix or list to identify the set of observation pairs that
#' have to be used as pseudo-observations. Note that the matrix and list should be either nameless,
#' or have the (col)names ’L’ and ’R’. If any other names are used, these are ignored and the first
#' column/element is considered to be ’L’. See also \code{\link{new.pim.poset}} in \code{\link{pim}} package.
#' 
#' It’s possible to store the model matrix and psuedo responses in the resulting object. By default this
#' is not done (keep.data = FALSE) as this is less burden on the memory and the gwfm object
#' contains all information to reconstruct both the model matrix and the pseudo responses. If either
#' the model matrix or the pseudo responses are needed for further calculations, setting keep.data to
#' TRUE might reduce calculation time for these further calculations.
#' 
#' @return A list with the following elements
#' \describe{
#'  \item{summary}{point estimation, SE and P value of parameter coefficients in gwfm}
#'  \item{P_tie}{tied probability}
#'  \item{coef}{point estimation of parameter coefficients in gwfm}
#'  \item{vcov}{a numeric matrix with the variance-covariance matrix for the coefficients}
#'  \item{fitted}{a numeric vector with the raw fitted values}
#'  \item{flag}{a flag indicating nleqslv is converged or not, its value is equal to termcd of nleqslv}
#'  \item{penv}{a gwfm environment from a model or formula}
#'  \item{link}{link function used in gwfm}
#'  \item{estim}{a list with two components named \code{coef} and \code{vcov}
#'  containing information on the used estimators for both.}
#'  \item{model.matrix}{a design matrix for a gwfm model}
#'  \item{response}{pseudo-observed response for win function. 1 for win result, 0 for loss or tied comparison result}
#'  \item{keep.data}{a logical value indicating whether the model matrix should be saved in the object}
#'  \item{model}{the chosen model in a gwfm}
#' }
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' @references Cao Z., Fang X. and Li F. Generalized win fraction regression for composite survival endpoints. 
#' under review. 2026;0(0):1-35.
#' 
#' @importFrom stats terms
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom pim new.pim.env
#' @importFrom pim new.pim.formula
#' 
#' @examples
#' # example 1: using gwfm with logit link
#' library(gumbel)
#' library(MASS)
#' library(pim)
#' library(copula)
#' library(nleqslv)
#' library(survival)
#' library(mvtnorm)
#' n = 200  
#' beta = c(-0.3,0.3) 
#' betac = c(-0.6,0.4) 
#' lambda_d = 0.2
#' lambda_h = 2
#' alpha = 2
#' set.seed(123456)
#' mydata1 = sim_gwfm_logit(n,beta,betac,lambda_d,lambda_h,alpha)
#' 
#' #no restriction event time, default L=Inf
#' res1 = gwfm(formula=y1~x1+x2,data=mydata1,SurvTime="y1",Status="delta1",
#' covariateC=c("x1","x2"),IPCW=FALSE, priority = 1:2, 
#' model = "difference",link = "logit")
#' print(res1$summary)
#' 
#' # restriction event time and Broyden method in nleqslv
#' res2 = gwfm(formula=y1~x1+x2,data=mydata1,SurvTime="y1",Status="delta1",
#' covariateC=c("x1","x2"), L = 1, IPCW=FALSE, priority = 1:2,
#'  model = "difference",link = "logit", method="Broyden")
#' print(res2$summary)
#' print(res2$WO)
#' 
#' # example 2: using gwfm with probit link
#' n = 200  
#' alpha = 1
#' betac = -1
#' u = 6
#' sig = 1
#' set.seed(123456)
#' mydata2 = sim_gwfm_probit(n,alpha,betac,u,sig)
#' 
#' # NO IPCW 
#' xtol = 1e-6; ftol = 1e-6
#' btol = 1e-3; maxit = 50
#' cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)
#' res3 = gwfm(formula=y1~x,data=mydata2,SurvTime="y1",Status="delta1",
#' covariateC=c("x"), priority = 1:2, L = 7, IPCW = FALSE,
#' model = "difference",link = "probit",control = cntl)
#' print(res3$summary)
#' 
#' # IPCW for restriction event time L=7
#' res4 = gwfm(formula=y1~x,data=mydata2,SurvTime="y1",Status="delta1",
#' covariateC=c("x"),priority = 1:2, L = 7, model = "difference",
#' link = "probit",control = cntl)
#' print(res4$summary)
#' 
#' # IPCW for restriction event time L=8
#' res5 = gwfm(formula=y1~x,data=mydata2,SurvTime="y1",Status="delta1",
#' covariateC=c("x"), priority = 1:2, L = 8, model = "difference",
#' link = "probit", control = cntl)
#' print(res5$summary)
#' 
#' # example 3: using gwfm with identity link
#' n = 200  
#' beta_d = 10 
#' beta_h = 8
#' lambda_d = 0.15
#' lambda_h = 0.3 
#' betac = -5.7
#' set.seed(123456) 
#' mydata3 = sim_gwfm_identity(n,beta_d,beta_h,lambda_d,lambda_h,betac)
#' 
#' # NO IPCW 
#' res6 = gwfm(formula=y1~x,data=mydata3,SurvTime="y1",Status="delta1",
#' covariateC=c("z"), priority = 1:2, L = 15, IPCW = FALSE, 
#' model = "marginal",link = "identity",control = cntl)
#' print(res6$summary)
#' 
#' # IPCW for restriction time L=15
#' res7 = gwfm(formula=y1~x,data=mydata3,SurvTime="y1",Status="delta1",
#' covariateC=c("z"), priority = 1:2, L = 15, model = "marginal",
#' link = "identity",control = cntl)
#' print(res7$summary)
#' 
#' # IPCW for restriction time L=12
#' res8 = gwfm(formula=y1~x,data=mydata3,SurvTime="y1",Status="delta1",
#' covariateC=c("z"), priority = 1:2, L = 12, model = "marginal",
#' link = "identity",control = cntl)
#' print(res8$summary)
#' 

gwfm <- function(formula, data, SurvTime, Status, covariateC = NULL, priority = 1:2, 
                       L = Inf, t_val = 0.01, IPCW = TRUE, start = NULL, 
                       control = list(), link = c("logit", "probit", "identity"),
                       model = c("difference", "marginal", "customized"), 
                       compare = if (model == "marginal") "all" else "unique", 
                       na.action = getOption("na.action"), keep.data = FALSE, ...){
  
  model = match.arg(model)
  if (is.character(compare)) {
    if (!compare %in% c("unique", "all")) 
      stop("compare should have the value 'unique' or 'all' in case it's a character value.")
  }
  
  # sort data by observed survival time
  index = order(data[,SurvTime])
  data.sort = data[index,]
  data = data.sort
  data = as.data.frame(data)
  link = match.arg(link)
  if (is.null(na.action)) na.action = "na.fail"
  if (!is.character(na.action)) na.action = deparse(substitute(na.action))
  
  f.terms = terms(formula, simplify = TRUE)
  vars = all.vars(formula)
  # prepare covariate matrix
  penv = new.pim.env(data, compare = compare, vars = vars, env = parent.frame())
  ff = new.pim.formula(formula, penv)
  x = model.matrix.gwfm(ff, na.action = na.action, model = model)
    
  # next, we prepare comparison results between composite survival endpoints 
  # length of endpoints
  n_ep = length(priority)
  # since the data set format is id, y1, ..., y_nep, delta1, ..., delta_nep, z1,...,zp
  datap = as.matrix(data[,c(1:(1+2*n_ep))])
  # construct restriction event times and data if L != Inf
  y_o = as.matrix(data[,2:(n_ep+1)]) #original observed endpoints
  delta_o = as.matrix(data[,(n_ep+2):(1+2*n_ep)]) #original observed censoring indicators
  y_n = pmin(y_o,L) #restriction event times
  delta_n = ifelse(y_n == L, 1, delta_o)
  datap_n = as.matrix(cbind(data[,1], y_n, delta_n))
  #merge data similar to compare treatment and control groups
  trt_con = merge_csd(datap = datap_n, model = model, n_ep = n_ep)
  #next, we compute pairwise comparison result based on Dong et al's (2020) method
  win_status = win_strategy(trt_con = trt_con, priority = priority)
  #find out win function is determined by which endpoint comparison
  determined_flag = NULL
  for(q in seq(3,(2+2*n_ep-1),2)){
    win_qth_res = win_status[,q:(q+1)]
    qth_flag = apply(as.matrix(win_qth_res,ncol = 2),1,sum)
    determined_flag = cbind(determined_flag,qth_flag)
  }
  #find out index of ties comparison (i.e., ties_flag=0)
  ties_flag = apply(determined_flag,1,sum)
  #y is the win result of W(Y_i,Y_j)
  y = apply(as.matrix(win_status[,seq(3,(2+2*n_ep-1),2)],ncol = n_ep),1,sum)
  #make those ties comparison results do not contribute in score function
  x = as.matrix(ties_flag*x)
  p = dim(x)[2]
  p_tie = sum(ties_flag==0)/length(y)  #tied probability
  #when model = "difference" or "marginal", penv is also different
  penv = as.environment(penv@poset)
  if(is.null(start)){  #using zero as default start points
    start = rep(0, p)
  }
  if(IPCW == TRUE){ 
    tdata = data.frame(datap_n,data[,covariateC]) #truncated
    colnames(tdata) = c(colnames(data[,c(1:(1+2*n_ep))]),covariateC)
    censor_cox = censoring_cox_est(Data=tdata,SurvTime=SurvTime,Status=Status,covariateC=covariateC)
    weight_res = ipcw_weight(tdata,trt_con,censor_cox,win_status,priority,determined_flag,L,n_ep,t_val=t_val)
    weight_est = weight_res$weight_est
    weight_var = weight_res$weight_var
  }else{ #NO IPCW
    weight_est = weight_var = NULL
  }
  #fit score and sandwich variance estimation
  est_res = gwfm_fit(x=x, y=y, link = link, start = start, control = control, 
                       weight_est = weight_est, weight_var = weight_var, penv = penv, ...)
  #summary estimation results
  names(est_res$coef) = colnames(x)
  if (!keep.data) {x = y = NULL}
  flag = est_res$flag
  Estimate = est_res$coef
  Std.Error = sqrt(diag(est_res$vcov))
  if(sum(is.na(Std.Error)>0) | sum(is.infinite(Std.Error))>0){
    flag = 0;
    Std.Error = rep(1,p)
  }
  z.value = Estimate/Std.Error
  p.value = rep(0,p)
  for(i in 1:p){
    z.valuei = z.value[i]
    p.value[i] = 2*pnorm(-abs(z.valuei))
  }
  para_est_res = data.frame(Estimate,Std.Error,z.value,p.value)
  final_res = list(summary = para_est_res, P_tie = p_tie, coef = est_res$coef, 
                   vcov = est_res$vcov, fitted = est_res$fitted, flag = flag, 
                   penv = penv, link = link, estimators = est_res$estim, model.matrix = x, 
                   response = y, keep.data = keep.data, model = model)
  return(final_res)
}
