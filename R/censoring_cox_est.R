#' Compute survival function and other related results of censoring based on Cox proportional hazards model
#' 
#' This function can obtain survival function of censoring, G, U and Omega results 
#' based on Cox model. Their concrete formulas are shown in Section 3 of paper (Cao et al., 2026).
#' 
#' @param Data A data frame includes survival time, censoring indicator and related covariates for censoring
#' @param SurvTime name of survival time for fatal endpoint (e.g., death) in Data
#' @param Status name of status for fatal endpoint in Data
#' @param covariateC a vector of covariates used in Cox model for censoring
#' 
#' @return A list with the following elements:
#' \describe{
#'  \item{Kij}{a matrix for survival functions of censoring, and the element in (i,j) refers to the survival function of censoring for ith subject at observed time t_j}
#'  \item{Omega}{a variance-covariance matrix for coefficients in Cox model for cesoring}
#'  \item{G}{value of G for subject i at its observed time t_i, a n times pc (dimension of covariateC) matrix}
#'  \item{Psi}{value of martingle integral for subject i, a n times pc (dimension of covariateC) matrix}
#' }
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @references Cao Z., Fang X. and Li F. Generalized win fraction regression for composite survival endpoints. 
#' under review. 2026;0(0):1-35.
#' 
#' 
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
#' @importFrom survival coxph
#' 


censoring_cox_est <- function(Data,SurvTime,Status,covariateC){
  #Data has been sorted by SurvTime
  censor.formula = as.formula(paste("Surv(",SurvTime,",I(1-",Status,"))","~",
                                    paste(covariateC,collapse="+"),sep="")) 
  xc = model.matrix(censor.formula, data=Data)
  pc = dim(xc)[2]
  n = dim(xc)[1]
  if(pc == 1){    #in this case, no covariate in censoring model, since x from model.matrix is a n*1 matrix 
    xc = xc       #(intercept) with all values=1, then estimate of x (intercept) will be NA in censored model 
  }else{
    xc = matrix(xc[,-1],ncol=pc-1) #delete the first column of x, since that case is intercept and no intercept in Cox model
  }
  pc = dim(xc)[2]
  colnames(xc) = paste("xc",1:dim(xc)[2],sep="")
  time = as.numeric(Data[,SurvTime])
  status = as.numeric(Data[,Status])
  datac = as.data.frame(cbind(time,status,xc))
  colnames(datac) = c(SurvTime,Status,colnames(xc))
  #update censoring formula
  censor.formula = as.formula(paste("Surv(",SurvTime,",I(1-",Status,"))","~",
                                    paste(colnames(xc),collapse = "+"),sep="")) 
  cen.model = coxph(censor.formula, data = datac)
  
  ###the following part is to compute Sc_i(t_j),i.e., survival probability
  ###of individual i at observed time t_j
  bc_est = cen.model$coef
  delta = Data[,Status]
  if(pc == 1){ #only one covariate
    if(is.na(bc_est)){  #then no covariate is included in Cox model
      bc_est = matrix(0,1,1)  #reset estimated parameter coefficient's value to 0
    }
  }
  s2 = cumsum(exp(as.matrix(xc[n:1,])%*%bc_est)) #using as.matrix is to make sure 
  #when p=1, x[n:1,] is a matrix (since x[n:1,] is a numeric if we do not use as.matrix)
  dlc = (1-delta)/s2[n:1]
  dlc[(1:length(dlc))[is.nan(dlc)]] = 0
  ssc = matrix(0, n, n)
  xtbc_est = as.numeric(exp(xc%*%bc_est))
  ssc[,1] = dlc[1]*xtbc_est
  otc = outer(xtbc_est, dlc, "*") #otc=exp(x_1*b)*dlc[1],exp(x_1*b)*dlc[2],exp(x_1*b)*dlc[3]  
  otc[,1] = rep(0, n)             #exp(x_2*b)*dlc[1],exp(x_2*b)*dlc[2],exp(x_2*b)*dlc[3]   
  otc_cum = t(apply(otc,1,cumsum))#exp(x_3*b)*dlc[1],exp(x_3*b)*dlc[2],exp(x_3*b)*dlc[3] 
  ssc_1 = outer(ssc[,1],rep(1,n)) #otc_cum=row cumsum
  cumhcij = ssc_1 + otc_cum
  #for patient i with covariates x[i], compute estimated survival function 
  #of censoring at observed time t[j]      
  Kij = exp(-cumhcij) #is a n*n matrix
  #################
  
  ###The next part is to compute R_C^{(0)}, R_C^{(1)}, R_C^{(2)} and \Omega(\beta_C)
  #inrev = rev(index)
  xrev = as.matrix(xc[n:1,]) #xrev is equal to x[n:1,]
  ss0r = as.numeric(exp(xrev%*%bc_est))
  ss0rs = cumsum(ss0r) #a vector, length is n
  Rc0 = ss0rs[n:1] #return to original order (i.e., from the smallest to largest)
  
  ss1r = xrev*ss0r
  ss1rs = apply(ss1r,2,cumsum)
  Rc1 = as.matrix(ss1rs[n:1,]) #n*pc matrix
  
  xsrev = apply(xrev,1,function(y) y%*%t(y))
  xsrev = t(xsrev)
  ss2r = xsrev*ss0r
  ss2rs = apply(ss2r,2,cumsum)
  if(pc == 1){
    Rc2 = as.matrix(ss2rs[n:1])
  }else{
    Rc2 = ss2rs[n:1,] #n*(p^2) matrix
  }
  
  xbar = Rc1/Rc0
  if(pc == 1){
    xbars = xbar^2
  }else{
    xbars = apply(xbar,1,function(y) y%*%t(y))
    xbars = t(xbars)
  }
  omegai = (1-delta)*(Rc2/Rc0-xbars)
  omegam = apply(omegai,2,mean) 
  Omega = matrix(omegam,pc,pc) #estimated Omega
  ###########
  
  #he next part is to compute G_i(t_i), i.e., G's value of individual i at observed time t_i
  G1 = xc*(xtbc_est*cumsum(dlc))  #n*pc matrix
  xbarbcum = xbar*dlc  
  xbarbcumc = apply(xbarbcum,2,cumsum) 
  G2 = xtbc_est*xbarbcumc
  G = G1-G2
 
  
  ###The next part is to compute average of U_l(\beta_c), which is a p*1 matrix
  U1 = (xc-xbar)*(1-delta) #x, xbar and delta has been sorted by observed time
  #l = 1:n
  #U2 = matrix(0,n,pc)
  #for(i in 1:n){
  #  xtemp = xc[i,]
  #  xtm = matrix(rep(xtemp,n),byrow = T,ncol=pc)
  #  #compute 1:i terms
  #  Iil = 1*(i>=l)
  #  xic = xtbc_est[i]*Iil*dlc
  #  xtmm = (xtm-xbar)*xic
  #  U2[i,] = apply(xtmm,2,sum)
  #} 
  #note: U2=G by calculation formula
  U2 = G
  Psi = U1-U2  #n*p matrix
  
  res = list(Kij=Kij, Omega=Omega, G=G, Psi=Psi)
  return(res)
}
