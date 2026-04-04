#' Calculate IPCW weights in score function and sandwich variance estimation for a generalized win fraction regression model
#' 
#' This function can compute IPCW weights used in score function and variance estimation for 
#' a generalized win fraction regression model. Their detailed expressions are shown in Section 2.3 and Section 3 of paper (Cao et al., 2026).
#' 
#' @param tdata A data frame includes survival time, censoring indicator and related covariates for censoring
#' @param trt_con Given n subjects, there are n(n-1)/2 pairs for non-marginal structure and n(n-1) pairs for marginal structure, each row represents a pair. The analysis dataset trt_con contains the following variables:
#' \describe{
#'  \item{pid_trt}{A vector for the subject id of the individuals from the "treatment" group}
#'  \item{Delta_j_trt}{A vector for the event status of the jth endpoint (1=event, 0=censored) for the individuals from the "treatment" group}
#'  \item{Y_j_trt}{A vector for the outcome of the jth endpoint for the individuals from the "treatment" group}
#'  \item{pid_con}{A vector for the subject id of the individuals from the "control" group}
#'  \item{Delta_j_con}{A vector for the event status of the jth endpoint (1=event, 0=censored) for the individuals from the "control" group}
#'  \item{Y_j_con}{A vector for the outcome of the jth endpoint for the individuals from the "control" group}
#' }
#' @param censor_cox a list results from \code{\link{censoring_cox_est}}
#' @param win_status A data.frame for the win status of each pair for each endpoint, which is from \code{\link{win_strategy}}
#' @param priority Priority order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 1, Endpoint 2, and Endpoint 3, input priority = c(1,2,3).
#' @param determined_flag a value to indicate pairwise comparison is determined or not (1 refers to the comparison is determined; 0 refers to a tied comparison finally)
#' @param L restriction event time
#' @param n_ep the length of multiple endpoints
#' @param t_val truncation value for IPCW weights, which is a small value users can set
#' 
#' @return A list with the IPCW weights used in score function and variance estimation
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @references Cao Z., Fang X. and Li F. Generalized win fraction regression for composite survival endpoints. 
#' under review. 2026;0(0):1-35.
#' 
#' @importFrom MASS ginv
#' 

ipcw_weight <- function(tdata,trt_con,censor_cox,win_status,priority,determined_flag,L,n_ep,t_val){
  #first, we compute double weights for estimation
  n = dim(tdata)[1] 
  # Setup survival object (based on the Death time)
  y1 = tdata[,2]
  u = c(0,y1)
  #obtain survival estimation at each observed time, from the smallest to the largest
  Kij = diag(censor_cox$Kij)
  Kij_ori = c(1,Kij)
  L.tpt = length(u)
  npim = dim(win_status)[1]
  #the following is used to compute weight_var
  Omega = censor_cox$Omega
  Omega_inv = ginv(Omega)  #generalized inverse
  Psi = censor_cox$Psi
  mPsi = apply(Psi,2,mean)
  pc = length(mPsi) #number of covariates in censoring model
  G = censor_cox$G  #is a n*pc matrix
  pc = dim(G)[2]
  G_ori = rbind(rep(0,pc),G)
  #comparison index
  con_id = win_status$pid_con
  trt_id = win_status$pid_trt
  IPCW_W = epislon_W = matrix(0,npim,n_ep)
  Ieqif = rep(1,npim)  #for the 1th endpoint   
  for(q in 1:n_ep){
    ind.prior = priority[q]
    eventq = tdata[,(1+n_ep+ind.prior)]
    if(q == 1){ #no matter restriction or not
      cen_timeq = tdata[,(1+ind.prior)]
    }else{ #q>=2
      trt_centimq_1 = trt_con[,1+n_ep+ind.prior-1] #observed survival time of the q-1 th endpoint
      con_centimq_1 = trt_con[,2+3*n_ep+ind.prior-1]
      if(L == Inf) { #no restriction
        cen_timeq = tdata[,(1+ind.prior)]
        Ieqqth = 1*(trt_centimq_1==con_centimq_1)
      }else{ #truncated
        cen_timeq = rep(L,n)
        Ieqqth = 1*(trt_centimq_1==L & con_centimq_1==L)
      }
      Ieqif = Ieqif*Ieqqth
    }
    
    csurv.indx = apply(cen_timeq >= t(array(rep(u, n), c(L.tpt, n))), 1, sum) #find index
    trt_index = csurv.indx[trt_id]
    con_index = csurv.indx[con_id]
    surv_trt = Kij_ori[trt_index]
    surv_con = Kij_ori[con_index]
    deltai = eventq[trt_id] #or trt_con[,1+ind.prior]
    deltaj = eventq[con_id] #or trt_con[,6+ind.prior]
    double_weight = surv_trt*surv_con
    double_weight[double_weight<t_val] = t_val    #truncation value
    IPCW_W[,q] = Ieqif*(deltai*deltaj)/(double_weight)
    G.i = G_ori[trt_index,]
    G.j = G_ori[con_index,]
    G.ij = as.matrix(G.i+G.j)
    epislon_W[,q] = as.numeric(G.ij%*%Omega_inv%*%matrix(mPsi,ncol=1))
  }
  weight_est = apply(IPCW_W*determined_flag,1,sum)
  weight_var0 = apply(epislon_W*determined_flag,1,sum)
  weight_var = weight_est*(1+weight_var0)
  res = list(weight_est = weight_est,weight_var=weight_var)
  return(res)
}
