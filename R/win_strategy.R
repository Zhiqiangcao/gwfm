#' The win strategy function
#' 
#' An function to determine the win status for each pair based on the win strategy of Dong et al.(2022).
#' Specifically, one compares each subject in the "treatment" group with every subject in the "control"
#' group to determine the win status. 
#' 
#' @param trt_con Given n subjects, there are n(n-1)/2 pairs for non-marginal structure and n(n-1) pairs for marginal structure, each row represents a pair. The analysis dataset trt_con contains the following variables:
#' \describe{
#'  \item{pid_trt}{A vector for the subject id of the individuals from the "treatment" group}
#'  \item{Delta_j_trt}{A vector for the event status of the jth endpoint (1=event, 0=censored) for the individuals from the "treatment" group}
#'  \item{Y_j_trt}{A vector for the outcome of the jth endpoint for the individuals from the "treatment" group}
#'  \item{pid_con}{A vector for the subject id of the individuals from the "control" group}
#'  \item{Delta_j_con}{A vector for the event status of the jth endpoint (1=event, 0=censored) for the individuals from the "control" group}
#'  \item{Y_j_con}{A vector for the outcome of the jth endpoint for the individuals from the "control" group}
#' }
#' @param priority Priority order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 1, Endpoint 2, and Endpoint 3, input priority = c(1,2,3).
#'       
#' @return A data.frame for the win status of each pair for each endpoint.
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' 
#' @references Cui, Y. and Huang, B. Package "WINS", CRAN R Repositary, 2025.
#' 
#' Dong, G., Huang, B., Wang, D., Verbeeck, J., Wang, J. and Hoaglin, D. Adjusting
#' win statistics for dependent censoring. Pharm Stat, 2021, 20(3):440--450.

win_strategy <- function(trt_con, priority){
  n_ep = length(priority)
  #### Obtain the indicator of the first endpoint for treatment and control
  colname.trt_con = colnames(trt_con)
  
  ind.delta1.trt = which(colname.trt_con=="Delta_1_trt")
  ind.delta1.con = which(colname.trt_con=="Delta_1_con")
  
  ind.time1.trt = which(colname.trt_con=="Y_1_trt")
  ind.time1.con = which(colname.trt_con=="Y_1_con")
  
  
  win.status0 = NULL
  
  #### For time-to-event outcome: Denote the observed survival time as Y_trt and Y_con, and event status
  #### as Delta_trt and Delta_con. There is a win for the treatment group if we have:
  #### Delta_con = 1 and Y_trt > Y_con,
 
  for(l in priority){
    delta_l_trt = trt_con[,ind.delta1.trt+l-1]; delta_l_con = trt_con[,ind.delta1.con+l-1]
    Y_l_trt = trt_con[,ind.time1.trt+l-1]; Y_l_con = trt_con[,ind.time1.con+l-1]
    
    win.temp1 = ((delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt > Y_l_con) |
                   (delta_l_trt < 1 & delta_l_con == 1 & Y_l_trt > Y_l_con))
    win.temp2 = ((delta_l_trt == 1 & delta_l_con == 1 & Y_l_con > Y_l_trt) |
                   (delta_l_trt == 1 & delta_l_con < 1 & Y_l_con > Y_l_trt))
    
    win.status0 = cbind(win.status0, win.temp1, win.temp2)
  }
  
  #### prioritize: once a winner is determined, then all the subsequent is set to zero
  win_status = t(apply(win.status0, 1, func<-function(x){
    if(sum(x)>1){
      temp = x; temp[min((2*min(ceiling(which(x==1)/2))+1),(2*n_ep-1)):(2*n_ep)] = 0
      return(temp)
    }else{
      return(as.numeric(x))
    }
  }))
  
  colnames(win_status) = paste0(rep(c("Trt","Con"),n_ep),"_Endpoint",rep(priority,each=2))
  win_status = as.data.frame(win_status)
  win_status = data.frame(pid_trt = trt_con[,1], pid_con = trt_con[,2+2*n_ep], win_status)
  return(win_status)
}
