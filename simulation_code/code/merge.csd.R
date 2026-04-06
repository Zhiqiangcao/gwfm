#merge treatment and control groups data, only tidy y and corresponding delta
merge.csd <- function(datap,model,n_ep){
  n = nrow(datap)
  contrast.ij <- function(x){
    di = datap[x[1],] #data of ith individual
    dj = datap[x[2],] #data of jth individual
    
    yi = di[2:(1+n_ep)] #y_1,y_2,...,y_ep
    deltai = di[(2+n_ep):(1+2*n_ep)] #delta_1,delta_2,...,delta_ep
    
    yj = dj[2:(1+n_ep)] #y_1,y_2,...,y_ep
    deltaj = dj[(2+n_ep):(1+2*n_ep)] #delta_1,delta_2,...,delta_ep
    return(c(deltai,yi,deltaj,yj))
  }
  if(model != "marginal"){
    value = combn(1:n, 2, FUN = contrast.ij)
    index = combn(1:n, 2)
  }else{
    index_m = rbind(rep(1,n-1),2:n)
    for(i in 2:(n-1)){
      index_i = rbind(rep(i,n-1),c(1:(i-1),(i+1):n))
      index_m = cbind(index_m,index_i)
    }
    index_n = rbind(rep(n,n-1),1:(n-1))
    index = cbind(index_m,index_n)
    value = apply(index, 2, contrast.ij)
  }
  outcome = t(rbind(index, value))
  outcome = outcome[,c(1,3:(2+2*n_ep),2,(3+2*n_ep):(2+4*n_ep))]
  delta_trt = paste0("Delta_",1:n_ep,"_trt",sep ="")
  colnames(outcome) = c("pid_trt",paste0("Delta_",1:n_ep,"_trt",sep =""),
                        paste0("Y_",1:n_ep,"_trt",sep =""),
                        "pid_con",paste0("Delta_",1:n_ep,"_con",sep =""),
                        paste0("Y_",1:n_ep,"_con",sep =""))
  return(outcome)
}
