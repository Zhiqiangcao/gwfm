ipcw.weight.m1 = function(tdata,trt_con,censoring_weight,win_status,priority,determined_flag,L,n_ep){
  #first, we compute double weights for estimation
  n = dim(tdata)[1] 
  # Setup survival object (based on the Death time)
  y1 = tdata[,2]
  u = c(0,y1)
  #obtain survival estimation at each observed time, from the smallest to the largest
  Kij = diag(censoring_weight$Kij)
  Kij_ori = c(1,Kij)
  L.tpt = length(u)
  npim = dim(win_status)[1]
  #the following is used to compute weight_var
  Omega = censoring_weight$Omega
  Omega_inv = ginv(Omega)  #generalized inverse
  U = censoring_weight$U
  mU = apply(U,2,mean)
  pc = length(mU) #number of covariates in censoring model
  G = censoring_weight$G  #is a n*pc matrix
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
    if(q == 1){ #no matter truncated or not
      cen_timeq = tdata[,(1+ind.prior)]
    }else{ #q>=2
      trt_centimq_1 = trt_con[,1+n_ep+ind.prior-1] #observed survival time of the q-1 th endpoint
      con_centimq_1 = trt_con[,2+3*n_ep+ind.prior-1]
      if(L == Inf) { #no truncated
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
    surv_trt = Kij_ori[con_index] #based on control group
    surv_con = Kij_ori[con_index] #based on control group
    deltai = eventq[con_id] 
    deltaj = eventq[con_id] 
    double_weight = surv_trt*surv_con
    double_weight[double_weight<0.01] = 0.01    #truncated 
    IPCW_W[,q] = Ieqif*(deltai*deltaj)/(double_weight)
    G.i = G_ori[trt_index,]
    G.j = G_ori[con_index,]
    G.ij = as.matrix(G.i+G.j)
    epislon_W[,q] = as.numeric(G.ij%*%Omega_inv%*%matrix(mU,ncol=1))
  }
  weight_est = apply(IPCW_W*determined_flag,1,sum)
  #weight_est[ties_flag==0] = 1 #let weight for tie situation be zero
  weight_var0 = apply(epislon_W*determined_flag,1,sum)
  #since epsilon corresponds residuals of \hat{W}_{ij}-W_{ij},
  #so for ties comparison, we let weight_var0 = 0
  weight_var = weight_est*(1+weight_var0)
  res = list(weight_est = weight_est,weight_var=weight_var)
  return(res)
}
