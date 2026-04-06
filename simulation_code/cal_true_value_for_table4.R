#This function shows how to calculate true values of \beta_L for Table 4 of main text
rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)
library(survival)
library(mvtnorm)

#You should reset the path below to run the code
setwd("E:/")
source("CreateScoreFun.R")
source("win.strategy.R")
source("merge.csd.R")
source("model.matrix.pim.csd.R")
source("censoring.weight.est.R")



gendata = function(n,beta_d,beta_h,lambda_d,lambda_h,C){
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
  data_set = cbind(id=1:n,time1= Do, time2 = Ho, delta1 = status1,
                   delta2 = status2, x)
  colnames(data_set) = c("id","y1","y2","delta1","delta2","x")
  #note: y1 is death and y2 is hospitalization
  return(data_set)
}


###censoring rate chosen functions
censor_choice = function(case){
  if(case == 1){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=11.5 #about 50-55% quantile of D
  }else if(case == 2){
    lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=11.5 #about 55-65% quantile of D
  }else if(case == 3){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=11.5 #about 50-55% quantile of D
  }else if(case == 4){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=15  #about 70-75% quantile of D
  }else if(case == 5){
    lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=15  #about 75-80% quantile of D
  }else if(case == 6){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=15  #about 70-75% quantile of D
  }else if(case == 7){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=26  #about 95% quantile of D
  }else if(case == 8){
    lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=22  #about 95% quantile of D
  }else if(case == 9){
    lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=27  #about 95% quantile of D
  }
  res = list(lambda_d=lambda_d, lambda_h = lambda_h, L=L,
             beta_d = beta_d, beta_h = beta_h, betac=betac)
  return(res)
}


#main simulation
#general settings
n = 500  #to save time, here we use n=500
N = 1000 #simulation times
link = "identity"
model = "marginal"
na.action = getOption("na.action")
formula = y1~x
compare = if (model == "marginal") "all" else "unique"
priority = 1:2
C = 10^4
true_val = NULL
case_set = 1:9
for(case in case_set){
  cc_res = censor_choice(case)
  lambda_d = cc_res$lambda_d
  lambda_h = cc_res$lambda_h
  beta_d = cc_res$beta_d
  beta_h = cc_res$beta_h
  L = cc_res$L
  est_beta = rep(0,N)
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(10000+i)
    gdata = gendata(n,beta_d,beta_h,lambda_d,lambda_h,C)
    data = gdata
    
    #create x
    if (is.null(na.action)) 
      na.action <- "na.fail"
    if (!is.character(na.action)) 
      na.action <- deparse(substitute(na.action))
    
    f.terms = terms(formula, simplify = TRUE)
    vars = all.vars(formula)
    
    penv = new.pim.env(data, compare = compare, vars = vars, env = parent.frame())
    ff = new.pim.formula(formula, penv)
    x = model.matrix.pim.csd(ff, na.action = na.action, model = model)
    
    # next, we prepare comparison results between composite survival endpoints 
    # length of endpoints
    n_ep = length(priority)
    # since the data set format is id, y1, ..., y_nep, delta1, ..., delta_nep, z1,...,zp
    datap = as.matrix(data[,c(1:(1+2*n_ep))])
    # construct truncated event times and data if L != Inf
    y_o = as.matrix(data[,2:(n_ep+1)]) #observed endpoints
    delta_o = as.matrix(data[,(n_ep+2):(1+2*n_ep)]) #observed censoring indicators
    y_n = pmin(y_o,L) #truncated event times
    delta_n = ifelse(y_n == L, 1, delta_o)
    datap_n = as.matrix(cbind(data[,1], y_n, delta_n))
    #merge data similar to compare treatment and control groups
    trt_con = merge.csd(datap = datap_n, model = model, n_ep = n_ep)
    #next, we compute pairwise comparison result based on Dong et al's (2020) method
    win_status = win.strategy(trt_con = trt_con, priority = priority)
    #find out win function is determined by which endpoint comparison
    determined_flag = NULL
    for(q in seq(3,(2+2*n_ep-1),2)){
      win_qth_res = win_status[,q:(q+1)]
      qth_flag = apply(as.matrix(win_qth_res,ncol = 2),1,sum)
      determined_flag = cbind(determined_flag,qth_flag)
    }
    #find out index of ties comparison (i.e., ties_flag=0)
    flag2 = which(determined_flag[,2]==1)  #find out index determined by the comparison of the 2th endpoint
    index = rep(1,length(flag2))
    index[trt_con[flag2,4]==trt_con[flag2,5] & trt_con[flag2,5]!=L] = 0 #index=0 means the above tie situation
    index[trt_con[flag2,9]==trt_con[flag2,10] & trt_con[flag2,10]!=L] = 0
    determined_flag[,2][flag2] = index  #set the above situation to ties
    ties_flag = apply(determined_flag,1,sum)
    #y is the win result of W(Y_i,Y_j)
    y = apply(as.matrix(win_status[,seq(3,(2+2*n_ep-1),2)],ncol = n_ep),1,sum)
    #make those ties comparison results do not contribute in score function
    x = as.numeric(ties_flag*x)
    lmreg = lm(y~x-1)
    est_beta[i] = coef(lmreg)
  }
  truebeta = mean(est_beta)
  true_val = c(true_val,truebeta)
}
true_val
#results: for cases 1-3
#L=11.5 0.7348549 0.7422573 0.7201865
#results: for cases 4-6
#L=15   0.7049795 0.7208126 0.7023488
#results: for cases 7-9
#L=c(26,22,27) 0.6951691 0.7170939 0.6950332

##determine different censoring rate
cal_rate = function(n,beta_d,beta_h,lambda_d,lambda_h,betac){
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
  m = 1000
  cenr1 = cenr2 = rep(0,m)
  for(i in 1:m){
    C = -log(runif(n))/exp(z*betac) #censoring time
    status1 = 1*(D<=C)
    Do = pmin(C,D)
    status2 = rep(0,n)
    #if simulated HFH time is after death then 
    #we sill not observe the HFH. So we censor at 
    #the earliest of C years or time of death
    status2[H<=C & H<=Do] = 1
    status2[H>Do] = 0
    cenr1[i] = sum(status1==0)/n #censoring rate (1-death rate)
    cenr2[i] = mean(status2)     #nonfata-rvent rates
  }
  res = apply(cbind(cenr1,cenr2),2,mean)
  return(res)
}

#parameter settings
if(case == 1){
  lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.2; L=11.5 #about 50% quantile of D
}else if(case == 2){
  lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.0; L=11.5 #about 50% quantile of D
}else if(case == 3){
  lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-4.2; L=11.5 #about 50% quantile of D
}

cal_rate(n=1000,beta_d,beta_h,lambda_d,lambda_h,betac)
#case1: betac=-5.7， 0.201860 0.724189； 
#case2: betac=-5.7， 0.212009 0.669496； 
#case3: betac=-5.7， 0.203514 0.769667； 

#case1: betac=-4.2， 0.404088 0.596327 
#case2: betac=-4.2， 0.403926 0.547022 
#case3: betac=-4.2， 0.407885 0.617966 