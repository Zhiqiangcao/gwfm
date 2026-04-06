#This function shows how to calculate true values of \beta_L for Table S6 in Web Appendix
rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)
library(survival)
library(mvtnorm)


setwd("C:/Users/82655/Dropbox/research/pim/R_code/simulation_code/R_programs")
source("CreateScoreFun.R")
source("win.strategy.R")
source("merge.csd.R")
source("model.matrix.pim.csd.R")


##generate data
gendata = function(n,alpha,u,sig,betac,C){
  x = seq(2,u,length=n)
  y = alpha*x
  mv.NE = mvdc(normalCopula(0.75), c("norm", "norm"),
               list(list(mean = 0, sd =sig), list(mean = 0, sd = sig)))
  x.sample = rMvdc(n, mv.NE)  ##error term
  D = y + x.sample[,1]
  H = y + x.sample[,2]
  #C = -log(runif(n))/exp(x*betac) #censoring time
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
                   delta2 = status2,x)
  colnames(data_set) = c("id","y1","y2","delta1","delta2","x") 
  #note: y1 is death and y2 is hospitalization
  return(data_set)
}


###censoring rate chosen functions
censor_choice = function(case){
  if(case == 1){
    alpha = 1; u = 5; sig = 0.4; L=4; betac = -0.79;
  }else if(case == 2){
    alpha = 1; u = 5; sig = 0.4; L=4.7; betac = -0.79;
  }else if(case == 3){
    alpha = 1; u = 5; sig = 0.4; L=5; betac = -0.79; 
  }else if(case == 4){
    alpha = 1; u = 5; sig = 0.4; L=Inf; betac = -0.79; 
  }else if(case == 5){
    alpha = 1; u = 10; sig = 0.8; L=8;  betac = -0.58; 
  }else if(case == 6){
    alpha = 1; u = 10; sig = 0.8; L=9.2; betac = -0.58;
  }else if(case == 7){
    alpha = 1; u = 10; sig = 0.8; L=9.8;  betac = -0.58; 
  }else if(case == 8){
    alpha = 1; u = 10; sig = 0.8; L=Inf; betac = -0.58; 
  }else if(case == 9){
    alpha = 2; u = 10; sig = 1.5; L=16; betac = -0.74; 
  }else if(case == 10){
    alpha = 2; u = 10; sig = 1.5; L=18.5; betac = -0.74;
  }else if(case == 11){
    alpha = 2; u = 10; sig = 1.5; L=19.5; betac = -0.74;
  }else if(case == 12){
    alpha = 2; u = 10; sig = 1.5; L=Inf; betac = -0.74; 
  }
  res = list(alpha=alpha, u=u, sig=sig, betac=betac, L=L)
  return(res)
}

#
#main simulation
#general settings
n = 1000 #sample size 200 or 500
N = 1000 #simulation times
link = "probit"
model = "difference"
na.action = getOption("na.action")
formula = y1~x
compare = if (model == "marginal") "all" else "unique"
priority = 1:2
C = 10^4
flag = matrix(0,N,18)
true_val = NULL
case_set = 1:12
for(case in case_set){
  cc_res = censor_choice(case)
  alpha = cc_res$alpha
  u = cc_res$u
  sig = cc_res$sig
  betac = cc_res$betac
  L = cc_res$L
  est_beta = rep(0,N)
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(123456+i)
    gdata = gendata(n,alpha,u,sig,betac,C)
    data = gdata
    mint = min(gdata[,2],gdata[,3])
    if(mint<0) flag[i,case] = 1
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
    #pid_trt Delta_1_trt Delta_2_trt Y_1_trt Y_2_trt pid_con Delta_1_con Delta_2_con Y_1_con  Y_2_con
    #      17           0           0    1.81    1.81      32           1           1    2.14 1.77
    #      18           0           0    1.84    1.84      32           1           1    2.14 1.77
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
    #probit link
    lm1 = glm(y~x-1,family = binomial(link = "probit"))
    est_beta[i] = coef(lm1)
  }
  truebeta = mean(est_beta[flag[,case]==0])
  true_val = c(true_val,truebeta)
}
colSums(flag)
true_val
#Note, for L=Inf, we can apply the relationship between PIM and normal linear regression
#to obtain true value of beta_L, that is,
#beta = alpha/(sqrt(2)*sig) 

##determine different censoring rate
cal_rate = function(n,alpha,u,sig,betac){
  x = seq(2,u,length=n)
  y = alpha*x
  mv.NE = mvdc(normalCopula(0.75), c("norm", "norm"),
               list(list(mean = 0, sd =sig), list(mean = 0, sd = sig)))
  x.sample = rMvdc(n, mv.NE)  ##error term
  D = y + x.sample[,1]
  H = y + x.sample[,2]
  m = 1000
  cenr1 = cenr2 = rep(0,m)
  for(i in 1:m){
    C = -log(runif(n))/exp(x*betac) #censoring time
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


n=1000
case=6
if(case == 1){
  alpha = 1; u = 5; sig = 0.4; L=4; betac = -0.79
}else if(case == 2){
  alpha = 1; u = 10; sig = 0.8; L=8; betac = -0.58
}else if(case == 3){
  alpha = 2; u = 10; sig = 1.5; L=16; betac = -0.74
}else if(case == 4){
  alpha = 1; u = 5; sig = 0.4; L=4; betac = -0.54  
}else if(case == 5){
  alpha = 1; u = 10; sig = 0.8; L=8; betac = -0.4 
}else if(case == 6){
  alpha = 2; u = 10; sig = 1.5; L=16; betac = -0.53  #90% quantile of D
}
cal_rate(n,alpha,u,sig,betac)