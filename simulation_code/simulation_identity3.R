#calculate true value of pim based on identity link with truncated time
#calculation true values based on identity link
rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)
library(survival)
library(mvtnorm)

setwd("C:/Users/82655/Dropbox/research/pim/R_code/simulation_paper")
source("CreateScoreFun.R")
source("win.strategy.R")
source("merge.csd.R")
source("model.matrix.pim.csd.R")
source("censoring.weight.est.R")
source("pim.csd.fit.R")
source("ipcw.weight.m2.R")
source("pim.csd.m2.R")
source("pim.csd.m2fff.R")

gendata = function(n,beta_d,beta_h,lambda_d,lambda_h,betac){
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
  C = -log(runif(n))/exp(z*betac) #censoring time
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
                   delta2 = status2,x,z)
  colnames(data_set) = c("id","y1","y2","delta1","delta2","x","z") 
  #note: y1 is death and y2 is hospitalization
  return(data_set)
}


###censoring rate chosen functions
censor_choice = function(choice,case){
  if(choice==1){ #20% censoring
    if(case == 1){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=11.5; beta=0.7348549 #about 50% quantile of D
    }else if(case == 2){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.3; L=11.5; beta=0.7422573 #about 50% quantile of D
    }else if(case == 3){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=11.5; beta = 0.7201865  #about 50% quantile of D
    }else if(case == 4){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=15; beta=0.7049795  #about 70% quantile of D
    }else if(case == 5){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.3; L=15; beta = 0.7208126  #about 70% quantile of D
    }else if(case == 6){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=15; beta = 0.7023488  #about 70% quantile of D
    }else if(case == 7){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=26; beta= 0.6951691
    }else if(case == 8){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.3; L=22; beta = 0.7170939
    }else if(case == 9){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=27; beta = 0.6950332
    }else if(case == 10){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.7; L=Inf; beta= 0.6947183 
    }else if(case == 11){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-5.3; L=Inf; beta = 0.7167179   
    }else if(case == 12){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-5.7; L=Inf; beta = 0.6947183   
    }
  }else{   #40% censoring
    if(case == 1){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.2; L=11.5; beta=0.7348549 #about 50% quantile of D
    }else if(case == 2){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4; L=11.5; beta=0.7422573 #about 50% quantile of D
    }else if(case == 3){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-4.2; L=11.5; beta = 0.7201865  #about 50% quantile of D
    }else if(case == 4){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.2; L=15; beta=0.7049795  #about 70% quantile of D
    }else if(case == 5){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4; L=15; beta = 0.7208126  #about 70% quantile of D
    }else if(case == 6){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-4.2; L=15; beta = 0.7023488  #about 70% quantile of D
    }else if(case == 7){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.2; L=26; beta= 0.6951691
    }else if(case == 8){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4; L=22; beta = 0.7170939
    }else if(case == 9){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-4.2; L=27; beta = 0.6950332
    }else if(case == 10){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4.2; L=Inf; beta= 0.6947183 #about 70% quantile of D
    }else if(case == 11){
      lambda_d = 0.2; lambda_h = 0.3; beta_d = 10; beta_h = 8; betac=-4; L=Inf; beta = 0.7167179   #about 70% quantile of D
    }else if(case == 12){
      lambda_d = 0.15; lambda_h = 0.3; beta_d = 10; beta_h = 6; betac=-4.2; L=Inf; beta = 0.6947183   #about 70% quantile of D
    }
  }
  
  res = list(lambda_d=lambda_d, lambda_h = lambda_h, L=L,
             beta_d = beta_d, beta_h = beta_h, betac=betac, beta = beta)
  return(res)
}
#                    cen_rate hosp_rate              cen_rate hosp_rate
#case1: betac=-5.7， 0.201860 0.724189； betac=-3.8  0.493797 0.550627 
#case2: betac=-5.3， 0.212009 0.669496； betac=-3.5  0.501573 0.474838 
#case3: betac=-5.7， 0.203514 0.769667； betac=-3.8  0.506140 0.597237
#case1: betac=-4.2， 0.404088 0.596327 
#case2: betac=-4.0， 0.403926 0.547022 
#case3: betac=-4.2， 0.407885 0.617966 

##parameter settings
n = 200 #sample size 200 or 500
N = 1000 #simulation times

#you can set convergence condidition in pim.co, or using default settings
xtol = 1e-6; ftol = 1e-6
btol = 1e-3; maxit = 50
cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)

case_set = 1:12
qna = qnorm(0.975)
choice = 1
for(case in case_set){
  cc_res = censor_choice(choice,case)
  lambda_d = cc_res$lambda_d
  lambda_h = cc_res$lambda_h
  beta_d = cc_res$beta_d
  beta_h = cc_res$beta_h
  #true coefficients
  L = cc_res$L
  betac = cc_res$betac
  beta = cc_res$beta
  flag1 = flag2 = rep(0,N)
  beta_est = se_est = cova = matrix(0, N, 2)
  kk = i = 1
  while(i<=N){
    cat("iter=",i,"\n")
    set.seed(10000+kk)
    gdata = gendata(n,beta_d,beta_h,lambda_d,lambda_h,betac)
    mint = min(gdata[,2],gdata[,3])
    if(mint<0){ #to make sure generated death time and hospitalization time is >0
      i=i; kk=kk+1 
    }else{
      #method 1
      res1 = pim.csd.m2(formula=y1~x,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("z"),
                        priority = 1:2, L = L, adjusted = "unadjusted", model = "marginal",
                        link = "identity",method = "Newton",control = cntl)
      flag1[i] = res1$flag
      #method 2
      res2 = pim.csd.m2(formula=y1~x,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("z"),
                        priority = 1:2, L = L,adjusted = "IPCW", model = "marginal",
                        link = "identity",method = "Newton",control = cntl)
      flag2[i] = res2$flag
      
      vare1 = res1$vcov; vare2 = res2$vcov;  
      flagg1 = as.numeric(is.na(vare1)|is.na(vare2))
      flagg2 = vare1<0|vare2<0
      flagg3 = is.infinite(vare1)|is.infinite(vare2)
      flagg4 = flag1[i]!=1 | flag2[i]!=1 
      if(!flagg1 & !flagg2 & !flagg3 & !flagg4){
        #summary results
        #method 1
        beta_est[i,1] = res1$coef
        se_est[i,1] = sqrt(diag(vare1))
        low1 = beta_est[i,1]-qna*se_est[i,1]
        high1 = beta_est[i,1]+qna*se_est[i,1]
        if(low1 <= beta & beta <= high1) cova[i,1] = 1 
        #method 2
        beta_est[i,2] = res2$coef
        se_est[i,2] = sqrt(diag(vare2))
        low2 = beta_est[i,2]-qna*se_est[i,2]
        high2 = beta_est[i,2]+qna*se_est[i,2]
        if(low2 <= beta & beta <= high2) cova[i,2] = 1
        i=i+1; kk=kk+1 
      }else{
        i=i; kk=kk+1
      }
    }
  }
  #summary estimation results
  index = flag1*flag2
  #estimatiion
  est1 = colMeans(beta_est[index==1,])
  #SD
  sd1 = apply(beta_est[index==1,],2,sd)
  #SE
  se1 = colMeans(se_est[index==1,])
  #CP
  cp1 = colMeans(cova[index==1,])
  resu = data.frame(est=est1,sd=sd1,se=se1,cp=cp1)
  res = round(resu,3)
  write.csv(res,paste0("pim_identity3_choice_",choice,"_case_",case,"_n_",n,".csv"))
  cat("case=",case,"\n")
}


