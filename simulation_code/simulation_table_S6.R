#This program can reproduce simulation results in Web Table S6 in Appendix
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
source("censoring.weight.est.R")
source("pim.csd.fit.R")
source("ipcw.weight.m2.R")
source("pim.csd.m2.R")
source("pim.csd.m2fff.R")

#generate time
gendata = function(n,alpha,u,sig,betac){
  x = seq(2,u,length=n)
  y = alpha*x
  mv.NE = mvdc(normalCopula(0.75), c("norm", "norm"),
               list(list(mean = 0, sd =sig), list(mean = 0, sd = sig)))
  x.sample = rMvdc(n, mv.NE)  ##error term
  D = y + x.sample[,1]
  H = y + x.sample[,2]
  C = -log(runif(n))/exp(x*betac) #censoring time
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


censor_choice = function(choice,case){
  if(choice==1){ #20% censoring
    if(case == 1){
      alpha = 1; u = 5; sig = 0.4; L=4; betac = -0.79; beta = 1.8364276
    }else if(case == 2){
      alpha = 1; u = 5; sig = 0.4; L=4.7; betac = -0.79; beta = 1.775185  #90% quantile of D
    }else if(case == 3){
      alpha = 1; u = 5; sig = 0.4; L=5; betac = -0.79; beta=1.7697528
    }else if(case == 4){
      alpha = 1; u = 5; sig = 0.4; L=Inf; betac = -0.79; 
    }else if(case == 5){
      alpha = 1; u = 10; sig = 0.8; L=8;  betac = -0.58; beta = 0.9056841
    }else if(case == 6){
      alpha = 1; u = 10; sig = 0.8; L=9.2; betac = -0.58; beta = 0.8878915   #90% quantile of D
    }else if(case == 7){
      alpha = 1; u = 10; sig = 0.8; L=9.8;  betac = -0.58; beta=0.8849669 
    }else if(case == 8){
      alpha = 1; u = 10; sig = 0.8; L=Inf; betac = -0.58; 
    }else if(case == 9){
      alpha = 2; u = 10; sig = 1.5; L=16; betac = -0.74; beta = 0.9664529
    }else if(case == 10){
      alpha = 2; u = 10; sig = 1.5; L=18.5; betac = -0.74; beta = 0.9465471  #90% quantile of D
    }else if(case == 11){
      alpha = 2; u = 10; sig = 1.5; L=19.5; betac = -0.74; beta=0.9438353
    }else if(case == 12){
      alpha = 2; u = 10; sig = 1.5; L=Inf; betac = -0.74; 
    }
  }else{ #40% censoring
    if(case == 1){
      alpha = 1; u = 5; sig = 0.4; L=4; betac = -0.54; beta = 1.8364276
    }else if(case == 2){
      alpha = 1; u = 5; sig = 0.4; L=4.7; betac = -0.54; beta = 1.775185  #90% quantile of D
    }else if(case == 3){
      alpha = 1; u = 5; sig = 0.4; L=5; betac = -0.54; beta=1.7697528
    }else if(case == 4){
      alpha = 1; u = 5; sig = 0.4; L=Inf; betac = -0.54; 
    }else if(case == 5){
      alpha = 1; u = 10; sig = 0.8; L=8;  betac = -0.4; beta = 0.9056841
    }else if(case == 6){
      alpha = 1; u = 10; sig = 0.8; L=9.2; betac = -0.4; beta = 0.8878915   #90% quantile of D
    }else if(case == 7){
      alpha = 1; u = 10; sig = 0.8; L=9.8;  betac = -0.4; beta=0.8849669 
    }else if(case == 8){
      alpha = 1; u = 10; sig = 0.8; L=Inf; betac = -0.4; 
    }else if(case == 9){
      alpha = 2; u = 10; sig = 1.5; L=16; betac = -0.53; beta = 0.9664529
    }else if(case == 10){
      alpha = 2; u = 10; sig = 1.5; L=18.5; betac = -0.53; beta = 0.9465471  #90% quantile of D
    }else if(case == 11){
      alpha = 2; u = 10; sig = 1.5; L=19.5; betac = -0.53; beta=0.9438353
    }else if(case == 12){
      alpha = 2; u = 10; sig = 1.5; L=Inf; betac = -0.53; 
    }
  }
  res = list(alpha=alpha, u=u, sig=sig, betac=betac, L=L, beta = beta)
  return(res)
}


##parameter settings
n = 200 #sample size 200 or 500
N = 1000 #simulation times

#you can set convergence condidition in pim.co, or using default settings
xtol = 1e-6; ftol = 1e-6
btol = 1e-3; maxit = 50
cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)

case_set = 1:12;
choice = 1
qna = qnorm(0.975)
for(case in case_set){
  cc_res = censor_choice(choice,case)
  alpha = cc_res$alpha
  u = cc_res$u
  sig = cc_res$sig
  betac = cc_res$betac
  #true coefficients
  L = cc_res$L
  if(L == Inf){
    beta = alpha/(sqrt(2)*sig) #relationship between alpha and beta
  }else{
    beta = cc_res$beta
  }
  flag1 = flag2 = rep(0,N)
  beta_est = se_est = cova = matrix(0, N, 2)
  kk = i = 1
  while(i<=N){
    cat("iter=",i,"\n")
    set.seed(10000+kk)
    gdata = gendata(n,alpha,u,sig,betac)
    mint = min(gdata[,2],gdata[,3])
    if(mint<0){ #to make sure generated death time and hospitalization time is >0
      i=i; kk=kk+1 
    }else{
      #method 1
      res1 = pim.csd.m2(formula=y1~x,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("x"),
                        priority = 1:2, L = L, adjusted = "unadjusted", model = "difference",
                        link = "probit",method = "Newton",control = cntl)
      flag1[i] = res1$flag
      #method 2
      res2 = pim.csd.m2(formula=y1~x,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("x"),
                        priority = 1:2, L = L,adjusted = "IPCW", model = "difference",
                        link = "probit",method = "Newton",control = cntl)
      flag2[i] = res2$flag
      
      vare1 = res1$vcov; vare2 = res2$vcov;  #vare3 = res3$vcov;  
      flagg1 = as.numeric(is.na(vare1)|is.na(vare2))  #|is.na(vare3))
      flagg2 = vare1<0|vare2<0  #|vare3<0
      flagg3 = is.infinite(vare1)|is.infinite(vare2)  #|is.infinite(vare3)
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
  index = flag1*flag2#*flag3
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
  write.csv(res,paste0("tables6_choice_",choice,"_case_",case,"_n_",n,".csv"))
  cat("case=",case,"\n")
}
