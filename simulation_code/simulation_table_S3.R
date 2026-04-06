#This program can reproduce simulation results in Web Table S3 in Appendix
#by set different restriction event time in line 71 of the program.
rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)
library(survival)

setwd("C:/Users/82655/Dropbox/research/pim/R_code/simulation_code")
source("CreateScoreFun.R")
source("win.strategy.R")
source("merge.csd.R")
source("model.matrix.pim.csd.R")
source("censoring.weight.est.R")
source("pim.csd.fit.R")
source("ipcw.weight.m2.R")
source("pim.csd.m2.R")

#Generate data
gendata_full = function(n,beta,betac,lambda_h,lambda_d,alpha){
  z1 = rnorm(n)
  z1 = z1*(-1 <= z1 & z1 <= 1) - 1*(z1 < -1) + 1*(z1>1)
  z2 = 2*rbinom(n,1,0.5)-1
  outcome = rgumbel(n,alpha=alpha,dim=2)
  Lambda_D = as.numeric(lambda_d*exp(-cbind(z1,z2)%*%beta))
  Lambda_H = as.numeric(lambda_h*exp(-cbind(z1,z2)%*%beta))
  D = -log(outcome[,1])/Lambda_D #survival time
  H = -log(outcome[,2])/Lambda_H #nonfatal event time
  #x1 = runif(n)
  #x2 = rnorm(n,mean=0,sd=1.5)
  x1 = z1   #using the same covariates as PIM model
  x2 = z2
  C = -log(runif(n))/exp(cbind(x1,x2)%*%betac) #censoring time
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<C & H<Do] = 1;
  status2[H>pmin(Do,C)] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = cbind(id=1:n,time1= Do, time2 = Ho, delta1 = status1,
                   delta2 = status2, z1, z2, x1, x2)
  
  colnames(data_set) = c("id","y1","y2","delta1","delta2","z1","z2","x1","x2") 
  #note: y1 is death and y2 is hospitalization
  return(data_set)
}

#main simulation
N = 1000
n = 200  #sample size 
lambda_d = 0.2
lambda_h = 2
start = rep(0, 2) #start points
#you can set convergence condition in pim.co, or using default settings
xtol = 1e-6; ftol = 1e-6
btol = 1e-3; maxit = 50
cntl = list(xtol=xtol, ftol=ftol, btol=btol, maxit=maxit)
qna = qnorm(0.975)
betac = c(-1,0.5) #parameter coefficient in censoring model

case_set = 1:6;
L = 0.5  #restriction event
flag = flag1 = matrix(0,N,length(case_set))
for(case in case_set){
  if(case == 1){         #L=0.5,1,2 are about 50%,70-75%,90-95%
    alpha = 1 #or 2, 
    beta = c(-0.5,0.5) 
  }else if(case == 2){   #L=0.5,1,2 are about 50%,70-75%,85-90%
    alpha = 1 
    beta = c(0,0)     
  }else if(case == 3){   #L=0.5,1,2 are about 50%,65-70%,80-85%
    alpha = 1 
    beta = c(0.5,-0.5) 
  }else if(case == 4){  #L=0.5,1,2 are about 50%,70-75%,90-95%
    alpha = 2 
    beta = c(-0.5,0.5)
  }else if(case == 5){  #L=0.5,1,2 are about 50%,70-75%,85-90%
    alpha = 2 
    beta = c(0,0)
  }else if(case == 6){   #L=0.5,1,2 are about 50%,65-70%,80-85%
    alpha = 2 
    beta = c(0.5,-0.5)
  }
  beta_est = se_est = cova = matrix(0, N, 2)
  beta_est1 = se_est1 = cova1 =matrix(0, N, 2)
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(123456+i)
    #generate data
    gdata = gendata_full(n,beta,betac,lambda_h,lambda_d,alpha)
    #logit link
    res = pim.csd.m2(formula=y1~z1+z2,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("x1","x2"),
                  L = L, adjusted="unadjusted", priority = 1:2, model = "difference",start = start,
                  link = "logit", method = "Newton")
    flag[i,case] = res$flag
    #probit link
    res1 = pim.csd.m2(formula=y1~z1+z2,data=gdata,SurvTime="y1",Status="delta1",covariateC=c("x1","x2"),
                   L = L, adjusted="unadjusted",priority = 1:2, model = "difference",start = start,
                   link = "probit", method = "Newton")
    
    flag1[i,case] = res1$flag
    beta_est[i,] = res$coef
    vare = res$vcov
    se_est[i,] = sqrt(diag(vare))
    low = beta_est[i,]-qna*se_est[i,]
    high = beta_est[i,]+qna*se_est[i,]
    if(low[1] <= beta[1] & beta[1] <= high[1]) cova[i,1] = 1 #modify
    if(low[2] <= beta[2] & beta[2] <= high[2]) cova[i,2] = 1 #modify
    #transform back by using delta method  
    beta_a = res1$coef
    pbeta = pnorm(beta_a)
    dbeta = dnorm(beta_a) 
    beta_est1[i,] = log(pbeta/(1-pbeta))
    vare_a = res1$vcov
    #y=ln(Phi(x)/(1-Phi(x)))--->y'=phi(x)/(Phi(x)*(1-Phi(x)))
    se_est1[i,] = abs(dbeta/(pbeta*(1-pbeta)))*sqrt(diag(vare_a)) #delta method
    low1 = beta_est1[i,]-qna*se_est1[i,]
    high1 = beta_est1[i,]+qna*se_est1[i,]
    if(low1[1] <= beta[1] & beta[1] <= high1[1]) cova1[i,1] = 1 #modify
    if(low1[2] <= beta[2] & beta[2] <= high1[2]) cova1[i,2] = 1 #modify
  }
  #summary estimation results
  index = flag[,case]*flag1[,case] 
  est0 = colMeans(beta_est[index==1,])
  est1 = colMeans(beta_est1[index==1,])
  
  sd0 = apply(beta_est[index==1,],2,sd)
  sd1 = apply(beta_est1[index==1,],2,sd)
  se0 = colMeans(se_est[index==1,])
  se1 = colMeans(se_est1[index==1,])
  cp0 = colMeans(cova[index==1,])
  cp1 = colMeans(cova1[index==1,])
  resu = data.frame(est=est0,sd=sd0,se=se0,cp=cp0,
                    est1=est1,sd1=sd1,se1=se1,cp1=cp1)
  res = round(resu,3)
  write.csv(res,paste0("pim_logit_L",L,"_case_",case,"_n_",n,".csv"))
  cat("case=",case,"\n")
}

