#This program can reproduce simulation results in Web Table S2 in Appendix
rm(list=ls())
library(gumbel)
library(WR)
library(MASS)
library(pim)
library(copula)
library(nleqslv)
library(survival)

setwd("C:/Users/82655/Dropbox/research/pim/R_code/simulation_code/R_programs")
source("CreateScoreFun.R")
source("win.strategy.R")
source("merge.csd.R")
source("model.matrix.pim.csd.R")
source("censoring.weight.est.R")
source("pim.csd.fit.R")
source("ipcw.weight.m2.R")
source("pim.csd.m2.R")

#Generate data for our setup
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

##generate data for Mao and Wang's format
gendata = function(n,beta,betac,lambda_h,lambda_d,alpha){
  z1 = rnorm(n)
  z1 = z1*(-1 <= z1 & z1 <= 1) - 1*(z1 < -1) + 1*(z1>1)
  z2 = 2*rbinom(n,1,0.5)-1
  outcome = rgumbel(n,alpha=alpha,dim=2)
  Lambda_D = as.numeric(lambda_d*exp(-cbind(z1,z2)%*%beta))
  Lambda_H = as.numeric(lambda_h*exp(-cbind(z1,z2)%*%beta))
  D = -log(outcome[,1])/Lambda_D #survival time
  H = -log(outcome[,2])/Lambda_H #nonfatal event time
  #use the same covariates for both censoring and failure time
  x1 = z1   #using the same covariates as PIM model
  x2 = z2
  C = -log(runif(n))/exp(cbind(x1,x2)%*%betac) #censoring time
  status1 = 1*(D<=C)
  Do = pmin(C,D)
  status2 = rep(0,n)
  #if simulated HFH time is after death then 
  #we sill not observe the HFH. So we censor at 
  #the earliest of C years or time of death
  status2[H<=C & H<=Do] = 1;
  status2[H>Do] = 0
  H[H>Do] = Do[H>Do]
  Ho = H
  #Ho = pmin(C,H)
  data_set = NULL
  for(i in 1:n){
    if(status2[i]==1){
      row1 = c(i,Ho[i],2,z1[i],z2[i])
      row2 = c(i,Do[i],status1[i],z1[i],z2[i])
      data_set = rbind(data_set,row1,row2)
    }else{
      row2 = c(i,Do[i],status1[i],z1[i],z2[i])
      data_set = rbind(data_set,row2)
    }
  }
  colnames(data_set) = c("id","time","status","z1","z2")
  data_set = as.data.frame(data_set)
  return(data_set)
}

#Note that, under the same seed, the above two data generations are the same

#main simulation
N = 1000
n = 200  #sample size or 500
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
L = Inf  #no restriction event time
flag0 = flag1 = flag2 = matrix(0,N,length(case_set))
for(case in case_set){
  if(case == 1){         
    alpha = 1 
    beta = c(-0.5,0.5) 
  }else if(case == 2){  
    alpha = 1 
    beta = c(0,0)     
  }else if(case == 3){   
    alpha = 1 
    beta = c(0.5,-0.5) 
  }else if(case == 4){  
    alpha = 2 
    beta = c(-0.5,0.5)
  }else if(case == 5){  
    alpha = 2 
    beta = c(0,0)
  }else if(case == 6){   
    alpha = 2 
    beta = c(0.5,-0.5)
  }
  beta_est0 = se_est0 = cova0 = matrix(0, N, 2)
  beta_est1 = se_est1 = cova1 =matrix(0, N, 2)
  beta_est2 = se_est2 = cova2 =matrix(0, N, 2)
  for(i in 1:N){
    cat("iter=",i,"\n")
    set.seed(123456+i)
    #generate data
    gdata1 = gendata(n,beta,betac,lambda_h,lambda_d,alpha)
    set.seed(123456+i)
    gdata2 = gendata_full(n,beta,betac,lambda_h,lambda_d,alpha)
    #PWFR
    ID = gdata1$id
    time = gdata1$time
    status = gdata1$status
    Z = cbind(z1=gdata1$z1,z2=gdata1$z2)
    pwreg.obj = pwreg(ID=ID,time=time,status=status,Z=Z,rho=0)
    beta_est0[i,] = pwreg.obj$beta
    vare0 = pwreg.obj$Var
    se_est0[i,] = sqrt(diag(vare0))
    low0 = beta_est0[i,]-qna*se_est0[i,]
    high0 = beta_est0[i,]+qna**se_est0[i,]
    if(low0[1] <= beta[1] & beta[1] <= high0[1]) cova0[i,1] = 1
    if(low0[2] <= beta[2] & beta[2] <= high0[2]) cova0[i,2] = 1
    #GWFM with logit link
    res1 = pim.csd.m2(formula=y1~z1+z2,data=gdata2,SurvTime="y1",Status="delta1",covariateC=c("x1","x2"),
                     L = L, adjusted="unadjusted", priority = 1:2, model = "difference",start = start,
                     link = "logit", method = "Newton")
    flag1[i,case] = res1$flag
    beta_est1[i,] = res1$coef
    vare1 = res1$vcov
    se_est1[i,] = sqrt(diag(vare1))
    low1 = beta_est1[i,]-qna*se_est1[i,]
    high1 = beta_est1[i,]+qna*se_est1[i,]
    if(low1[1] <= beta[1] & beta[1] <= high1[1]) cova1[i,1] = 1 
    if(low1[2] <= beta[2] & beta[2] <= high1[2]) cova1[i,2] = 1 
    #GWFM with probit link
    res2 = pim.csd.m2(formula=y1~z1+z2,data=gdata2,SurvTime="y1",Status="delta1",covariateC=c("x1","x2"),
                      L = L, adjusted="unadjusted",priority = 1:2, model = "difference",start = start,
                      link = "probit", method = "Newton")
    
    flag2[i,case] = res2$flag
    #transform back by using delta method  
    beta_a = res2$coef
    pbeta = pnorm(beta_a)
    dbeta = dnorm(beta_a) 
    beta_est2[i,] = log(pbeta/(1-pbeta))
    vare_a = res2$vcov
    #y=ln(Phi(x)/(1-Phi(x)))--->y'=phi(x)/(Phi(x)*(1-Phi(x)))
    se_est2[i,] = abs(dbeta/(pbeta*(1-pbeta)))*sqrt(diag(vare_a)) #delta method
    low2 = beta_est2[i,]-qna*se_est2[i,]
    high2 = beta_est2[i,]+qna*se_est2[i,]
    if(low2[1] <= beta[1] & beta[1] <= high2[1]) cova2[i,1] = 1 
    if(low2[2] <= beta[2] & beta[2] <= high2[2]) cova2[i,2] = 1 
  }
  #summary estimation results
  index = flag1[,case]*flag2[,case] 
  est0 = colMeans(beta_est0[index==1,])
  est1 = colMeans(beta_est1[index==1,])
  est2 = colMeans(beta_est2[index==1,])
  
  #MCSD
  sd0 = apply(beta_est0[index==1,],2,sd)
  sd1 = apply(beta_est1[index==1,],2,sd)
  sd2 = apply(beta_est2[index==1,],2,sd)
  
  #AESE
  se0 = colMeans(se_est0[index==1,])
  se1 = colMeans(se_est1[index==1,])
  se2 = colMeans(se_est2[index==1,])
  
  #cp
  cp0 = colMeans(cova0[index==1,])
  cp1 = colMeans(cova1[index==1,])
  cp2 = colMeans(cova2[index==1,])
  
  #combine all results
  res0 = data.frame(est=est0,sd=sd0,se=se0,cp=cp0)
  res1 = data.frame(est=est1,sd=sd1,se=se1,cp=cp1)
  res2 = data.frame(est=est2,sd=sd2,se=se2,cp=cp2)
  
  resu = rbind(res0,res1,res2)
  res = round(resu,3)
  write.csv(res,paste0("PWFM_GWFM_case_",case,"_n_",n,".csv"))
  cat("case=",case,"\n")
}

