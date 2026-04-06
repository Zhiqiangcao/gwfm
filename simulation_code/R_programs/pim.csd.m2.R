#Probability index model for composite survival data

#parameter description:
#formula: An object of class formula (or one that can be coerced to that class): 
#         A symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#data: data set (matrix format) including id, multiple endpoints (e.g., y1, y2)
#      corresponding censoring status (e.g., delta1, delta2) and covariates
#SurvTime: survival time of fatal event (e.g., death)
#Status: status of fatal event
#covariateC: covariates used in Cox regression model for IPCW adjusted method
#            if adjusted = "unadjusted", then this setting is useless (or default NULL)
#priority: composite endpoints ranked in descending order of importance.
#          (e.g., priority = 1:2 means y1 is more important than y2)
#L: prespecified time horizon to construct truncated event endpoints (default Inf)
#adjusted: unadjusted or using IPCW adjusted method
#start: start value for fitting pim model for composite survival data (default NULL)
#control: control settings for nleqslv program when fitting pim model for composite survival data (default list())
#The other parameter settings are the same as pim program in R package pim 
#using type 2 weight
pim.csd.m2 = function(formula, data, SurvTime, Status, covariateC = NULL, priority = 1:2, 
                       L = Inf, t_val = 0.01, adjusted = c("unadjusted","IPCW"), start = NULL, 
                       control = list(), link = c("logit", "probit", "identity"),
                       model = c("difference", "marginal", "customized"), 
                       compare = if (model == "marginal") "all" else "unique", 
                       na.action = getOption("na.action"), keep.data = FALSE, ...){
  
  model = match.arg(model)
  if (is.character(compare)) {
    if (!compare %in% c("unique", "all")) 
      stop("compare should have the value 'unique' or 'all' in case it's a character value.")
  }
  
  # sort data by observed survival time
  index = order(data[,SurvTime])
  data.sort = data[index,]
  data = data.sort
  data = as.data.frame(data)
  link = match.arg(link)
  if (is.null(na.action)) na.action = "na.fail"
  if (!is.character(na.action)) na.action = deparse(substitute(na.action))
  
  f.terms = terms(formula, simplify = TRUE)
  vars = all.vars(formula)
  # prepare covariate matrix
  penv = new.pim.env(data, compare = compare, vars = vars, env = parent.frame())
  ff = new.pim.formula(formula, penv)
  x = model.matrix.pim.csd(ff, na.action = na.action, model = model)
    
  # next, we prepare comparison results between composite survival endpoints 
  # length of endpoints
  n_ep = length(priority)
  # since the data set format is id, y1, ..., y_nep, delta1, ..., delta_nep, z1,...,zp
  datap = as.matrix(data[,c(1:(1+2*n_ep))])
  # construct truncated event times and data if L != Inf
  y_o = as.matrix(data[,2:(n_ep+1)]) #original observed endpoints
  delta_o = as.matrix(data[,(n_ep+2):(1+2*n_ep)]) #original observed censoring indicators
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
  ties_flag = apply(determined_flag,1,sum)
  #y is the win result of W(Y_i,Y_j)
  y = apply(as.matrix(win_status[,seq(3,(2+2*n_ep-1),2)],ncol = n_ep),1,sum)
  #make those ties comparison results do not contribute in score function
  x = as.matrix(ties_flag*x)
  p = dim(x)[2]
  #when model = "difference" or "marginal", penv is also different
  penv = as.environment(penv@poset)
  if(is.null(start)){  #using zero as default start points
    start = rep(0, p)
  }
  if(adjusted == "unadjusted"){ 
    #no adjusted for win statistics or don't need to adjust
    weight_est = weight_var = NULL
  }else{ #IPCW adjusted
    tdata = data.frame(datap_n,data[,covariateC]) #truncated
    colnames(tdata) = c(colnames(data[,c(1:(1+2*n_ep))]),covariateC)
    censoring_weight = censoring.weight.est(Data=tdata,SurvTime=SurvTime,Status=Status,covariateC=covariateC)
    weight_res = ipcw.weight.m2(tdata,trt_con,censoring_weight,win_status,priority,determined_flag,L,n_ep,t_val=t_val)
    weight_est = weight_res$weight_est
    weight_var = weight_res$weight_var
  }
  #fit score and sandwich variance estimation
  est_res = pim.csd.fit(x=x, y=y, link = link, start = start, control = control, 
                       weight_est = weight_est, weight_var = weight_var, penv = penv, ...)
  #summary estimation results
  names(est_res$coef) = colnames(x)
  if (!keep.data) {
    x = matrix(nrow = 0, ncol = 0)
    y = numeric(0)
  }
  flag = est_res$flag
  Estimate = est_res$coef
  Std.Error = sqrt(diag(est_res$vcov))
  if(sum(is.na(Std.Error)>0) | sum(is.infinite(Std.Error))>0){
    flag = 0;
    Std.Error = rep(1,p)
  }
  z.value = Estimate/Std.Error
  p.value = rep(0,p)
  for(i in 1:p){
    z.valuei = z.value[i]
    p.value[i] = 2*pnorm(-abs(z.valuei))
  }
  coefficients = data.frame(Estimate,Std.Error,z.value,p.value)
  final_res = list(summary = coefficients, coef = est_res$coef, vcov = est_res$vcov, fitted = est_res$fitted, 
                   flag = flag, penv = penv, link = link, estimators = est_res$estim, 
                   model.matrix = x, response = y, keep.data = keep.data, model = model)
  return(final_res)
}
