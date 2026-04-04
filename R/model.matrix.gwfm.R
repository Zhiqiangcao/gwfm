#' Create a model matrix for a generalized win fraction regression model
#' 
#' This function creates a model matrix for use in a generalized win fraction regression model (gwfm). 
#' This model matrix can be passed to \code{\link{gwfm_fit}}.
#' 
#' @param object a gwfm.formula object that contains the formula necessary for constructing the model matrix. 
#' @param data an optional argument specifying the data frame for which the model matrix should be constructed. 
#' @param model a single character value with possible values "difference" 
#' (the default), "marginal" or "customized". See also \code{\link{gwfm}}.
#' @param ... extra arguments passed to or from other methods.
#' This is currently only implemented in concordance with the generic \code{\link[stats]{model.matrix}} function.
#' 
#' @return a design matrix for a gwfm model
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}
#' @importFrom stats update
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats formula
#' @importFrom pim has.specials
#' @importFrom pim has.intercept
#' @importFrom pim poset
#' 

model.matrix.gwfm <- function(object, data, 
           model = c("difference","marginal","customized"), ...){
  model = match.arg(model)
  if(missing(data)) data = object@penv
  tt = object@terms
  specials = has.specials(object)
  
  if(specials){
    if(model != "customized") model = "customized"
      #warning("Argument model is changed to 'customized'")
    tt[[2]] = object@lhs
    tt = terms(formula(tt, env=data))
  }else if(model == "marginal"){
    tt = terms(update(tt, ~ . - 1))
  }
  mm = model.matrix(tt,data, ...)
  if(!specials){
    pos = poset(data, as.list=TRUE)
    if(model == "difference"){
      mm = mm[pos$L,,drop = FALSE] - mm[pos$R,,drop=FALSE]
    } else if(model == "marginal"){
      mm = mm[pos$L,,drop = FALSE]
    }
  }
  if(has.intercept(object)){
    if((id = match("(Intercept)",colnames(mm),0L)) > 0L){
      mm[,id] = 1
    } else {
      mm = cbind(mm, "(Intercept)" = 1)
    }
  }else{
    if((id = match("(Intercept)", colnames(mm),0L)) > 0L){
      mm = mm[,-id, drop=FALSE]
    }
  }
  return(mm)
}
