model.matrix.pim.csd <- function(object, data, 
           model = c("difference","marginal",
                     "regular","customized"), ...){
  model <- match.arg(model)
  if(missing(data)) data <- object@penv
  tt <- terms(object)
  specials <- has.specials(object)
  
  if(specials){
    if(model != "customized") model <- "customized"
      #warning("Argument model is changed to 'customized'")
    tt[[2]] <- object@lhs
    tt <- terms(formula(tt, env=data))
  }else if(model == "marginal"){
    tt <- terms(update(tt, ~ . - 1))
  }
  mm <- model.matrix(tt,data, ...)
  if(!specials){
    pos <- poset(data, as.list=TRUE)
    if(model == "difference"){
      mm <- mm[pos$L,,drop = FALSE] - mm[pos$R,,drop=FALSE]
    } else if(model == "marginal"){
      mm <- mm[pos$L,,drop = FALSE]
    }
  }
  if(has.intercept(object)){
    if((id <- match("(Intercept)",colnames(mm),0L)) > 0L){
      mm[,id] <- 1
    } else {
      mm <- cbind(mm, "(Intercept)" = 1)
    }
  }else{
    if((id <- match("(Intercept)",colnames(mm),0L)) > 0L){
      mm <- mm[,-id, drop=FALSE]
    }
  }
  mm
}