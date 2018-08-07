ic_penalty <- function(glmnet_in, type){
  # This function aims to calculate information criteria PENALTY factor (which then needs to be ADDEd to loglik(= minus loss))
  ## INPUTS:
  # glmnet_in: a model fitted through glmnet()
  # type: the sort of information criterion wanted, can be: AIC, BIC at the moment
  
  ## Error Catching
  # Checking if model is glmnet
  #if(class(glmnet_in) != c("lognet", "glmnet")){
  #  stop("The input for ic_penalty() should be a glmnet model")
  #}
  # Checking if type is supported
  if(type != "AIC" & type != "BIC"){
    stop("Only information criteria supported: AIC, BIC")
  }
  
  # Extracting parameters
  model <- glmnet_in # Renaming for convenience
  #alpha <- model$call$alpha # Extracting alpha DOESN'T WORK
  alpha <- 1 # NEW: FORCING TO 1 JUST FOR THE MOMENT
  h <- model$nobs
  coefs_nonzero <- model$df
  
  
  # If alpha = 1 (LASSO) then the degrees of freedom (df) is equal to amount of nonzero coefs
  if(alpha == 1){
    df <- coefs_nonzero # For LASSO: simple
  }
  
  # Calculating the information criteria
  if(type == "AIC"){
    penalty <- 2 * df/h
  }
  
  if(type == "BIC"){
    penalty <- (df * log(n))/h
  }
  return(penalty)
}