ic_penalty <- function(glmnet_in, type, X, alpha, intercept){
  # This function aims to calculate information criteria PENALTY factor (which then needs to be ADDEd to loglik(= minus loss))
  
  ## INPUTS:
  # glmnet_in: a model fitted through glmnet()
  # type: the sort of information criterion wanted, can be: AIC, BIC at the moment
  # X: the data
  # alpha: because for some stupid reason one cannot extract it from the glmnet function.
  # intercept: if intercept was used to fit the logit-glmnet (family = "binomial") model: TRUE.
  
  ## OUPUTS:
  # Numeric value denoting the effective degrees of freedom
  
  ## Error Catching
  # Checking if model is glmnet
  #if(class(glmnet_in) != c("lognet", "glmnet")){
  #  stop("The input for ic_penalty() should be a glmnet model")
  #}
  
  # Checking if type is supported
  if(type != "AIC" & type != "BIC" & type != "EBIC"){
    stop("Only information criteria supported: AIC, BIC, EBIC")
  }
  
  # Extracting parameters
  model <- glmnet_in # Renaming for convenience
  h <- model$nobs
  coefs_nonzero <- model$df # Intercept or not?
  p <- length(coefficients(model))
  
  
  # If alpha = 1 (LASSO) then the degrees of freedom (df) is equal to amount of nonzero coefs
  if(alpha == 1){
    df <- coefs_nonzero # For LASSO: simple
  } else { # ONLY enetLTS() MODELS SUPPORTED!!!
    df <- logit_df(logit = glmnet_in, X)
    
  }
  
  ## Calculating the information criteria penalty
  # Akaike Information Criterion (AIC)
  if(type == "AIC"){
    penalty <- (2 * df)/h
  }
  
  # Bayesian Information Criterion (BIC)
  if(type == "BIC"){
    penalty <- (df * log(h))/h
  }
  
  # Extended Bayesian Information Criterion (EBIC)
  if(type == "EBIC"){
    sigma <- 0.25 # Some default value!
    penalty <- (coefs_nonzero * log(h) + 2 * coefs_nonzero * sigma * log(p))/h
  }
  
  # Generalized Information Criterion (GIC)
  if(type == "GIC"){
    penalty <- 0
  }
  
  return(penalty)
}
