logit_df <- function(logit, X){
  ### logit_df() FUNCTION: calculating the degrees of freedom in a potentially penalized regression 
  ## INPUTS
  # logit: a logit model (needs to be specified further)
  # X: the data matrix WITHOUT INTERCEPT
  
  model <- logit
  X <- cbind(1, X) # Adding intercept
  names(X) <- c("X0", names(X)) # Giving intercept column X0 name
  # QUID: CAN WE GET X OUT OF THE MODEL SOMEHOW?
  
  ### Extracting some parameters
  ## enetLTS models
  if(class(model) == "enetLTS"){
    # Predictions
    p_hat <- model$fitted.values # These are the reweigted fitted values!
    
    # Coefficients
    coefs <- coefficients(model) # Gets coefficients (INCLUDING INTERCEPT)
    coefs_nonzero <- coefs[coefs != 0]
    index_coefs_nonzero <- as.numeric(names(coefs_nonzero)) # COLUMN NUMBERS of nonzero coefs
    # Note the column numbers start from 1! and e.g. go to 51
    
    ## Getting subset out
    wt <- model$raw.wt == 1
    
    X <- X[wt, ]
    p_hat <- p_hat[wt]
    
    
    
    # Alpha, Lambda(reweighted)
    lambda <- model$lambdaw
    alpha <- model$alpha
  } else {  ## Other kind of models
    p_hat <- predict(object = model, type = "response")
    lambda <- model$lambda # UNDER THE ASSUMPTION ONLY A SINGLE LAMBDA IS GIVEN
    
    
    
    
    
  }
  
  # Extracting the active set
  X <- X[, index_coefs_nonzero]
  
  # Calculating W matrix
  W <- diag(p_hat * (1 - p_hat))
  
  # Amount of obs
  nobs <- length(p_hat)
  
  # Calculating Hat (H) Matrix
  H <- sqrt(W) %*% X %*% solve((t(X) %*% W %*% X)  + lambda * (1 - alpha)/2) %*% t(X) %*% sqrt(W)
  
  # Calculating df 
  H_trace <- sum(diag(H))
  
  # OUTPUT
  return(H_trace)
  
}
