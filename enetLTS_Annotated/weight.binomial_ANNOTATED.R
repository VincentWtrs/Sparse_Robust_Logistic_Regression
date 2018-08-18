weight.binomial_ANNOTATED <- function (x, y, beta, intercept, del) 
{
  ### weight.binomial(): Function that calculates the binary weight function for the reweighting step for binomial regression
  
  ## INPUTS
  # x: Robustly standardized (MAD scale + med center) predictors, intercept NOT included. Note: the whole dataset (n obs).
  # y: Just the outcomes 0/1, nothing done with it. Note: the whole dataset (n obs)
  # beta: Fitted coefficients from the regression, if intercept = TRUE, should include beta_hat_0 for the intercept.
  # intercept: TRUE/FALSE stating if an intercept was fitted or not.
  # del: The 1-quantile to flag something as an outlier. e.g. 0.0125 = 0.975 quantile and above are given 0-weights = outlier.
  
  ## OUPUT
  # A binary vector (0/1) of length n with weights for all observations. 1: OK, 0: NOK, outlier!
  
  ## Workings
  # 1. Fitted values based on the raw inputted coefs are calculated, if necessary intercept is removed or not because...
  #    it is always called with an added intercept from enetLTS() anyway. 
  # 2. Then Pearson-type residuals for all n are calculated as (y - pi_hat)/(sqrt(pi_hat * (1 - p_hat)))
  # 3. These are approx. Normal, hence their absolute values are compared with Normal percentile using a cutoff...
  #    values that are above the cutoff are given a 0-weight. The result is a n-vector of 0/1 weights
  
  # Case: with intercept
  if (intercept == TRUE) {
    # Getting predictions (prob hat)
    pi <- exp(cbind(1, x) %*% beta)/(1 + exp(cbind(1, x) %*%  beta))
    
    # Getting Pearson Deviances
    res <- (y - pi)/sqrt(pi * (1 - pi))
    # NOTE: in KHF (2017) paper wrong: they forgot the square root
    # NOTE: Agresti (2015, p.137) state that it would be better to divide res by: sqrt(1 - h_ii) (H: hat)
    
  # Case: without intercept
  } else {
    pi <- exp(x %*% beta[-1])/(1 + exp(x %*% beta[-1])) # Yes the intercept is inputted but just removed again
    res <- (y - pi)/sqrt(pi * (1 - pi))
  }
  
  # Binary weight function
  we <- as.integer(abs(res) <= qnorm(1 - del)) # If abs(res) < qnorm(1-del) -> 1
  # If above qnorm(1-del) e.g. if abs(res) > qnorm(0.975) -> assign 0 -> OUTLIER
  
  # OUTPUT
  return(we) # A binary vector of length n (not h!)
}
