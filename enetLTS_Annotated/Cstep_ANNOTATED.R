CStep_ANNOTATED <- function (x, y, family, indx, h, hsize, alpha, lambda, scal) 
{
  n <- nrow(x)
  
  ## Data Preparation
  if (scal) {
    # Nonrobust Scaling (Data going in already Normalized robust!)
    # Why not robust: the data should be outlier free already (In the beginning it's even only 3/4 obs)
    scl <- enetLTS:::prepara(x, y, family, indx, robu = 0) 
    xs <- scl$xnor # Predictors normalized
    ys <- scl$ycen # Outcome centered (BUT NOT FOR BINOMIAL (would make no sense))
    
    #### Fitting the glmnet
    ### Binomial Case
    if (family == "binomial") {
      fit <- glmnet(x = xs[indx, ], # Data is subsetted using the supplied indx
                    y = ys[indx],  # Same for outcome
                    family = family, 
                    alpha = alpha, 
                    lambda = lambda, 
                    standardize = FALSE,  # Already done manually
                    intercept = FALSE) # No because we standardized the data ??? THIS DOES HOLD FOR THE ELEMENTAL SUBSETS AS THEY ARE 50/50 FOR 0/1 OUTCOMES BY CONSTRUCTION
      # THIS FITTING IS ACTUALLY THROWING THE ERROR WHEN DOING THE FITTING FOR THE ELEMENTAL SUBSETS (OF SIZE 3/4), THIS IS NOT PROBLEMATIC
      # Warning: In lognet(x, is.sparse, ix, jx, y, weights, offset, alpha, nobs, : one multinomial or binomial class has fewer than 8  observations; dangerous ground
      
      ## Calculating Residuals ON ALL (STANDARDIZED, BUT POSSIBLY NONROBUST) DATA (i = 1, ..., n)
      beta <- matrix(fit$beta) # Extracting beta (Here it's just a single value so the fit$beta object is just a matrix)
      resid <- -(ys * xs %*% beta) + log(1 + exp(xs %*% beta)) # Deviance of obs i = 1, ..., n (not h!)
      
      # If all betas are 0
      if (all(beta == 0)) { # In case beta == 0 (?)
        return(list(object = -Inf, 
                    index = indx, 
                    residu = resid, 
                    beta = beta))
      }
      
      # Sorting residuals
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE) # Decreasing! (Other for Gaussian)
      # By setting index.return = TRUE, the returned object is a LIST with $x: the sorted object and $ix: the indices of orders of the ORIGINAL object! (i.e. resid, not resid.sort)
      
      # Creating h0 and h1 (for 0/1 outcomes)
      h0 <- floor((length(y[y == 0]) + 1) * hsize) # h0: h for 0-outcomes
      h1 <- h - h0 # h1: h for 1-outcomes (h0 + h1 = h)
      
      # Getting indices for 0 and 1-outcomes respecitvely!
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0] 
      # y[resid.sort$ix] creates logical (T/F) vector of 0 outcomes (length of this logical = n!)
      # resid.sort$ix[y[resid.sort$ix] == 0]: selects the residual INDICES!! only for the 0 outcomes so length is n0
      # resid.sort$ix[y[resid.sort$ix] == 0][1:h0]: only picks the h0 smallest of these
      index1 <- resid.sort$ix[y[resid.sort$ix] == 1][1:h1] # Same for 1-outcomes
      indxnew <- c(index0, index1) # Combining, to get length h (75)
      
      ## GAUSSIAN CASE
    } else if (family == "gaussian") {
      fit <- glmnet(x = xs[indx, ], 
                    y = ys[indx], 
                    family = family, 
                    alpha = alpha, 
                    lambda = lambda, 
                    standardize = FALSE, 
                    intercept = FALSE) # For linear models this seems reasonable
      beta <- matrix(fit$beta)
      resid <- ys - predict(fit, xs, exact = TRUE) # Increasing!
      resid.sort <- sort(abs(resid), index.return = TRUE)
      indxnew <- resid.sort$ix[1:h]
    }
    
    ## Calculating the value of the objective function
    obj <- enetLTS:::Objval(xs, ys, family, beta, indxnew, alpha, lambda) #  This gives a single value!!!!
    
  # IF NO SCALING REQUIRED! (Don't really see why this should be so different)  (basically the whole thing goes from xs to x, redundant, just redefine in beginning???)
  } else { # NO SCALING (?)
    if (family == "binomial") {
      fit <- glmnet(x = x[indx, ], 
                    y = y[indx], 
                    family = family, 
                    alpha = alpha, 
                    lambda = lambda, 
                    standardize = FALSE, 
                    intercept = FALSE) # This seems inappropriate for binomial!
      beta <- matrix(fit$beta)
      resid <- -(y * x %*% beta) + log(1 + exp(x %*% beta))
      if (all(beta == 0)) {
        return(list(object = -Inf, index = indx, residu = resid, 
                    beta = beta))
      }
      resid.sort <- sort(resid, decreasing = FALSE, index.return = TRUE)
      h0 <- floor((length(y[y == 0]) + 1) * hsize)
      h1 <- h - h0
      index0 <- resid.sort$ix[y[resid.sort$ix] == 0][1:h0]
      index1 <- resid.sort$ix[y[resid.sort$ix] == 1][1:h1]
      indxnew <- c(index0, index1)
    } # Gaussian
    else if (family == "gaussian") {
      fit <- glmnet(x = x[indx, ], 
                    y = y[indx], 
                    family = family, 
                    alpha = alpha, 
                    lambda = lambda, 
                    standardize = FALSE, 
                    intercept = FALSE) # For linear models its appropriate after standardizing
      beta <- matrix(fit$beta)
      resid <- y - predict(fit, x, exact = TRUE)
      resid.sort <- sort(abs(resid), index.return = TRUE)
      indxnew <- resid.sort$ix[1:h]
    }
    obj <- enetLTS:::Objval(x, y, family, beta, indxnew, alpha, lambda)
  }
  # OUTPUT
  return(list(object = obj, 
              index = indxnew, 
              residu = resid, 
              beta = beta))
  # So the output is basically saying: e.g. in the beginnning: for the elemental subset you supplied, I fitted a glmnet, got the residuals and robust objective and extracted the best 75 observations
}