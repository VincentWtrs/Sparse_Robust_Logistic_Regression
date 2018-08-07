### ANNOTATING calc_evalCrit function
# This was split off from the cv.enetLTS() function

calc_evalCrit <- function(rowind, combis_ind, alphas, lambdas, index, xx, yy, nfold, repl){
  
  ## INPUT:
  # combis_ind: the result from expand.grid
  # rowind: a apply-type loop will go over all the row indices
  
  ## SO IF I'M RIGHT, THIS IS CALLED FOR A SINGLE ALPHA/LAMBDA VALUE
  
  # Extracting indices out of the expand.grid
  i <- combis_ind[rowind, 1] # Column 1: lambda INDEX
  j <- combis_ind[rowind, 2] # Column 2: alpha INDEX
  # So one would start at 1-1
  
  # Extracting the lambda-alpha VALUE combination based on the i,j from above
  lambda <- lambdas[i] # Looking up in lambdas vector, position from earlier
  alpha <- alphas[j] # Same, but for alpha
  
  # Printing message which values are being used now
  print(paste("cross-validating for alpha: ", alpha, " and lambda :", 
              lambda), sep = "")
  
  # index: comes from the warmCsteps functions (see lvl 1)
  if (is.null(index)) { # So this should not be the case
    x <- xx
    y <- yy
  }
  # Indeed, indexes (amount = h) should be received from warmCstep(), i.e. getting the indices in set H of size h (?)
  else {
    x <- xx[index[, i, j], ] # Gives one dataset of size h (X)
    y <- yy[index[, i, j]]  # Same (y)
  }
  
  
  #### CROSS VALIDATION STRUCTURE:
  ### REPLICATION LEVEL (REPEATED CV)
  evalCritl <- rep(NA, repl) # Creating NA vector with element for each replication
  for (l in 1:repl) {
    
    ## For Binomial - logistic
    if (family == "binomial") {
      # Folds for "0" outcomes (Need to have similar 0/1 balance)
      folds0 <- cvFolds(length(y[y == 0]), # cvTools::cvFolds
                        K = nfold,  
                        R = 1, # 1 replication, the replications are done manually in the loop that is defined in this function
                        type = "random") # Default: random folds
      # outcome of cvFolds(): $subset an integer vector (matrix with ncols = amount of reps, here = 1) containing a permutation of indices
      # $which = an integer vector stating to which fold each permuted observation belongs (so here generally filled with integers from 1 to 5 (K))
      # So $subset is basically a permutation (shuffled) of the indices
  
      # Folds for "1" outcomes
      folds1 <- cvFolds(length(y[y == 1]), 
                        K = nfold, 
                        R = 1, 
                        type = "random")
      # Initiating losses for 0/1 seperately
      loss0 <- rep(NA, sum(y == 0))
      loss1 <- rep(NA, sum(y == 1))
    }
    
    ## For Gaussian (Linear models)
    else if (family == "gaussian") {
      folds <- cvFolds(length(y), K = nfold, R = 1, 
                       type = "random")
      loss <- rep(NA, nrow(x))
    }
    
    ### FOLDS LEVEL (whats actually happening here? are the folds passed on?)
    for (f in 1:nfold) {
      
      ## For Binomial - Logistic
      if (family == "binomial") {
        # For the y = 0 outcomes
        xtrain0 <- x[y == 0, ][folds0$subsets[folds0$which != f, 1], ]
        ytrain0 <- y[y == 0][folds0$subsets[folds0$which != f, 1]]
        xtest0 <- x[y == 0, ][folds0$subsets[folds0$which == f, 1], ]
        ytest0 <- y[y == 0][folds0$subsets[folds0$which == f, 1]]
        # The [folds0$subsets[folds0$which == f, 1], ] is The internal one basically says: pick the observations that are equal to the fold we're currently at in the loop
        # the [, 1] is a bit redundant here I think, can drop it but guess it's a bit safer here.
        
        # For the y = 1 outcomes
        xtrain1 <- x[y == 1, ][folds1$subsets[folds1$which != f, 1], ]
        ytrain1 <- y[y == 1][folds1$subsets[folds1$which != f, 1]]
        xtest1 <- x[y == 1, ][folds1$subsets[folds1$which == f, 1], ]
        ytest1 <- y[y == 1][folds1$subsets[folds1$which == f, 1]]
        
        # Combining 0/1 into 1 data structure
        xtrain <- rbind(xtrain0, xtrain1)
        ytrain <- c(ytrain0, ytrain1)
        xtest <- rbind(xtest0, xtest1)
        ytest <- c(ytest0, ytest1)
      }
      
      ## Gaussian (Linear)
      else if (family == "gaussian") {
        xtrain <- x[folds$subsets[folds$which != f, 1], ]
        ytrain <- y[folds$subsets[folds$which != f, 1]]
        xtest <- x[folds$subsets[folds$which == f, 1], ]
        ytest <- y[folds$subsets[folds$which == f, 1]]
      }
      
      ## FITTING THE GLMNET MODEL
      res <- tryCatch({
        hpen <- length(ytrain) # Extracting size of H (h)
        
        # Fitting glmnet (not cv.glmnet!)
        trainmod <- glmnet(xtrain, 
                           ytrain, 
                           family, 
                           alpha = alpha, 
                           #lambda = lambda/hpen, # The lambda is scaled! (Probably different definition ?) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           # BETTER WOULD BE TO NOT SPECIFY ANY LAMBDA AND THEN GET THE COEFFICIENTS WITH THE LAMBDA WANTED
                           standardize = FALSE, 
                           intercept = FALSE)
        # LOOKS LIKE THEY DON'T USE THE LAMBDA WARM STARTS BECAUSE IT'S CALLED AGAIN EACH TIME
      }, error = function(err) { # Defining error function or so
        error <- TRUE
        return(error)
      })
      # When error is occuring the result is a LOGICAL value (?)
      if (is.logical(res)) {
        print(paste("CV broke off for alpha=", alpha, 
                    "and lambda=", lambda))
      }
      else {
        trainmod <- res
        
        # Binomial Loss:
        if (family == "binomial") {
          # Loss is negative loglikelihood in form of - y*eta + log(1 + exp(eta))
          # NOTE: trainmod$beta is a "dgCmatrix" object of size nobs (amount of obs) x nvars (amount of predictors)
          loss0[folds0$which == f] <- -(ytest0 * xtest0 %*% matrix(trainmod$beta)) + log(1 + exp(xtest0 %*% matrix(trainmod$beta))) # ORIGINAL
          loss1[folds1$which == f] <- -(ytest1 * xtest1 %*% matrix(trainmod$beta)) + log(1 + exp(xtest1 %*% matrix(trainmod$beta))) # ORIGINAL

          # NOTE2: because the glmnet model is fitted with a single alpha-lambda combination, the dgCmatrix simplifies to an ordinary one
        }
        
        # Gaussian Loss:
        else if (family == "gaussian") 
          # Loss is also negative loglik = squared loss y - y_hat (y_hat = X*beta_hat)
          #loss[folds$which == f] <- ytest - xtest %*% matrix(trainmod$beta) # ORIGINAL
          loss[folds$which == f] <- ytest - xtest %*% matrix(coef(traindmod, s = lamda/h)) # UPDATED
          
      }
    } # END OF FOLD LEVEL
    
    ## Going from loss to evaluation criterion (= avg. loss)
    # Binomial:
    if (family == "binomial") {
      loss <- c(loss0, loss1) # Combining loss for 0 and 1 outcomes
      # Taking Mean loglikelihood (evalCrit is just the MEAN LOSS)
      evalCritl[l] <- mean(loss, na.rm = TRUE) # This is the mean PER REPLICATION (I think only over h (????))
    }
    # Gaussian:
    else if (family == "gaussian") {
      evalCritl[l] <- sqrt(mean(loss^2))
    }
  }
  return(list(lambda_ind = i, alpha_ind = j, evalCritl = evalCritl))
}