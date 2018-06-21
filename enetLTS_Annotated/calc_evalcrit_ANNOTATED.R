### ANNOTATING calc_evalCrit function
# This was split off from the cv.enetLTS() function

calc_evalCrit <- function(rowind, combis_ind, alphas, lambdas, index, xx, yy, nfold, repl){
  
  ## INPUT:
  # combis_ind: the result from expand.grid
  # rowind: a apply-type loop will go over all the row indices
  
  # Extracting indices out of the expand.grid
  i <- combis_ind[rowind, 1] # Column 1: lambda
  j <- combis_ind[rowind, 2] # Column 2: alpha
  
  # Using the indices to get the actual VALUE of lambda and alpha
  lambda <- lambdas[i] # Looking up in lambdas vector, position from earlier
  alpha <- alphas[j] # Same, but for alpha
  
  # Printing message which values are being used now
  print(paste("cross-validating for alpha: ", alpha, " and lambda :", 
              lambda), sep = "")
  
  # index: comes from the warmCsteps functions (see lvl 1)
  if (is.null(index)) {
    x <- xx
    y <- yy
  }
  else {
    x <- xx[index[, i, j], ] # WEIRD SUBSETTING MUST BE SOMEKIND OF 3 level strucure! [x, y, z]?
    y <- yy[index[, i, j]]
  }
  
  
  ### CROSS VALIDATION STRUCTURE:
  ## REPLICATION LEVEL (REPEATED CV)
  evalCritl <- rep(NA, repl)
  for (l in 1:repl) {
    
    # For logistic
    if (family == "binomial") {
      # Folds for "0" outcomes
      folds0 <- cvFolds(length(y[y == 0]), # cvTools::cvFolds
                        K = nfold,  
                        R = 1, 
                        type = "random")
      # Folds for "1" outcomes
      folds1 <- cvFolds(length(y[y == 1]), 
                        K = nfold, 
                        R = 1, 
                        type = "random")
      # Initiating losses for 0/1 seperately
      loss0 <- rep(NA, sum(y == 0))
      loss1 <- rep(NA, sum(y == 1))
    }
    
    # For Gaussian (Linear models)
    else if (family == "gaussian") {
      folds <- cvFolds(length(y), K = nfold, R = 1, 
                       type = "random")
      loss <- rep(NA, nrow(x))
    }
    
    ## FOLDS LEVEL
    for (f in 1:nfold) {
      
      # Logistic
      if (family == "binomial") {
        xtrain0 <- x[y == 0, ][folds0$subsets[folds0$which != f, 1], ]
        ytrain0 <- y[y == 0][folds0$subsets[folds0$which != f, 1]]
        xtest0 <- x[y == 0, ][folds0$subsets[folds0$which == f, 1], ]
        ytest0 <- y[y == 0][folds0$subsets[folds0$which == f, 1]]
        
        xtrain1 <- x[y == 1, ][folds1$subsets[folds1$which != f, 1], ]
        ytrain1 <- y[y == 1][folds1$subsets[folds1$which != f, 1]]
        xtest1 <- x[y == 1, ][folds1$subsets[folds1$which == f, 1], ]
        ytest1 <- y[y == 1][folds1$subsets[folds1$which == f, 1]]
        
        xtrain <- rbind(xtrain0, xtrain1)
        ytrain <- c(ytrain0, ytrain1)
        xtest <- rbind(xtest0, xtest1)
        ytest <- c(ytest0, ytest1)
      }
      
      # Gaussian (Linear)
      else if (family == "gaussian") {
        xtrain <- x[folds$subsets[folds$which != f, 1], ]
        ytrain <- y[folds$subsets[folds$which != f, 1]]
        xtest <- x[folds$subsets[folds$which == f, 1], ]
        ytest <- y[folds$subsets[folds$which == f, 1]]
      }
      
      # 
      res <- tryCatch({
        hpen <- length(ytrain)
        trainmod <- glmnet(xtrain, ytrain, family, 
                           alpha = alpha, lambda = lambda/hpen, standardize = FALSE, 
                           intercept = FALSE)
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
          loss0[folds0$which == f] <- -(ytest0 * xtest0 %*% matrix(trainmod$beta)) + log(1 + exp(xtest0 %*% matrix(trainmod$beta)))
          loss1[folds1$which == f] <- -(ytest1 * xtest1 %*% matrix(trainmod$beta)) + log(1 + exp(xtest1 %*% matrix(trainmod$beta)))
        }
        
        # Gaussian Loss:
        else if (family == "gaussian") 
          # Loss is also negative loglik = squared loss y - y_hat (y_hat = X*beta_hat)
          loss[folds$which == f] <- ytest - xtest %*% matrix(trainmod$beta)
      }
    }
    ## Going from loss to evaluation criterion
    # Binomial:
    if (family == "binomial") {
      loss <- c(loss0, loss1)
      # Taking Mean loglikelihood
      evalCritl[l] <- mean(loss, na.rm = TRUE)
    }
    # Gaussian:
    else if (family == "gaussian") {
      evalCritl[l] <- sqrt(mean(loss^2))
    }
  }
  return(list(lambda_ind = i, alpha_ind = j, evalCritl = evalCritl))
}