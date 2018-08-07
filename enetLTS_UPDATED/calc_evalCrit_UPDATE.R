# CACL_EVALCRIT_UPDATE: Allowing nfold = 1 FOR INFORMATION CRITERION-BASED APPROACHES
calc_evalCrit_UPDATE <- function(rowind, combis_ind, alphas, lambdas, 
                          index, xx, yy, nfold, repl, family, ic_type = NULL) {
  # family argument defined as well, because it was defined within cv.enetLTS and calc_evalCrit just used it there as well
  
  if(!is.null(ic_type)){
    ic <- TRUE
  } else if(is.null(ic_type)){
    ic <- FALSE
  }
  
  i <- combis_ind[rowind, 1]
  j <- combis_ind[rowind, 2]
  lambda <- lambdas[i]
  alpha <- alphas[j]
  print(paste("cross-validating for WEJOOOOWABUDABI alpha: ", alpha, " and lambda :", 
              lambda), sep = "")
  if (is.null(index)) {
    x <- xx
    y <- yy
  }
  else {
    x <- xx[index[, i, j], ]
    y <- yy[index[, i, j]]
  }
  evalCritl <- rep(NA, repl)
  for (l in 1:repl) {
    
    ######### NEW: if(nfold > 1): Keep old functionality
    if(nfold > 1){ # BEGIN if(nfold > 1)
      
      
      if (family == "binomial") {
        folds0 <- cvFolds(length(y[y == 0]), K = nfold, 
                          R = 1, type = "random")
        folds1 <- cvFolds(length(y[y == 1]), K = nfold, 
                          R = 1, type = "random")
        loss0 <- rep(NA, sum(y == 0))
        loss1 <- rep(NA, sum(y == 1))
      }
      else if (family == "gaussian") {
        folds <- cvFolds(length(y), K = nfold, R = 1, 
                         type = "random")
        loss <- rep(NA, nrow(x))
      }
      for (f in 1:nfold) {
        if (family == "binomial") {
          xtrain0 <- x[y == 0, ][folds0$subsets[folds0$which != f, 1], ]
          ytrain0 <- y[y == 0][folds0$subsets[folds0$which != f, 1]]
          xtest0 <- x[y == 0, ][folds0$subsets[folds0$which == f, 1], ]
          ytest0 <- y[y == 0][folds0$subsets[folds0$which == f, 1]]
          
          xtrain1 <- x[y == 1, ][folds1$subsets[folds1$which != f, 1], ]
          ytrain1 <- y[y == 1][folds1$subsets[folds1$which !=  f, 1]]
          xtest1 <- x[y == 1, ][folds1$subsets[folds1$which == f, 1], ]
          ytest1 <- y[y == 1][folds1$subsets[folds1$which == f, 1]]
          
          xtrain <- rbind(xtrain0, xtrain1)
          ytrain <- c(ytrain0, ytrain1)
          xtest <- rbind(xtest0, xtest1)
          ytest <- c(ytest0, ytest1)
        }
        else if (family == "gaussian") {
          xtrain <- x[folds$subsets[folds$which != f, 1], ]
          ytrain <- y[folds$subsets[folds$which != f, 1]]
          xtest <- x[folds$subsets[folds$which == f, 1], ]
          ytest <- y[folds$subsets[folds$which == f, 1]]
        }
        res <- tryCatch({
          hpen <- length(ytrain)
          trainmod <- glmnet(xtrain, 
                             ytrain, 
                             family, 
                             lambda = lambda/hpen,
                             alpha = alpha, 
                             standardize = FALSE, 
                             intercept = FALSE)
        }, error = function(err) {
          error <- TRUE
          return(error)
        })
        if (is.logical(res)) {
          print(paste("CV broke off for alpha=", alpha, 
                      "and lambda=", lambda))
        }
        else {
          trainmod <- res
          if (family == "binomial") {
            loss0[folds0$which == f] <- -(ytest0 * xtest0 %*% matrix(trainmod$beta)) + log(1 + exp(xtest0 %*% matrix(trainmod$beta))) # ORIGINAL
            loss1[folds1$which == f] <- -(ytest1 * xtest1 %*% matrix(trainmod$beta)) + log(1 + exp(xtest1 %*% matrix(trainmod$beta))) # ORIGINAL
          }
          else if (family == "gaussian") 
            loss[folds$which == f] <- ytest - xtest %*% 
              matrix(trainmod$beta)
        }
      }
      if (family == "binomial") {
        loss <- c(loss0, loss1)
        evalCritl[l] <- mean(loss, na.rm = TRUE)
      }
      else if (family == "gaussian") {
        evalCritl[l] <- sqrt(mean(loss^2))
      }
    } # END OF if(nfold > 1)
    
    ########## NEW: if(nfold = 1):
    if(nfold == 1){ # BEGINNING if(nfold == 1)
      loss <- rep(NA, nrow(x))
      xtrain <- x
      ytrain <- y
      
      # Fitting (within an error catching structure)
      res <- tryCatch({
        hpen <- length(ytrain)
        trainmod <- glmnet(xtrain, 
                           ytrain, 
                           family, 
                           alpha = alpha, 
                           lambda = lambda/hpen, # NEEDS TO BE lambda/hpen!
                           standardize = FALSE, 
                           intercept = FALSE)
      }, error = function(err) {
        error <- TRUE
        return(error)
      })
      if (is.logical(res)) {
        print(paste("Fitting broke off for alpha=", alpha, 
                    "and lambda=", lambda))
      }
      else {
        trainmod <- res
        # Binomial Case
        if (family == "binomial") {
          # NEW
          loss <- -(ytrain * xtrain %*% matrix(trainmod$beta)) + log(1 + exp(xtrain %*% matrix(trainmod$beta)))
        }
        
        # Gaussian Case
        else if (family == "gaussian") 
          loss <- ytrain - xtrain %*% matrix(trainmod$beta)
      }
      
      ## Calculating evaluation criteria
      # Binomial
      if (family == "binomial") {
        evalCritl[l] <- mean(loss, na.rm = TRUE)
        if(ic == TRUE){
          evalCritl[l] <- 2 * mean(loss, na.rm = TRUE) + ic_penalty(trainmod, type = ic_type)
            #evalCritl[l] + ic_penalty(trainmod, type = ic_type)
        }
      } # Gaussian
      else if (family == "gaussian") {
        evalCritl[l] <- sqrt(mean(loss^2))
        if(ic == TRUE){
          evalCritl[l] <- 2 * evalCritl[l] + ic_penalty(trainmod, type = ic_type) # Minus loss + Penalty
        }
      }
      
    } # END OF if(nfold == 1)
    
    
  } # END OF REPL. LOOP
  return(list(lambda_ind = i, alpha_ind = j, evalCritl = evalCritl))
}