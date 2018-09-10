warmCsteps <- function(x, y, h, n, p, family, alphas, lambdas, hsize, nsamp, s1, csteps, nfold, para, ncores, tol, scal, seed) 
{
  
  ### warmCsteps() FUNCTION
  ## Called by enetLTS with following inputs: x = x, y = y, h = h, n = n, p = p, family = family, alphas = alphas, 
  # lambdas = lambdas, hsize = hsize, nsamp = nsamp, s1 = s1, csteps = ncsteps, nfold = nfold, para = para, ncores = ncores,
  # tol = tol, scal = scal, seed = seed
  
  # Allocating STARTING alpha,lambda combination
  alpha <- alphas[1] # Allocating first alphas value (will be the SMALLEST)
  lambda <- lambdas[1] # Allocating first lambdas value (will be the LARGEST!)
  
  # Allocating residuals and index ARRAYS!
  residall <- array(NA, dim = c(n, length(lambdas), length(alphas))) # Array of n x #lambdas x #alphas
  indexall <- array(NA, dim = c(h, length(lambdas), length(alphas))) # Same but now for h: h x #lambdas x #alphas
  
  #### Beginning C-step
  ### Case: singla alpha-lambda combination (must be user specified in that case)
  if (length(alphas) == 1 & length(lambdas) == 1) {
    beginning.Cstep.with500 <- enetLTS:::beginningCstep(x = x, 
                                                        y = y, 
                                                        family = family, 
                                                        h = h, 
                                                        hsize = hsize, 
                                                        alpha = alpha, 
                                                        lambda = lambda, 
                                                        nsamp = nsamp, 
                                                        s1 = s1, 
                                                        ncores = ncores, 
                                                        csteps = csteps, 
                                                        tol = tol, 
                                                        scal = scal, 
                                                        para = para, 
                                                        seed = seed)
    residall[, 1, 1] <- beginning.Cstep.with500$resid
    indexall[, 1, 1] <- beginning.Cstep.with500$index
    return(list(residall = residall, indexall = indexall))
    
    ### Case: multiple alpha-lamda combinations
  } else { 
    beginning.Cstep.with500 <- enetLTS:::beginningCstep(x = x, 
                                                        y = y, 
                                                        family = family, 
                                                        h = h, 
                                                        hsize = hsize, 
                                                        alpha = alpha, 
                                                        lambda = lambda, 
                                                        nsamp = nsamp, 
                                                        s1 = s1, 
                                                        ncores = ncores, 
                                                        csteps = csteps, 
                                                        tol = tol, 
                                                        scal = scal, 
                                                        para = para, 
                                                        seed = seed)
    # Gathering results (initiation)
    index1_al <- beginning.Cstep.with500$index # h indices for alphas (# weird)
    index1_la <- beginning.Cstep.with500$index # h indices for lambdas (# weird)
    resid1_al <- beginning.Cstep.with500$resid # n residuals for alphas (# weird)
    resid1_la <- beginning.Cstep.with500$resid # n residuals for lambdas (# weird)
    
    ### ITERATIVELY (WARM) STARTING h-SUBSET SEARCH FOR DIFFERENT LAMBDA/ALPHA COMBINATIONS (This structure is quite complex)
    ## For each alpha
    for (al in 1:length(alphas)) {
      alpha <- alphas[al] # Get current alpha value (al = some index: 1:length(alphas))
      index1_la <- index1_al # Some redefining because for each new alpha, one "starts over" in some way
      resid1_la <- resid1_al # Some redefining because for each new alpha, one "starts over" in some way
      
      # Case: single lambda (maybe skip)
      if (length(lambdas) == 1) {
        newindex_la <- index1_la 
        objbest <- tol
        cstep.mod <- CStep(x = x, 
                           y = y, 
                           family = family, 
                           indx = newindex_la,
                           h = h,
                           hsize = hsize,
                           alpha = alpha,
                           lambda = lambda/h, 
                           scal = scal)
        
        countloop <- 0
        # while loop for running C-steps
        while ((cstep.mod$object > objbest) & (countloop < csteps)) {
          countloop <- countloop + 1 # Keeping track of loops
          objbest <- cstep.mod$object # Objective function value
          newindex_la <- cstep.mod$index # h indices
          newresid_la <- cstep.mod$residu # n residuals
          cstep.mod <- CStep(x = x, 
                             y = y, 
                             family = family, 
                             indx = newindex_la, 
                             h = h, 
                             hsize = hsize, 
                             aplha = alpha, 
                             lambda = lambda/h, 
                             scal = scal)
          
          index1_la <- newindex_la
        } # End of while-loop
        
        indexall[, , al] <- newindex_la # Save in indexall for the current alpha
        residall[, , al] <- newresid_la # Save in the residuals for the current alpha
      
        
      ## CASE: MULTIPLE LAMBDAS! (Interesting case)
      } else {
        # Constructing matrices, rows = h, columns is lambdas 
        IndexMatrix <- matrix(NA, nrow = h, ncol = (length(lambdas) -  1)) # lambda - 1 because the first lambda value (default: lambda0) was used for the warm start already
        ResidMatrix <- matrix(NA, nrow = n, ncol = (length(lambdas) -  1)) # lambda - 1 because the first lambda value (default: lambda0) was used for the warm start already
        
        # For each lambda
        for (la in 1:(length(lambdas) - 1)) { # (Alternatively, and clearer: for(la in 2:length(lambda))) 
          lambda <- lambdas[la + 1] # Extracting lambda to use for this loop
          newindex_la <- index1_la # Using the last arrived-upon h observation indices
          objbest <- tol
          
          # Running a single C-step (to have at least 1 run before the while, actually to just initiate the whole thing if the while-loop doesn't even need to be ran)
          cstep.mod <- CStep(x = x, 
                             y = y, 
                             family = family, 
                             indx = newindex_la, 
                             h = h , 
                             hsize = hsize, 
                             alpha = alpha, 
                             lambda = lambda/h, 
                             scal = scal)
          # Note: the fact that a single C-step is done here, probably just to be sure that at least one is done (?)
          
          # Running (potentially) many C-steps
          countloop <- 0 # Initiating the loop counter
          while ((cstep.mod$object > objbest) & (countloop <  csteps)) {
            countloop <- countloop + 1
            objbest <- cstep.mod$object # Saving objective function value from the LAST run 
            newindex_la <- cstep.mod$index # Saving h observation indices for this lambda from the LAST run (Because they (POTENTIALLY) changed as a result from the Cstep!)
            newresid_la <- cstep.mod$residu # Saving n residuals from this lambda from the LAST run ...
            # Note: for the first run these are the results from the single run above ...
            
            # Running many C-steps potentially
            cstep.mod <- CStep(x = x, 
                               y = y, 
                               family = family, 
                               indx = newindex_la, 
                               h = h, 
                               hsize = hsize, 
                               alpha = alpha, 
                               lambda = lambda/h, 
                               scal = scal)
            index1_la <- newindex_la # This makes sure that the last C-step result will be used for the following lambda!
          } # End of while-loop
          
          # Saving the results after the many C-Steps
          IndexMatrix[, la] <- newindex_la # The result is also assigned  (after the while-loop) as the final h subset for that lambda
          ResidMatrix[, la] <- newresid_la # Same for the n residuals
        }
        # For the first lambda (Don't fully understand why another series of C-steps is done, but in any case, doesn't do any harm!)
        lambda <- lambdas[1] # First lambdas (often lambda_0)
        newindex_al <- index1_al # Why with the index from the last one, this is kind of weird... (?)
        objbest <- tol
        cstep.mod <- CStep(x = x, 
                           y = y, 
                           family = family, 
                           indx = newindex_al, 
                           h = h, 
                           hsize = hsize, 
                           alpha = alpha, 
                           lambda = lambda/h, 
                           scal = scal)
        countloop <- 0
        while ((cstep.mod$object > objbest) & (countloop < csteps)) {
          countloop <- countloop + 1
          objbest <- cstep.mod$object
          newindex_al <- cstep.mod$index
          newresid_al <- cstep.mod$residu
          cstep.mod <- CStep(x = x, 
                             y = y, 
                             family = family, 
                             indx = newindex_al, 
                             h = h, 
                             hsize = hsize, 
                             alpha = alpha, 
                             lambda = lambda/h, 
                             scal = scal)
          index1_al <- newindex_al
        }
        # Combining the results from the first lambda with the other lambdas
        IndexMatrix <- cbind(newindex_al, IndexMatrix)
        ResidMatrix <- cbind(newresid_al, ResidMatrix)
        indexall[, , al] <- IndexMatrix
        residall[, , al] <- ResidMatrix
      }
    }
  }
  # OUTPUT
  return(list(indexall = indexall, residall = residall))
}