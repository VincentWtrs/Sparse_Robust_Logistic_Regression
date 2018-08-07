warmCsteps <- function(x, y, h, n, p, family, alphas, lambdas, hsize, nsamp, 
                       s1, csteps, nfold, para, ncores, tol, scal, seed) 
{
  
  ### warmCsteps() FUNCTION
  ## Called by enetLTS with following inputs: x = x, y = y, h = h, n = n, p = p, family = family, alphas = alphas, 
  # lambdas = lambdas, hsize = hsize, nsamp = nsamp, s1 = s1, csteps = ncsteps, nfold = nfold, para = para, ncores = ncores,
  # tol = tol, scal = scal, seed = seed
  
  alpha <- alphas[1] # Allocating first alphas value (?)
  lambda <- lambdas[1] # Allocating first lambdas value (?)
  
  residall <- array(NA, dim = c(n, length(lambdas), length(alphas))) # Array of n x #lambdas x #alphas
  indexall <- array(NA, dim = c(h, length(lambdas), length(alphas))) # Same but now for h!
  
  ## Beginning C step
  # If only a single lambda and alpha value were given (must be user specified in that case)
  if (length(alphas) == 1 & length(lambdas) == 1) {
    beginning.Cstep.with500 <- beginningCstep(x, 
                                              y, 
                                              family, 
                                              h, 
                                              hsize, 
                                              alpha, 
                                              lambda, 
                                              nsamp, 
                                              s1, 
                                              ncores, 
                                              csteps, 
                                              tol, 
                                              scal, 
                                              para, 
                                              seed)
    residall[, 1, 1] <- beginning.Cstep.with500$resid
    indexall[, 1, 1] <- beginning.Cstep.with500$index
    return(list(residall = residall, indexall = indexall))
  }
  else {
    beginning.Cstep.with500 <- beginningCstep(x, 
                                              y, 
                                              family, 
                                              h, 
                                              hsize, 
                                              alpha, 
                                              lambda, 
                                              nsamp, 
                                              s1, 
                                              ncores, 
                                              csteps, 
                                              tol, 
                                              scal, 
                                              para, 
                                              seed)
    index1_al <- beginning.Cstep.with500$index
    index1_la <- beginning.Cstep.with500$index
    resid1_al <- beginning.Cstep.with500$resid
    resid1_la <- beginning.Cstep.with500$resid
    for (al in 1:length(alphas)) {
      alpha <- alphas[al]
      index1_la <- index1_al
      resid1_la <- resid1_al
      if (length(lambdas) == 1) {
        newindex_la <- index1_la
        objbest <- tol
        cstep.mod <- CStep(x, y, family, newindex_la, 
                           h, hsize, alpha, lambda/h, scal)
        countloop <- 0
        while ((cstep.mod$object > objbest) & (countloop < 
                                               csteps)) {
          countloop <- countloop + 1
          objbest <- cstep.mod$object
          newindex_la <- cstep.mod$index
          newresid_la <- cstep.mod$residu
          cstep.mod <- CStep(x, y, family, newindex_la, 
                             h, hsize, alpha, lambda/h, scal)
          index1_la <- newindex_la
        }
        indexall[, , al] <- newindex_la
        residall[, , al] <- newresid_la
      }
      else {
        IndexMatrix <- matrix(NA, nrow = h, ncol = (length(lambdas) - 
                                                      1))
        ResidMatrix <- matrix(NA, nrow = n, ncol = (length(lambdas) - 
                                                      1))
        for (la in 1:(length(lambdas) - 1)) {
          lambda <- lambdas[la + 1]
          newindex_la <- index1_la
          objbest <- tol
          cstep.mod <- CStep(x, y, family, newindex_la, 
                             h, hsize, alpha, lambda/h, scal)
          countloop <- 0
          while ((cstep.mod$object > objbest) & (countloop < 
                                                 csteps)) {
            countloop <- countloop + 1
            objbest <- cstep.mod$object
            newindex_la <- cstep.mod$index
            newresid_la <- cstep.mod$residu
            cstep.mod <- CStep(x, y, family, newindex_la, 
                               h, hsize, alpha, lambda/h, scal)
            index1_la <- newindex_la
          }
          IndexMatrix[, la] <- newindex_la
          ResidMatrix[, la] <- newresid_la
        }
        lambda <- lambdas[1]
        newindex_al <- index1_al
        objbest <- tol
        cstep.mod <- CStep(x, y, family, newindex_al, 
                           h, hsize, alpha, lambda/h, scal)
        countloop <- 0
        while ((cstep.mod$object > objbest) & (countloop < 
                                               csteps)) {
          countloop <- countloop + 1
          objbest <- cstep.mod$object
          newindex_al <- cstep.mod$index
          newresid_al <- cstep.mod$residu
          cstep.mod <- CStep(x, y, family, newindex_al, 
                             h, hsize, alpha, lambda/h, scal)
          index1_al <- newindex_al
        }
        IndexMatrix <- cbind(newindex_al, IndexMatrix)
        ResidMatrix <- cbind(newresid_al, ResidMatrix)
        indexall[, , al] <- IndexMatrix
        residall[, , al] <- ResidMatrix
      }
    }
  }
  return(list(indexall = indexall, residall = residall))
}