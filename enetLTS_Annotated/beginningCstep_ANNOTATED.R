function (x, y, family, h, hsize, alpha, lambda, nsamp, s1, ncores, 
          csteps, tol, scal, para, seed) 
{
  
  ### beginningCstep() FUNCTION
  ## This function is called by warmCsteps() with following inputs: x = x, y = y, family = family, h = h, hsize = hsize,
  # alpha = alpha, lambda = lambda, nsamp = nsamp, s1 = s1, ncores = ncores, csteps = csteps, tol = tol, scal = scal, para = para, seed = seed
  
  ## From the KHF (2017) paper: this beginning Cstep is only done for the first combination of alpha, lambda
  
  # Selects best 10 
  H2 <- selectbest10(x, y, family, h, hsize, alpha, lambda, 
                     nsamp, s1, para, ncores, scal, seed)
  if (para) {
    lastbestindex <- mclapply(1:s1, function(zz, x, y, family, 
                                             h, hsize, alpha, lambda, H2) {
      indexsubbest <- H2$idxbest[[zz]]
      objbest <- tol
      cstep.mod <- CStep(x, y, family, indexsubbest, h, 
                         hsize, alpha, lambda/h, scal)
      countloop <- 0
      while ((cstep.mod$object > objbest) & (countloop < 
                                             csteps)) {
        countloop <- countloop + 1
        objbest <- cstep.mod$object
        newindex <- cstep.mod$index
        beta <- cstep.mod$beta
        cstep.mod <- CStep(x, y, family, newindex, h, 
                           hsize, alpha, lambda/h, scal)
      }
      return(list(lastindex = newindex, objbest = objbest, 
                  countloop = countloop, residu = cstep.mod$residu, 
                  beta = beta))
    }, x, y, family, h, hsize, alpha, lambda, H2, mc.cores = ncores)
  }
  else {
    lastbestindex <- lapply(1:s1, function(zz, x, y, family, 
                                           h, hsize, alpha, lambda, H2) {
      indexsubbest <- H2$idxbest[[zz]]
      objbest <- tol
      cstep.mod <- CStep(x, y, family, indexsubbest, h, 
                         hsize, alpha, lambda/h, scal)
      countloop <- 0
      while ((cstep.mod$object > objbest) & (countloop < 
                                             csteps)) {
        countloop <- countloop + 1
        objbest <- cstep.mod$object
        newindex <- cstep.mod$index
        beta <- cstep.mod$beta
        cstep.mod <- CStep(x, y, family, newindex, h, 
                           hsize, alpha, lambda/h, scal)
      }
      return(list(lastindex = newindex, objbest = objbest, 
                  countloop = countloop, residu = cstep.mod$residu, 
                  beta = beta))
    }, x, y, family, h, hsize, alpha, lambda, H2)
  }
  obj <- NULL
  for (i in 1:s1) {
    obj <- c(obj, lastbestindex[[i]]$objbest)
  }
  whichbestindex <- sort(obj, decreasing = TRUE, index.return = TRUE)$ix[1]
  index <- lastbestindex[[whichbestindex]]$lastindex
  resid <- lastbestindex[[whichbestindex]]$residu
  return(list(index = index, resid = drop(resid)))
}
<bytecode: 0x000000008518ba98>
  <environment: namespace:enetLTS>