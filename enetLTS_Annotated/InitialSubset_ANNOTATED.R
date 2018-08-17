InitialSubset_ANNOTATED <- function (x, y, family, h, hsize, alpha, lambda, nsamp, scal, 
                                     para, ncores, seed) 
{
  if (!is.null(seed)) 
    set.seed(seed)
  
  ## Creation of 500 elemental subsets
  # Binomial Case: Subsets of size 4 (2 for each outcome)
  if (family == "binomial") {
    index.subsets <- replicate(nsamp, 
                               c(sample(which(y == 1), 2), 
                                 sample(which(y == 0), 2)))
    # Sampling of ROWNUMBERS (indices), 2 for 0-outcomes and 2 for 1 outcomes
    # The result, index.subsets is a MATRIX (4 x nsamp) default: 4 x 500)
    # HENCE: EACH COLUMN IS AN ELEMENTAL SUBSET OF SIZE 4 (SEE LATER *)
    
    # Gaussian Case: Subsets of size 3 (outcome doesn't matter here)
  } else if (family == "gaussian") {
    index.subsets <- replicate(nsamp, 
                               sample.int(nrow(x), 3)) # Gaussian: only 3 observations needed
  }
  
  ############################ DEFINING twoCstep() #####################################################
  
  twoCstep <- function(c, x, y, family, h, hsize, alpha, lambda) {

    ## Running FIRST C-Step
    # Binomial Case
    if (family == "binomial") {
      Cstep1 <- enetLTS:::CStep(x = x, 
                                y = y, 
                                family = family, 
                                indx = index.subsets[, c],  # WHERE DOES c come from? -> IT'S THE 1:nsamp that is supplied in the lapply loop (SEE EARLIER ALSO *)
                                h = h, 
                                hsize = hsize, 
                                alpha = alpha, 
                                lambda = lambda/4, # Rescaling lambda because 4 observations in the elemntal subset
                                scal = FALSE)
      # CStep() function call returns: object (objective function), index (indices of length = h), residu (residuals of size n!!!), beta (of the fit used in the C-step)
      
      # Gaussian Case
    } else if (family == "gaussian") {
      Cstep1 <- enetLTS:::CStep(x = x, 
                                y = y, 
                                family = family, 
                                indx = index.subsets[, c], # WHERE DOES c come from? -> IT'S THE 1:nsamp that is supplied in the lapply loop
                                h = h, 
                                hsize = hsize, 
                                alpha = alpha, 
                                lambda =lambda/3, # Rescaling lambda because 3 observations in the elemental subset
                                scal = FALSE)
    }
    # Assigning results offirst C-step
    indx1 <- Cstep1$index # SERVES AS INPUT FOR C-STEP 2
    object1 <- Cstep1$object 
    
    ## Running SECOND C-Step (using the results from the first one!)
    Cstep2 <- enetLTS:::CStep(x = x, 
                              y = y, 
                              family = family, 
                              indx = indx1,  # So now we supply them with a new observation index set of size h, because that is was a CStep() outputs
                              h = h, 
                              hsize = hsize, 
                              alpha = alpha, 
                              lambda = lambda/h, 
                              scal = scal)
    # Assigning Results of second C-step
    indx2 <- Cstep2$index
    object <- Cstep2$object
    
    # OUTPUT (results of SECOND C-step are assigned)
    return(list(obj = object, indx = indx2))
  } # END OF twoCstep() FUNCTION
  
  #################################################################################
  
  # Case: parallel computing is set to TRUE
  if (para) {
    subsets <- parallel:::mclapply(1:nsamp, FUN = twoCstep, x = x, y = y,  # YES NOTE THAT c is indeed supplied as the first argument here implicitly (1:nsamp)
                                   family = family, h = h, hsize = hsize, alpha = alpha, 
                                   lambda = lambda, mc.cores = ncores)
    # YES THIS GIVES WARNIGNS ABOUT LOGNET WITH LESS THAN 8 OBS!
    
    # Case: nonparallel
  } else {
    # Baiscally what is done: supply a COLUMN of index.subsets (= MATRIX) to the twoCstep() function
    subsets <- lapply(1:nsamp, FUN = twoCstep, x = x, y = y, 
                      family = family, h = h, hsize = hsize, alpha = alpha, 
                      lambda = lambda)
    # YES THIS GIVES WARNIGNS ABOUT LOGNET WITH LESS THAN 8 OBS! (Because it calls CStep() with the INITIAL SUBSETS!)
    
    # subsets object is a list of length nsamp (default: 500) with $obj with the value of the objective function (?) and $indx of length h with the index numbers!
  }
  return(list(subsets = subsets, index.subsets = index.subsets))
}