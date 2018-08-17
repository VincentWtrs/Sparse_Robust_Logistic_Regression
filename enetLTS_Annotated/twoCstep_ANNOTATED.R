#### IGNORE THIS /// BETTER TO KEEP IT IN THE ORIGINAL FUNCTION (InitialSubset())



## Defining twoCstep() function
twoCstep <- function(c, x, y, family, h, hsize, alpha, lambda) {
  # I'M NOT CALLING IT "twoCstep_ANNOTATED" BECAUSE IT'S CALLED FROM ANOTHER FUNCTION
  
  ## Running FIRST C-Step
  # Binomial Case
  if (family == "binomial") {
    Cstep1 <- enetLTS:::CStep(x, 
                              y, 
                              family, 
                              index.subsets[, c],  # WHERE DOES c come from? -> from selectbest10
                              h, 
                              hsize, 
                              alpha, 
                              lambda/4, # Rescaling lambda because 4 observations in the elemntal subset
                              scal = FALSE)
    
    # Gaussian Case
  } else if (family == "gaussian") {
    Cstep1 <- enetLTS:::CStep(x, 
                              y, 
                              family, 
                              index.subsets[, c], 
                              h, 
                              hsize, 
                              alpha, 
                              lambda/3, # Rescaling lambda because 3 observations in the elemental subset
                              scal = FALSE)
  }
  # Assigning results of first C-step
  indx1 <- Cstep1$index
  object1 <- Cstep1$object
  
  ## Running SECOND C-Step (using the results from the first one!)
  Cstep2 <- CStep(x, 
                  y, 
                  family, 
                  indx1, 
                  h, 
                  hsize, 
                  alpha, 
                  lambda/h, 
                  scal)
  # Assigning Results of second C-step
  indx2 <- Cstep2$index
  object <- Cstep2$object
  
  # OUTPUT (results of SECOND C-step are assigned)
  return(list(obj = object, indx = indx2))
} # END OF twoCstep() FUNCTION