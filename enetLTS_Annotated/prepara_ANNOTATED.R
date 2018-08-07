enetLTS:::prepara <- function (x, y, family, index = NULL, robu = NULL) 
{
  ### DATA PREPARATION FUNCTION 
  ## This function is called by the main function enetLTS() with: x = xx, y = yy, family = family (binomial/gaussian), 
  # index not supplied and robu = 1
  
  # Setting robu = 0 if no value supplied
  if (is.null(robu)) 
    robu = 0
  
  ### If no index is supplied:
  if (is.null(index)) {
    ## No index AND robu > 0 (probably == 1)
    if (robu > 0) {
      # Binomial Case
      if (family == "binomial") {
        muy = y # Robust location estimate of = y itself
      }
      # Gaussian case
      else if (family == "gaussian") { # Robust location estimate of y = is med
        muy <- median(y)
      }
      # For the Covariates (X)
      mux <- apply(x, 2, median) # Location estimate: median for all columns (apply over cols)
      sigx <- apply(x, 2, mad) # Scale estimate: MAD for all colimns (apply over cols)
    }
    
    ## No index AND robu == 0 (classical location/scale estimates)
    else {
      if (family == "binomial") {
        muy = y
      }
      else if (family == "gaussian") {
        muy <- mean(y)
      }
      mux <- apply(x, 2, mean)
      sigx <- apply(x, 2, sd)
    }
  }
  
  ## Index is supplied (same but only do on indexed observations)
  else {
    if (robu > 0) {
      if (family == "binomial") {
        muy = y
      }
      else if (family == "gaussian") {
        muy <- median(y[index])
      }
      mux <- apply(x[index, ], 2, median)
      sigx <- apply(x[index, ], 2, mad)
    }
    else {
      if (family == "binomial") {
        muy = y
      }
      else if (family == "gaussian") {
        muy <- mean(y[index])
      }
      mux <- apply(x[index, ], 2, mean)
      sigx <- apply(x[index, ], 2, sd)
    }
  }
  
  # Normalizing Covariates (X) with the calculated location/scale
  xnor <- scale(x, mux, sigx) # Same for Gaussian or Binomial
  
  ## Normalizing outcome y
  # Binomial case
  if (family == "binomial") {
    ycen <- y # NOTHING IS DONE OBVIOUSLY
  }
  # Gaussian case
  else if (family == "gaussian") {
    ycen <- scale(y, muy, FALSE) # Only centering, not scaling (scale = FALSE)
  }
  
  # OUTPUT
  return(list(xnor = xnor, ycen = ycen, mux = mux, sigx = sigx, 
              muy = muy))
}
<environment: namespace:enetLTS>