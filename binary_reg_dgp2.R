binary_reg_dgp2 <- function(n, beta, beta0 = 1, sigma_in = NULL, dirty = 0, type = "latent", v_outlier = "VW", test, outlier_mean = 20, outlier_sd = 1){
  
  ## binary_reg_dgp2(): generates data following KHF (2017) based on Bernoulli, parametrized by structural eta
  ## INPUTS: 
  # n: sample size
  # beta: the true coefficients including the zeros
  # beta0: the intercept, default = 1
  # sigma: the covariance matrices (?list of covs?)
  # dirty: proportion (decimals) of contamination (default = 0)
  # type: type of DGP: "latent" is the one based on latent factors as in KHF, "bernoulli" is the one based on Bernoulli parameter parametrized by eta
  # v_outlier: if additional vertical outliers are added, KHF (2017) advocate this approach, I think this does not make much sense
  # test: if a test set needs to be created, this just overrides all settings to have a clean dataset
  # outlier_mean: mean of the outlying N(mean , sd)-distributed outlying points, default is 20
  # outlier_sd: standard deviation of the N(mean , sd)-distributed outlying points, default is 1
  
  # setting seed? Stil not 100% sure on this, it will allow us to fully reconstruct based on the seed...
  
  # to do: load necessary packages for rmvnorm
  
  # Extracting dimensionalities from the beta vector
  p <- length(beta) # Total dimensionality WITHOUT INTERCEPT
  p_a <- sum(beta != 0) # Dimensionality INFORMATIVE part of the DGP
  p_b <- p - p_a # Dimensionality of the UNINFORMATIVE part (is actually not part of the DGP)
  
  ## Error Catching
  # Unadmissable type of DGP (auto to bernoulli)
  if(type != "bernoulli" & type != "latent"){
    warning("No correct type of DGP specified (bernoulli, latent), type set to bernoulli")
    type <- "bernoulli"
  }
  
  # Inadmissable dirty%
  if(dirty > 0.5 | dirty < 0){
    stop("The amount of contamination (dirty) must be set at min. 0 or max. 0.5")
  }
  # v_outlier
  if(v_outlier != "VW" & v_outlier != "KHF"){
    stop("The outlier type should be VW or KHF")
  }
  
  #### CLEAN Predictor (X) data
  ### Informative predictors (X_a)
  ## Covariance Matrix (sigma_a)
  # Default sigma (KHF paper, rho = 0.9 decreasing)
  if(is.null(sigma_in)){   
    
    # Initializing
    sigma_a <- matrix(NA, nrow = p_a, ncol = p_a)
    rho_a <- 0.9 # Standard used in the paper
    
    # Filling up with a series of powers of rho based on "distance" between predictors
    for(i in 1:p_a){
      for(j in 1:p_a){
        sigma_a[i, j] <- rho_a^(abs(i - j))
      }
    }
  }
  
  # If sigma is supplied
  if(!is.null(sigma_in)){
    sigma_a <- sigma_in
  }
  
  # Generating Gaussian predictor data using sigma_a
  X_a <- rmvnorm(n = n,
                 mean = rep(0, p_a),
                 sigma = sigma_a)
  
  ### Uninformative predictors (x_b)
  X_b <- NULL
  if(p_b > 0){ # In case there are only informative predictors: this will be skipped
    sigma_b <- matrix(NA, nrow = p_b, ncol = p_b)
    rho_b <- 0.5
    
    # Filling up with a series of powers of rho based on "distance" between predictors
    for(i in 1:p_b){
      for(j in 1:p_b){
        sigma_b[i, j] <- rho_b^(abs(i - j))
      }
    }
    X_b <- rmvnorm(n = n,
                   mean = rep(0, p_b),
                   sigma = sigma_b)
    
  }
  
  # Combining both X_a and X_b
  X <- cbind(X_a, X_b) # Can bind X_b even if NULL
  
  ## Creating CLEAN outcomes
  # Creating linear predictor (that is, including intercept, given by beta0)
  eta <- beta0 + X %*% beta
  
  # y: BERNOULLI DGP TYPE
  if(type == "bernoulli"){
  # Creating prob (pi or p) parameter
  prob <- exp(eta) / (1 + exp(eta))
  
  # y: Creating 0/1 outcomes
  y <- rbinom(n = n,
              size = 1,
              prob = prob)
  }
  
  # y: LATENT DGP
  if(type == "latent"){
    # Creating Gaussian error
    e <- rnorm(n = n, 
               mean = 0, 
               sd = 1)
    
    # Creating latent continuous outcome
    y_cont <- eta + e
    
    # Dichotomizing
    y <- ifelse(y_cont > 0, 
                yes = 1, 
                no = 0) # BINARY OUTCOME
  }
  
  # Saving proportion of 0/1
  n_0 <- sum(y == 0)
  n_1 <- sum(y == 1)
  
  ### CONTAMINATION (only in informative part!)
  if(dirty > 0){
    
    ## Leverage points (X-outliers)
    # Amount of observations to contaminate
    n_contam <- ceiling(dirty * n_0) # had this set to floor originally, but got into trouble with 0 :(
    
    # Contaminating using N(20, 1) (Different case if single predictor versus multiple predictors)
    if(p == 1){
      if(beta > 0){
        # If beta positive then add outliers far right at y = 0
        X_a[y == 0][1:n_contam] <- rnorm(n = n_contam, 
                                         mean = outlier_mean, 
                                         sd = outlier_sd)
      }
      if(beta < 0){
        # If beta negative then add outliers far right at y = 1 because Sigmoid is mirrored
        X_a[y == 1][1:n_contam] <- rnorm(n = n_contam, 
                                         mean = outlier_mean, 
                                         sd = outlier_sd)
      }
    } else if(p != 1){
      # Creating indicator variable that matches with y
      outlier_side <- ifelse(beta > 0,
                             yes = 0,
                             no  = 1)
      # Because, given N(20, 1)-outliers, we need them put at y == 0 if beta_j > 0 and at y == 1 if beta_j < 0
      
      for(i in 1:p_a){
        X_a[y == outlier_side[i], i][1:n_contam] <- rnorm(n = n_contam, 
                                                          mean = outlier_mean, 
                                                          sd = outlier_sd)
      }
    }
    
    ## Vertical outliers following KHF (I don't recommend this)
    if(v_outlier == "KHF"){
      y[y == 0][1:n_contam] <- 1 
    }
    
    ## Combining both
    X <- cbind(X_a, X_b)
  }
  
  # OUTPUT
  data_frame <- data.frame(y = y, X)
  return(data_frame)
}