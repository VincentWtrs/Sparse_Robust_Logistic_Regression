binary_reg_dgp2 <- function(n, p, p_a, beta, sigma_in = NULL, dirty = 0, type = "latent", v_outlier = 0, test){
  ## binary_reg_dgp1(): generates data following KHF (2017) based on Bernoulli, parametrized by structural eta
  ## INPUTS: n: sample size
  # p: dimensionality EXCLUDING INTERCEPT!
  # p_a: dimensionality of the informative part (excluding intercept)
  # sigma: the covariance matrices (?list of covs?)
  # dirty: proportion (decimals) of contamination
  ## NOTE: INTERCEPT ASSUMED TO BE 1!
  
  # QUID: setting seed?
  
  p_b <- p - p_a # Dimensionality of the uninformative part
  
  ## Error Catching
  # Unadmissable dimensionality
  if(p_a > p){
    stop("The amount of nonzero coefficients (p_a) must be equal or smaller than the total dimensionality (p)")
  }
  
  # Unadmissable type of DGP (auto to bernoulli)
  if(type != "bernoulli" & type != "latent"){
    warning("No correct type of DGP specified (bernoulli, latent), type set to bernoulli")
    type <- "bernoulli"
  }
  
  # Inadmissable dirty%
  if(dirty > 0.5 | dirty < 0){
    stop("The amount of contamination (dirty) must be set at min. 0 or max. 0.5")
  }
  
  # Nonconformable beta with p_a
  if(length(beta) != p){
    stop("The dimensionality of the supplied beta vector is not equal to that of true supplied dimensionality (p_a)")
  }
  
  # v_outlier
  if(v_outlier != 0 & v_outlier != 1){
    stop("The outlier type should be 0 or 1 (= KHF method)")
  }
  
  #### CLEAN Predictor (X) data
  ### Informative predictors (X_a)
  ## Covariance Matrix (sigma_a)
  # Default sigma (KHF paper)
  if(is.null(sigma_in)){   
    
    # Initializing
    sigma_a <- matrix(NA, nrow = p_a, ncol = p_a)
    rho_a <- 0.9
    
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
  # Creating linear predictor (that is, including intercept, set to 1)
  eta <- 1 + X %*% beta
  
  # y: BERNOUILLI DGP
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
    e <- rnorm(n = n, mean = 0, sd = 1)
    
    # Creating latent continuous outcome
    y_cont <- eta + e
    
    # Dichotomizing
    y <- ifelse(y_cont > 0, yes = 1, no = 0) # BINARY OUTCOME
  }
  
  # Saving proportion of 0/1
  n_0 <- sum(y == 0)
  n_1 <- sum(y == 1)
  
  ### CONTAMINATION (only in informative part!)
  if(dirty > 0){
    
    ## Leverage points (X-outliers)
    # Amount of observations to contaminate
    n_contam <- ceiling(dirty * n_0) # had this set to floor originally, but got into trouble with 0 :(
    
    # Contaminating using N(20, 1)
    if(p == 1){
      X_a[y == 0][1:n_contam] <- rnorm(n = n_contam, mean = 20, sd = 1)
    } else if(p != 1){
      
      X_a[y == 0, ][1:n_contam, ] <- rmvnorm(n = n_contam, mean = rep(20, p_a))
      # (sigma default hence I-matrix (independent))
      # , drop= FALSE to maintain a matrix
    }
    
    ## Vertical outliers
    if(v_outlier == 1){
      y[y == 0][1:n_contam] <- 1
    }
    
    ## Combining both
    X <- cbind(X_a, X_b)
    
  }
  
  # OUTPUT
  data_frame <- data.frame(y = y, X)
  return(data_frame)
}