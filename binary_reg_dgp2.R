binary_reg_dgp2 <- function(n, beta, beta0 = 1, sigma_in = NULL, rho_a = 0.9, dirty = 0, type = "latent", v_outlier = "VW", test, outlier_mean = 20, outlier_sd = 1, output_outlier_rows = FALSE, verbose=FALSE){
  
  ## binary_reg_dgp2(): generates data following KHF (2017) based on Bernoulli, parametrized by structural eta
  # NEEDS UPDATING (DOCS)
  
  ## INPUTS: 
  # n: sample size
  # beta: the true coefficients including the zeros
  # beta0: the intercept, default = 1, this is needed to create the linear predictor but IT IS NOT PART OF THE CREATED DATASET!
  # sigma_in: the covariance matrices (?list of covs?)
  # rho_a: the correlation coefficient between different informative covariates (exponentially decreasing by distance over columns)
  # dirty: proportion (decimals) of contamination (default = 0)
  # type: type of DGP: "latent" is the one based on latent factors as in KHF, "bernoulli" is the one based on Bernoulli parameter parametrized by eta
  # v_outlier: if additional vertical outliers are added, KHF (2017) advocate this approach, I think this does not make much sense
  # test: if a test set needs to be created, this just overrides all settings to have a clean dataset
  # outlier_mean: mean of the outlying N(mean , sd)-distributed outlying points, default is 20
  # outlier_sd: standard deviation of the N(mean , sd)-distributed outlying points, default is 1
  # outlier_rows: if set to TRUE, the function outputs a list with the row numbers of the outlying observations
  # verbose: makes some additional prints, handy for debugging
  
  # setting seed? Stil not 100% sure on this, it will allow us to fully reconstruct based on the seed...
  
  ## OUTPUT
  # A dataframe with data (EXCLUDING INTERCEPT), hence run functions for fitting with intercept = TRUE
  
  # Loading library for constructing multivariate normal distributions
  library(mvtnorm)
  
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
  
  # Inadmissable dirty proportions
  if(dirty > 0.5 | dirty < 0){
    stop("The amount of contamination (dirty) must be set at min. 0 or max. 0.5 (50%)")
  }
  # v_outlier
  if(v_outlier != "VW" & v_outlier != "KHF"){
    stop("The outlier type should be VW or KHF")
  }
  
  # Inadmissable informative dimension p_a:
  if(p_a < 1){
    stop("The useful dimensionality should be 1 or higher")
  }
  
  #if not (all(beta > 0) | all(beta < 0)){
  #  warning("")
  #}
  
  #### CLEAN Predictor (X) data
  ### Informative predictors (X_a)
  ## Covariance Matrix (sigma_a)
  # Default sigma (KHF paper, rho = 0.9 decreasing)
  if(is.null(sigma_in)){   
    
    # Initializing
    sigma_a <- matrix(NA, nrow = p_a, ncol = p_a)
    # Filling up with a series of powers of rho based on "distance" between predictors
    for(i in 1:p_a){
      for(j in 1:p_a){
        sigma_a[i, j] <- rho_a^(abs(i - j)) # Exponentially decreasing over the "column distance"
      }
    }
  }
  
  # If sigma is supplied
  if(!is.null(sigma_in)){
    sigma_a <- sigma_in
  }
  
  # Generating Gaussian predictor data using sigma_a
  X_a <- rmvnorm(n = n,
                 mean = rep(0, p_a), # 0-mean vector
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
  
  # y: Creating 0/1 outcomes, with the Bernouilli parameter being the expression with eta
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
  
  # Printing outcome distribution
  if(verbose) {
    print(paste0("The proportion of 0-outcomes is: ", round(n_0/n, 2)))
    print(paste0("The proportion of 1-outomces is: ", round(n_1/n, 2)))
  }
  
  if (n_0/n < 0.3 | n_0/n > 0.7){
    warning(paste0("Unbalanced classes, the proportion of 0-outcomes is: ", round(n_0/n, 2)))
  }
  
  ### CONTAMINATION (only in informative part!)
  if(dirty > 0){
    
    ## Leverage points (X-outliers)
    # Amount of observations to contaminate
    n_contam_1 <- ceiling(dirty * n_1) # had this set to floor originally, but got into trouble with 0 :(
    n_contam_0 <- ceiling(dirty * n_0)
    
    n_contam = n_contam_0 + n_contam_1 # TEMP
    
    # ERROR CATCH:
    if(n_contam >= n_0){
      stop("There are not enough 0 outcomes to contaminate")
    }
      
    
    ## Contaminating using another normal or multivariate normal (with vastly differing location)
    # Single covariate case (rather trivial)
    # TODO (FIX THIS TO INCORPORATE n_contam_0 and n_contam_1)
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
    # Multivariate case
    } else if(p != 1) {
      # Creating indicator variable that will tell us on which side an outlier should lie for that informative covariate
      outlier_side_1 <- ifelse(beta > 0,
                               yes = 0,
                               no  = 1)
      outlier_side_0 <- ifelse(beta > 0,
                               yes = 1,
                               no = 0)
      print(outlier_side_1)
      # EXPLANATION: Given e.g. N(20, 1) outliers, we need them placed at y == 0 if beta_j > 0 and at y == 1 if beta_j < 0
      # Visualize this for yourself if not clear in a setting with a single covariate.
      
      # Looping over each informative covariate
      for(i in 1:p_a){
        X_a[y == outlier_side_1[i], i][1:n_contam_1] <- rnorm(n = n_contam_1, 
                                                              mean = outlier_mean, 
                                                              sd = outlier_sd)
        
        X_a[y == outlier_side_0[i], i][1:n_contam_0] <- rnorm(n = n_contam_0,
                                                              mean = - outlier_mean,
                                                              sd = outlier_sd)
      }
    }
    
    # Gathering rows that have outliers
    outlier_rows_0 <- which(y == outlier_side_0)[1:n_contam_0]
    outlier_rows_1 <- which(y == outlier_side_1)[1:n_contam_1]
    outlier_rows <- sort(c(outlier_rows_0, outlier_rows_1))
    
    if(verbose) {
      print(n_contam_0)
      print(n_contam_1)
      print(paste0("The row numbers of the observations with outliers are: ", paste0(outlier_rows, collapse=" ")))
    }

    ## Vertical outliers following KHF (I don't recommend this, mainly for testing purposes!)
    if(v_outlier == "KHF"){
      print("You are using an outdated version of the DGP, there may be some issues with this option...")
      y[y == 0][1:n_contam] <- 1 
    }
    
    ## Combining both informative predictor matrix and uninformative one
    X <- cbind(X_a, X_b)
  }
  
  # Putting y, X into a dataframe
  data_frame <- data.frame(y = y, X)
  
  # Sorting dataframe on y to have the outliers be the first observations (because they are the 0 outcomes that should have been 1s)
  #data_frame <- data_frame[order(data_frame$y), ]
  
  # Resetting row names to avoid confusion
  #rownames(data_frame) <- NULL
  
  # OUTPUT
  if(output_outlier_rows) {
    return(list(data_frame, outlier_rows))
  } else {
    return(data_frame)
  }
}