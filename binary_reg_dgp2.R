binary_reg_dgp2 <- function(n, beta, beta0 = 1, sigma_in = NULL, rho_a = 0.9, X_a_mean = 0, rho_b = 0.5, dirty = 0, type = "latent", v_outlier = "VW", test, outlier_mean = 20, outlier_sd = 1, output_probs = FALSE, output_outlier_rows = FALSE, verbose=FALSE) {
  
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
  
  # TODO: MAINLY TESTED NOW FOR BERNOULLI TYPE OF DGP
  
  ## OUTPUT
  # A dataframe with data (EXCLUDING INTERCEPT), hence run functions for fitting with intercept = TRUE
  
  ## The idea is to first make clean data and then contaminate parts of this clean data.
  
  
  ### PROCESS FLOW
  ## .... TODO
  
  
  # Loading library for constructing multivariate normal distributions
  library(mvtnorm)
  
  # Extracting dimensionalities from the beta vector
  p <- length(beta) # Total dimensionality WITHOUT INTERCEPT
  p_a <- sum(beta != 0) # Dimensionality INFORMATIVE part of the DGP
  p_b <- p - p_a # Dimensionality of the UNINFORMATIVE part (is actually not part of the DGP)
  
  ## Error Catching
  # Unadmissable type of DGP (auto to bernoulli)
  if (type != "bernoulli" & type != "latent") {
    warning("No correct type of DGP specified (bernoulli, latent), type forced to bernoulli")
    type <- "bernoulli"
  }
  # Inadmissable dirty proportions
  if (dirty > 0.5 | dirty < 0) {
    stop("The amount of contamination (dirty) must be set at min. 0 or max. 0.5 (50%)")
  }
  # v_outlier
  if (v_outlier != "VW" & v_outlier != "KHF") {
    warning("The vertical outlier type should be VW or KHF, v_outlier forced to VW method")
  }
  
  # Inadmissable informative dimension p_a:
  if (p_a < 1) {
    stop("The useful dimensionality should be 1 or higher!")
  }
  
  # All informative betas should have value != 0
  if (any(which(beta == 0) < p_a)) {
    stop("The informative part of beta seems to contain true zeros, this is not allowed currently!")
  }
  
  #### CLEAN Predictor (X) data
  ### Informative predictors (X_a)
  ## Covariance Matrix (sigma_a)
  # Default sigma format (KHF paper, rho = 0.9, decreasing)
  if (is.null(sigma_in)) {   
    sigma_a <- matrix(NA, nrow = p_a, ncol = p_a) # Initializing
    for(i in 1:p_a){
      for(j in 1:p_a){
        sigma_a[i, j] <- rho_a^(abs(i - j)) # Exponentially decreasing over the "column distance"
        # Explanation: x1-x2 have rho_a corr, x1 - x3 have corr rho_a^2, ... 
      }
    }
  }
  
  # If sigma is supplied, pass through
  if(!is.null(sigma_in)){
    sigma_a <- sigma_in
  }
  
  # Generating Gaussian predictor data using sigma_a
  X_a_mean <- rep(X_a_mean, times = p_a) # Making mean scalar into a vector!
  X_a <- rmvnorm(n = n,
                 mean = X_a_mean, # DEFAULT: 0-vector
                 sigma = sigma_a)
  
  ### Uninformative predictors (x_b)
  X_b <- NULL
  if (p_b > 0) { # In case there are only informative predictors: this will be skipped
    sigma_b <- matrix(NA, nrow = p_b, ncol = p_b)  # Initializing
    for (i in 1:p_b) {
      for (j in 1:p_b) {
        sigma_b[i, j] <- rho_b^(abs(i - j)) # Filling up with a series of powers of rho based on "distance" between predictors
        # Explanation: x1-x2 have rho_b corr, x1 - x3 have corr rho_b^2, ... 
      }
    }
    
    # Generating Gaussian uninformative variables using sigma_b
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
    # Creating prob (pi, p or prob) parameter
    prob <- exp(eta) / (1 + exp(eta))
    
    # y: Creating 0/1 outcomes, with the Bernouilli parameter being parametrized the expression with eta
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
    
    # Creating latent (=continuous) outcome
    y_cont <- eta + e
    
    # Dichotomizing to binary outcome
    y <- ifelse(y_cont > 0, 
                yes = 1, 
                no = 0)
  }
  
  # Saving proportion of 0/1 ()
  n_0 <- sum(y == 0)
  n_1 <- sum(y == 1)
  
  # Printing outcome distribution
  if(verbose) {
    print(paste0("The proportion of 0-outcomes is: ", round(n_0/n, 2)))
    print(paste0("The proportion of 1-outomces is: ", round(n_1/n, 2)))
  }
  
  # Checking if distribution too imbalanced: printing warning
  if (n_0/n < 0.3 | n_0/n > 0.7){
    warning(paste0("Unbalanced classes, the proportion of 0-outcomes is: ", round(n_0/n, 2)))
  }
  
  ### CONTAMINATION (only in informative part!)
  if(dirty > 0) {
    
    ## TODO: Could be two options: you just change some 1 -> 0 and vice versa or you take some extremes
    # Or you could just do it similarly to 
    
    ## Leverage points (X-outliers)
    # Calculating amount of observations to contaminate
    n_contam_0 <- ceiling(dirty * n_0) # Amount of (truly) 0 outcomes that need to become 1 as outlier
    n_contam_1 <- ceiling(dirty * n_1) # Amount of (truly) 1 outcomes that need to become 0 as outlier
    n_contam = n_contam_0 + n_contam_1
    # Used floor() first instead of ceiling() but gave issues...
    
    # If too many contaminated data points... # TODO CHECK FROM HERE
    if(n_contam >= n_0 | n_contam >= n_1){
      error("There are not enough 0 or 1 outcomes to contaminate, the outcome distribution might be too imbalanced") # Just a check, unlikely, unless in skewed situation
    }
    
    #### REPLACING THE COVARIATE VALUES WITH 
    if (v_outlier == "VW") {
      ### Single covariate case (rather trivial)
      # TODO (FIX THIS TO INCORPORATE n_contam_0 and n_contam_1) i.e. both 0 and 1 outliers
      if (p == 1) {
        ## Case: positive beta1
        if (beta > 0) {
          # Replacing predictor value by something too high for the original 0-outcome
          X_a[y == 0][1:n_contam_0] <- rnorm(n = n_contam_0, 
                                             mean = X_a_mean + outlier_mean, 
                                             sd = outlier_sd)
          
          # Replacing predictor value by something too low for the original 1-outcome
          X_a[y == 1][1:n_contam] <- rnorm(n = n_contam_1,
                                           mean = X_a_mean - outlier_mean,  # TAKING NEGATIVE
                                           sd = outlier_sd)
          
          ## Case: negative beta1 (mirrored situation)
        } else if (beta < 0) {
          # Replacing predictor value by something too high for the original 0-outcome
          X_a[y == 0][1:n_contam_0] <- rnorm(n = n_contam_0,
                                             mean = X_a_mean - outlier_mean, # TAKING NEGATIVE
                                             sd = outlier_sd) 
          
          # Replacing predictor value by something too low for the original 1-outcome
          X_a[y == 1][1:n_contam_1] <- rnorm(n = n_contam_1, 
                                             mean = X_a_mean + outlier_mean,
                                             sd = outlier_sd)
        }
        
        ### Multivariate case
      } else if(p > 1) {
        # Separating beta into informative and uninformative part
        beta_informative <- beta[beta != 0]
        
        # Making an indicator of positive beta (Positive: +1 / Negative: -1) for later use when placing outliers
        positive_beta_indicator <- ifelse(beta_informative > 0,
                                          yes = 1,
                                          no = - 1)
        
        ## Looping over each informative predictor
        for (i in 1:p_a) {
          # Extracting current indicator: for positive beta: 1 / for negative beta: -1
          indicator <- positive_beta_indicator[i]
          X_a[y == 0, i][1:n_contam_0] <- rnorm(n = n_contam_0,
                                                mean = X_a_mean + (indicator * outlier_mean),
                                                sd = outlier_sd)
          
          X_a[y == 1, i][1:n_contam_1] <- rnorm(n = n_contam_1,
                                                mean = X_a_mean - (indicator * outlier_mean), 
                                                sd = outlier_sd)
        }
      }
      
      # Gathering row numbers that have outliers
      outlier_rows_0 <- which(y == 0)[1:n_contam_0]
      outlier_rows_1 <- which(y == 1)[1:n_contam_1]
      outlier_rows <- sort(c(outlier_rows_0, outlier_rows_1))
      
      # Create n-vector boolean mask based on the outlier row numbers
      mask <- as.matrix(rep(FALSE, length = n)) # Forcing as matrix good for naming (vs. vectors :( )
      mask[outlier_rows] <- TRUE
      
      # Additional prints (debugging, or interested)
      if (verbose) {
        print(paste0("The amount of outliers in the 0-outcomes is: ", n_contam_0))
        print(paste0("The amount of outliers in the 1-outcomes is: ", n_contam_1))
        print(paste0("The row numbers of the observations with outliers are: ", paste0(outlier_rows, collapse=" ")))
      }
    }
    
    ## Vertical outliers following KHF (I don't recommend this, mainly for testing purposes!)
    if (v_outlier == "KHF") {
      warning("You are using an outdated version of the DGP, there may be some issues with this option...")
      y[y == 0][1:n_contam] <- 1 
    }
    
    ## Combining both informative predictor matrix and uninformative one
    X <- cbind(X_a, X_b)
  }
  
  # Putting y, X into a dataframe
  data_frame <- data.frame(y = y, X)
  
  
  # OUTPUT
  if (output_outlier_rows) {
    output <- list(data_frame, outlier_rows)
    if (output_probs) {
      output <- list("data" = data_frame, "outliers" = mask, "probs" = prob)
    }
  } else if (output_probs) {
    output <- list("data" = data_frame, "probs" = prob)
    
  }  else {
    output <- data_frame
    
  }
  # RETURN OUTPUT
  return(output)
}