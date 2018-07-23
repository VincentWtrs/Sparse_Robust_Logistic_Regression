### logit_sim_summary() function: to take the results from a logit_sim object and return the relevant statistics
## The main idea is that it will summarize over the runs, i.e. per estimator
## INPUT:
# logit_sim: logit_sim object from the logit_sim function, a list containing model information, ...
# NOTE: don't mix up beta dimensionalities here, this depends on what comes out of the beta
logit_sim_summary <- function(logit_sim){
  
  # Extracting/defining quantities from logit_sim input
  l <- length(logit_sim) # Amount of estimators used
  estimator_names <- names(logit_sim) # Getting names of estimators used
  
  runs <- logit_sim[[1]]$runs # Amount of runs done
  
  beta <- c(1, logit_sim[[1]]$beta) # Adding intercept (fixed at 1 by binary_dgp atm) # THIS CAN CHANGE!!!
  beta_dim <- length(beta) # beta dim (INCL. INTERCEPT)
  p_a <- sum(beta != 0) # True nonzero (Useful)
  p_b <- sum(beta == 0) # True zero
  
  # Initiating output
  output <- vector("list", length = l)
  
  output_stats <- vector("list", length = l)
  output2 <- vector("list", length = l)
  

  ### For each estimator
  for(i in 1:l){
    ## Initializing some variables
    coefs <- data.frame(matrix(NA, nrow = runs, ncol = beta_dim)) # Each row = run, each column = variable
    names(coefs) <-  c("(Intercept)", paste0("X", 1:(p-1))) # Giving column names
    
    fpr <- numeric(length = runs)
    fnr <- numeric(length = runs)
    precision <- numeric(length = runs)
    tables <- vector("list", length = runs)
    prop0 <- numeric(length = runs)
    

    ## For each run
    for(r in 1:runs){
      
      # Extracting current (run) model
      model <- logit_sim[[i]]$models[[r]] # Model i, run r
      
      # Extracting Coefficients
      coefs_current <- coef(model)
      coefs[i, ] <- coefs_current
      
  
      # Precision (RMSE of betas) DIFFERENT DEF BY TAKING MEAN TO NORMALIZE FOR p
      precision[r] <- sqrt(mean(coefs_current - beta)^2)
                           
      # False Positive Rate (FPR, lower is better)
      fpr[r] <- sum(coefs_current[(p_a + 1):beta_dim] != 0) / p_b

      # False Negative Rate (FNR, lower is better)
      fnr[r] <- sum(coefs_current[1:p_a] == 0) / p_a
      
      # Outcome (y) class frequency
      tables[[r]] <- table(logit_sim[[i]]$train[[r]]$y) # Freq table
      prop0[r] <- tables[[r]][1] / sum(tables[[r]]) # Getting proportion of 0
      # At the end we take summary of it to an idea
    }
    # Calculating averages
    avg_precision <- mean(precision)
    avg_fpr <- mean(fpr)
    avg_fnr <- mean(fnr)
    avg_prop0 <- mean(prop0)
    
    # Some other variables
    freq <- summary(prop0)
    
    # Putting it all in a one list
    summary <- list(precision = precision,
                    avg_precision = avg_precision,
                    fpr = fpr,
                    avg_fpr = avg_fpr,
                    fnr = fnr,
                    avg_fnr = avg_fnr,
                    prop0 = prop0,
                    avg_prop0 = avg_prop0)
    
    # The pure statistics to be computed
    stats <- list(precision, fpr, fnr, prop0)
    stats_names <- c("precision", "fpr", "fnr", "prop0")
    
    # Putting list per model into big list
    output[[i]] <- summary
    # For our new output structure
    output_stats[[i]] <- stats
  }
  # Giving names of estimators to output
  names(output) <- estimator_names
  
  # OUTPUT
  #return(output)
  
  
  # For each stat
  for(i in 1:length(stats)){
    # Init
    output2[[i]] <- data.frame(matrix(NA, nrow = runs, ncol = l))
    names(output2[[i]]) <- estimator_names # Giving estimator names to columns of each df
    # For each estimator
    for(j in 1:l){
      output2[[i]][, j] <- output_stats[[j]][[i]]
    }
  }
  names(output2) <- stats_names
  
  stats_avg <- vector("list", length = length(stats))
  names(stats_avg) <- paste0(stats_names, "_avg")
  for(i in 1:length(stats)){
    stats_avg[[i]] <- colMeans(output2[[i]])
  }
  output2_complete <- c(output2, stats_avg)
  
  
  return(output2_complete)
}

#### Okay think about it. The input is a logit_sim object that actually already contains a lot of things. Now we want an object. Ahaa now it should basically be defined per 
# but it should actually be able to work easily with an object with multiple fits together. So it is a list with: # levels per estimator 
# then for each estimator: the summary and the plots. There should be a plotting indicator and also it should be easy to plot the figures with the data later on.


# Rethinking the format, maybe at the end we should actually do the whole thing per topic e.g. precision, fpr, etc. because you want small dataframes for easy plotting
temp <- logit_sim_summary(sim2)

# Now that's very easy to plot!
boxplot(temp$precision)
boxplot(temp$fpr)
boxplot(temp$fnr)




