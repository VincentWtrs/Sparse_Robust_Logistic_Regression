### logit_sim_summary() function: to take the results from a logit_sim object and return the relevant statistics
## The main idea is that it will summarize over the runs, i.e. per estimator
## INPUT:
# logit_sim: A LIST OF logit_sim objects from the logit_sim function, a list containing model information, ...
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
  stats <- vector("list", length = l)
  
  # Defining the names of summary stats to be calculated (excluding their overall averages)
  stats_names <- c("loss", "misclass", "prec", "tpr", "fpr", "tnr", "fnr", "prop0_train", "prop0_test")
  
  ### For each estimator
  for(i in 1:l){
    ## Initializing
    # Coefficients
    coefs <- data.frame(matrix(NA, nrow = runs, ncol = beta_dim)) # Each row = run, each column = variable
    names(coefs) <-  c("(Intercept)", paste0("X", 1:(p-1))) # Giving column names
    
    tables <- vector("list", length = runs)
    
    # Initiating statistics
    loss <- numeric(length = runs) # Avg. Neg. Loglik loss
    misclass <- numeric(length = runs) # Misclass. rate
    precision <- numeric(length = runs) # Avg. RMSE on betas
    tpr <- numeric(length = runs) # True Positive Ratio
    fpr <- numeric(length = runs) # False Positive Ratio
    tnr <- numeric(length = runs) # True Negative Ratio
    fnr <- numeric(length = runs) # False Negative Ratio
    prop0_train <- numeric(length = runs) # The proportion of 0-outcomes (class imbalance)
    prop0_test <- numeric(length = runs) # Same, for test set

    ## For each run
    for(r in 1:runs){
      
      # Extracting model of current run
      model <- logit_sim[[i]]$models[[r]] # Estimator i, run r
      
      # Extracting Coefficients
      coefs_current <- coef(model) # For convenience, later on
      coefs[i, ] <- coefs_current
      
      # Loss
      loss[r] <- logit_sim[[i]]$loss[[r]]$avg_loss
      
      # Misclass 
      misclass[r] <- mean(logit_sim[[i]]$misclass[[r]])
      
      # Precision (Avg. RMSE of betas) DIFFERENT DEF BY TAKING MEAN TO NORMALIZE FOR p
      precision[r] <- sqrt(mean(coefs_current - beta)^2)
      
      # TPR: True Positive Rate (higher is better)
      tpr[r] <- sum(coefs_current[1:p_a] != 0) / p_a
                           
      # FPR: False Positive Rate (lower is better)
      fpr[r] <- sum(coefs_current[(p_a + 1):beta_dim] != 0) / p_b
      
      # TNR: True Negative Rate (higher is better)
      tnr[r] <- sum(coefs_current[(p_a + 1):beta_dim] == 0) / p_b

      # FNR: False Negative Rate (lower is better)
      fnr[r] <- sum(coefs_current[1:p_a] == 0) / p_a
      
      # Outcome (y) class frequency
      table_train <- table(logit_sim[[i]]$train[[r]]$y) # For training sets
      table_test <- table(logit_sim[[i]]$test[[r]]$y) # For test sets
      
      prop0_train[r] <- table_train[1] / sum(table_train)
      prop0_test[r] <- table_test[1] / sum(table_test)
    }
    
    # Putting it all in a one list
    stats[[i]]  <- list(loss,
                        misclass,
                        precision,
                        tpr,
                        fpr,
                        tnr,
                        fnr,
                        prop0_train,
                        prop0_test)
    names(stats[[i]]) <- stats_names
  }
  # Giving appropriate names to stats
  names(stats) <- estimator_names 
  # Now we have a list per estimator, for each estimator a numeric giving the stats
  
  # We'd like to get the output: list of dfs, each stats 2 df (one for stat, one for avg) with l columns (one for each estim.)
  
  ## OUTPUT
  output <- vector("list", length = length(stats_names))
  stats_avg <- vector("list", length = length(stats_names))
  names(stats_avg) <- paste0(stats_names, "_avg")
  
  # For each different stat
  for(i in 1:length(stats_names)){
    
    # Initializing a df with row for each run and col for each estimator
    output[[i]] <- data.frame(matrix(NA, nrow = runs, ncol = l))
    names(output[[i]]) <- estimator_names
    names(output) <- stats_names
    
    # For each estimator
    for(j in 1:l){
      output[[i]][, j] <- stats[[j]][[i]] # Col j (estim j) <- stats of estim j, stats i
      
      # Taking averages of the dfs we just have made
      stats_avg[[i]] <- colMeans(output[[i]])
      
    }
  }
  output <- c(output, stats_avg)
  return(output)
}