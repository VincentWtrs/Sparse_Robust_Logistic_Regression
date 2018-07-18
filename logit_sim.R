logit_sim <- function(beta, sigma_in = NULL, p, p_a, n, runs, seed = 1234, dirty = 0, type, dgp = "latent"){
  ### LOGIT SIMULATION FUNCTION
  ## INPUT: beta_vector: a vector of true betas
  # n: the sample size (integer)
  # runs: amount of simulation runs (integer)
  # seed: a seed number for being able to reproduce
  # dirty: proportion of outlying cases (by default no outliers, decimal)
  # type: type of estimator
  # dgp: latent or bernoulli DGP from binary_reg_dgp
  ## OUPUT: A list of models, one for each run
  
  
  set.seed(seed) # Setting seed for reproduceability
  
  training_data <- vector("list", length = runs)
  test_data <- vector("list", length = runs)
  model_list <- vector("list", length = runs)
  predictions <- vector("list", length = runs)
  
  for(i in 1:runs){
    training_data[[i]] <- binary_reg_dgp2(n = n,
                                          p = p,
                                          p_a = p_a,
                                          beta = beta,
                                          sigma_in = sigma_in,
                                          dirty = dirty,
                                          type = dgp,
                                          v_outlier = 0,
                                          test = FALSE)
    
    test_data[[i]] <- binary_reg_dgp2(n = n,
                                      p = p,
                                      p_a = p_a,
                                      beta = beta,
                                      sigma_in = sigma_in,
                                      dirty = 0,
                                      type = dgp,
                                      v_outlier = 0,
                                      test = TRUE)
  }
  
  
  # Making model formula (for glm function)
  #formula <- paste("y", paste0("X", 1:(p - 1), collapse = " + "), sep = " ~ ") # ? why not:
  formula <- paste("y", paste0("X", 1:p, collapse = " + "), sep = " ~ ") 
  # not including X0 (intercept), does this automatically
  
  # glm (Ordinary logit)
  if(type == "glm"){
    for(r in 1:runs){
      # Fit model
      model_list[[r]] <- glm(formula = formula, 
                             family = binomial(link = "logit"),
                             control = glm.control(maxit = 200),
                             data = training_data[[r]])
      
      # Predict
      predictions[[r]] <- predict(model_list[[r]], 
                                  newdata = test_data[[r]],
                                  type = "response")
    }
  }
  
  # glmrob (robust)
  if(type == "glmrob"){
    for(r in 1:runs){
      # Fit model
      model_list[[r]] <- glmrob(formula = formula, 
                                family = binomial, 
                                method = "WBY",
                                data = training_data[[r]]) 
      
      # Predict
      predictions[[r]] <- 1/(1 + exp(- predict(model_list[[r]], newdata = test_data[[r]])))
      # !!Need to re-transform the log-odds scale because just doing type = "response" gives ERROR :(
    }
  }
  
  # enetLTS (robust, variable selection)
  if(type == "enetLTS"){
    for(r in 1:runs){
      # Prepare data
      X <- as.matrix(training_data[[r]][, - 1, drop = FALSE]) # Drop = false s.t. they don't become vectors
      y <- training_data[[r]][, 1]
      
      X_test <- as.matrix(test_data[[r]][, -1, drop = FALSE])

      # Fit model
      model_list[[r]] <- enetLTS(xx = X,
                                 yy = y,
                                 family = "binomial",
                                 hsize = 0.75,
                                 nfold = 2,
                                 intercept = TRUE)
      # intercept = TRUE set because the model matrix doesn't have intercept atm
      
      # Predict on test sets
      predictions[[r]] <- unname(unlist(predict(model_list[[r]], 
                                                newX = X_test,
                                                type = "response")))
      # Need to unname/unlist the thing otherwise it's a mess
    }
  }
  
  # glmnet
  if(type == "glmnet"){
    for(r in 1:runs){
      # Prepare data
      X <- as.matrix(training_data[[r]][, - 1, drop = FALSE]) # Drop = false s.t. they don't become vectors
      y <- training_data[[r]][, 1]
      
      X_test <- as.matrix(test_data[[r]][, -1, drop = FALSE])
      
      # Fit model
      model_list[[r]] <- cv.glmnet(x = X,
                                   y = y,
                                   family = "binomial",
                                   alpha = 0.7,
                                   type.measure = "deviance")
      
      predictions[[r]] <- predict(model_list[[r]],
                                  newx = X_test,
                                  s = "lambda.1se",
                                  type = "response")
    }
  }
  
  # glmnet2: glmnet with alpha tuning (SLOW!)
  if(type == "glmnet2"){
    # Initializing some values
    lambda_opt <- numeric(runs)
    alpha_opt <- numeric(runs)
    
    # Defining trainControl object (requires caret!)
    trControl_glmnet <- trainControl(method = "repeatedcv",
                                     number = 5,
                                     repeats = 5,
                                     search = "grid",
                                     selectionFunction = "best")
    
    for(r in 1:runs){
      # Preparing data
      X <- as.matrix(training_data[[r]][, - 1, drop = FALSE]) # Drop = false s.t. they don't become vectors
      y <- training_data[[r]][, 1]
      
      X_test <- as.matrix(test_data[[r]][, -1, drop = FALSE])
      
      # Defining training grid lambda bounds (yes differs on data)
      lambda_range <- unique(unlist(lapply(seq(0, 1, 0.025), FUN = function(x){
        init <- glmnet(x = X,
                       y = y,
                       family = "binomial",
                       nlambda = 100,
                       alpha = X)
        lambda <- c(min(init$lambda),
                    max(init$lambda))
      })))
      
      # Filling actual grid
      tuneGrid_glmnet <- expand.grid(.alpha = seq(0, 1, 0.025),
                                     .lambda = seq(min(lambda_range), max(lambda_range)),
                                     length.out = 100)
      
      # Fitting model with created settings
      glmnet_train <- train(x = X,
                            y = factor(y),
                            method = "glmnet",
                            family = "binomial",
                            tuneGrid = tuneGrid_glmnet,
                            trControl = trControl_glmnet)
      
      # Saving optimal tuning params
      alpha_opt[r] <- glmnet_fit$bestTune$alpha
      lambda_opt[r] <- glmnet_fit$bestTune$lambda
      
      # Fitting with optimal settings
      model_list[[r]] <- glmnet(x = x,
                                y = y,
                                family = "binomial",
                                alpha = alpha_opt[r],
                                lambda_opt = lambda_opt[r])
      
      # Predictions
      predictions[[r]] <- predict(model_list[[r]], 
                                  newx = X, 
                                  type = "response")
    }
  }
    
  ## Calculating loss (Same for all models)
  # Negative Log-likelihood (Cross entropy) loss 
  loss <- vector("list", length = runs) # Initializing
  run_avg_loss <- numeric(runs) # Initializing
  
  for(r in 1:runs){ # Over all runs
    y_true <- test_data[[r]][, 1] # Extracting test outcomes
    preds <- predictions[[r]] # Extracting predictions (probabilities!)
    loss[[r]] <- logit_loss(y_true = y_true, prob_hat = preds) # Applying logit_loss function per run
    run_avg_loss[r] <- unlist(loss[[r]]$avg_loss) # Calculating avg per run
  }
  # Calculating avg (over-runs-average) -loglik loss 
  avg_loss <- mean(run_avg_loss) # Calculating avg over all runs
  
  # Classification error loss (same idea as above)
  misclass <- vector("list", lengt = runs)
  run_avg_misclass <- numeric(runs)
  for(r in 1:runs){
    y_true <- test_data[[r]][, 1]
    preds_class <- ifelse(predictions[[r]] > 0.5, yes = 1, no = 0)
    misclass[[r]] <- abs(y_true - preds_class)
    run_avg_misclass[r] <- mean(misclass[[r]])
  }
  
  # Calcualting avg (over-runs-average) missclass loss
  avg_misclass <- mean(run_avg_misclass)
  
  
  # OUPUT
  return(list(models = model_list, 
              train = training_data, 
              test = test_data, 
              preds = predictions, 
              loss = loss, 
              avg_loss = avg_loss, 
              misclass = misclass,
              avg_misclass = avg_misclass))
}
