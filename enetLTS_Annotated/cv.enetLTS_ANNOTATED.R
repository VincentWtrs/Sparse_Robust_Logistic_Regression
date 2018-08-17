##### THIS FILE CONTAINS AN ANNOTATED COPY OF enetLTS::cv.enetLTS() (version 0.1.0) #####

cv.enetLTS_ANNOTATED <- function(index = NULL, xx, yy, family, h, alphas, lambdas, nfold, 
                                 repl, ncores, plot = TRUE) 
{
  
  ### GOAL: Finding the optimal alpha-lambda combination using cross validation and ...
  # giving a plot of the performance for the alpha-lambda combinations. Outside of...
  # calc_evalCrit() and the plotting part at the end it basically only just finds the minimum
  
  ## INPUTS: 
  # index: Array (h x length(lambdas) x length(alphas)).
  # xx: Robustly normalized predictors.
  # yy: Robustly centered outcome (for binomial response: no centering!)
  # family: Response distribution, "binomial" or "gaussian".
  # h: Sample size (e.g. 75) # Note: hsize is a fractional: e.g. n = 100, hsize = 0.75 then h = 75.
  # alphas: A whole sequence (can be length =  1) of alphas, as constructed in enetLTS().
  # lambdas: A whole sequence (can be length = 1) of lambdas, as constructed in enetLTS().
  # nfold: Amount of folds for the cross-validation process.
  # repl: Number specifying amount of replications, how many times to repeat the cross-validation process.
  # ncores: Amount of cores, for parallel processing.
  # plot: Logical denoting if the visual plot for performance should be given.
  
  ## CALLED BY: enetLTS::enetLTS(), a single time for each time for 1 enetLTS() call.
  
  ## OUTPUTS:
  # evalCrit: Named matrix with values of all evaluation criteria for all lambda-alpha combinations
  # minevalCrit: Minimum of the evaluation criterion
  # indexbest: h Observation indices associated with minevalCrit
  # lambdaopt: Lambda associated with minevalCrit
  # alphaopt: Alpha associated with minevalCrit

  ## OUTPUT USED BY: enetLTS::enetLTS()
  
  ## INNER WORKINGS: for each combination of alpha and lambda, enetLTS::calc_evalCrit() is called (lapply loop), ...
  # and averaged over the replications, the optimal lambda and alpha are chosen based on the value of the evaluation criterion, ...
  # such as the mean squared error or the robustified deviance (Bianco-Yohai) in the local variable "evalCrit". If plot == TRUE...
  # a tile plot is given giving a visual representation of the performance.
  
  ## NOTES:
  # 1. Indeed, the function is called cv.enetLTS() but the inner logic does not actually perform cross validation, this is done by the...
  #    internally defined function calc_evalCrit(), which e.g. makes the folds and loop over them. I SPLIT OFF THIS FUNCTION! 
  #    cv.handler.enetLTS()...
  #    would be a more suitable name, with the calc_evalCrit() split-off. 
  # 2. The plotting part should ideally be split-off as well.
  
  # Setting all loss functions to NULL (only come into play for plotting, see below!)
  MNLL <- TMNLL <- RTMSPE <- RMSPE <- NULL
  # MNLL = Mean Negative Loglikelihood
  # TMNLL = Trimmed Mean Negative Loglikelihood
  # RTMSPE = Trimmed Mean Squared Prediction Error
  # RMSPE = Root Mean Squared Prediction Error
  
  # Extracting sample information
  n <- nrow(xx) # Sample size
  p <- ncol(xx) # Dimensionality including or excluding intercept?
  
  # Checking for impossible alpha values:
  wh <- (alphas < 0 | alphas > 1) # Returns TRUEs/FALSEs
  if (sum(wh) > 0) # If > 0, there are impossible values
    stop("alphas can take the values only between 0 and 1")
  # Some other scenarios for the tuning params:
  if (missing(alphas)) 
    stop("provide an alphas sequence")
  if (missing(lambdas)) 
    stop("provide an lambdas sequence")
  
  ## Initiating evalCrit-object matrix over lambdas/alphas
  # Filling with NA; row = lambda values / col = alpha values
  evalCrit <- matrix(NA, nrow = length(lambdas), ncol = length(alphas))
  dimnames(evalCrit) <- list(paste("lambdas", lambdas), paste("alpha", alphas)) # dimnames (colnames/rownames)

  # Making grid which has POSITIONS of the values i.e. 1-1, 2-1, 3-1, ..., 1-2, 2-2, ...
  combis_ind <- expand.grid(1:length(lambdas), 1:length(alphas)) # Somekind of table of indices!
  # So the object combis_ind is only dependent on the amount of alpha and lambda values not on the actual values itself
  
  # Creating sequence of 1, 2, ..., -> AMOUNT OF COMBINATIONS
  indcombi <- 1:nrow(combis_ind)
  
  # calc_evalCrit() function (This was defined here but I SPLIT IT OFF!)
  
  # Multicore applying calc_evalCrit() over ROW NUMBER of combis_ind
  temp_result <- parallel:::mclapply(1:nrow(combis_ind), FUN = calc_evalCrit, # calc_evalCrit now defined in separate .R file
                                     combis_ind = combis_ind, 
                                     alphas = alphas, 
                                     lambdas = lambdas, 
                                     index = index, 
                                     xx = xx, 
                                     yy = yy, 
                                     nfold = nfold, 
                                     repl = repl, 
                                     mc.cores = ncores, 
                                     mc.allow.recursive = FALSE)
  # So at the end of mclapply we have passed all alpha-lambda combinations but at each individual loop of the mclapply only one single
  # row (i.e. one single alpha-lambda) combination is passed through! However temp_Reults
  
  ## Restructuring results
  # Temporary matrix
  temp_result2 <- matrix(unlist(temp_result), 
                         ncol = repl + 2, # +2 Because in col 1 and 2 we save the lambda and alpha INDICES
                         byrow = TRUE)
  # Contains evalcrits: nrows = unique alpha-lambda combinations, ncols = amount of replications + 2 for INDICES of alpha-lambda
  
  # Simplifying: Mean (over repl)
  for (k in 1:nrow(temp_result2)) {
    i <- temp_result2[k, 1] # Getting lambda
    j <- temp_result2[k, 2] # Getting alpha
    
    evalCrit[i, j] <- mean(temp_result2[k, 3:(repl + 2)]) # Getting the means over the replications (first 2 cols are indices)
  }
  
  ## EXTRACTING OPTIMAL SOLUTION (wrt. tuning params)
  optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), arr.ind = TRUE)[1, ]
  minevalCrit <- evalCrit[optind[1], optind[2]]
  indexbest <- index[, optind[1], optind[2]]
  alphas <- round(alphas, 4)
  alpha <- alphas[optind[2]]
  lambdas <- round(lambdas, 4)
  lambda <- lambdas[optind[1]]
  
  
  ############################################### PLOTTING PART ###############################################################
  # Would be more suitable to have this in a separate function...
  
  ### PLOTTING PERFORMANCE GRID
  if (plot == TRUE) {
    print(paste("optimal model: lambda =", lambda, "alpha =", alpha))
    
    # Defining colour palette to be used for the tile fillings:
    lenCol <- length(alphas) * length(lambdas) # Amount of colours
    mycol.b <- colorRampPalette(c("black", "blue2", "purple", "orange", "yellow"))(lenCol) # Interpolates between this colours
    # from help file: These functions return functions that interpolate a set of given colors (...)
    
    # Creating object to be plotted (to a ggplot compliant format)
    rownames(ggmspe) <- lambdas
    colnames(ggmspe) <- alphas
    ggmspe <- melt(ggmspe)
    
    ## Plotting if no index (i.e. only 1 lambda/alpha combo????)
    if (is.null(index)) {
      # Plotting for Binomial
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "TMNLL")
        
        # ggplot
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), y = as.factor(alpha), fill = TMNLL)) + 
          geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + 
          theme(axis.text.x = element_text(angle = -90)) +
          ggtitle(paste0("TMNLL (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
        
        # added this to above:
        #mspeplot <- mspeplot + ggtitle(paste0("TMNLL (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
      
      # Plotting for Gaussian
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RTMSPE")
        
        # ggplot
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), y = as.factor(alpha), fill = RTMSPE)) + 
          geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + 
          theme(axis.text.x = element_text(angle = -90)) +
          ggtitle(paste0("RTMSPE (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")")) # added this immediately here
        
        # added this line to above:
        #mspeplot <- mspeplot + ggtitle(paste0("RTMSPE (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
    }
    
    ## Plotting if yes index (multiple lambda/alpha)
    else {
      # Plotting for Binomial
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "MNLL")
        
        # ggplot
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), y = as.factor(alpha), fill = MNLL)) + 
          geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + 
          theme(axis.text.x = element_text(angle = -90)) +
          ggtitle(paste0("MNLL (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
        
          #mspeplot <- mspeplot + ggtitle(paste0("MNLL (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
      # Plotting for Gaussian
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RMSPE")
        
        # ggplot
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), y = as.factor(alpha), fill = RMSPE)) + 
          geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + 
          theme(axis.text.x = element_text(angle = -90)) +
          ggtitle(paste0("RMSPE (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
        
        #mspeplot <- mspeplot + ggtitle(paste0("RMSPE (minimum at lambda=", lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
    }
    
    mspeplot <- mspeplot + 
      xlab("lambda") + 
      ylab("alpha")
    
    grid.newpage()
    pushViewport(viewport(layout = grid:::grid.layout(1, 1)))
    print(mspeplot, vp = grid:::viewport(layout.pos.row = 1, layout.pos.col = 1))
  }
  
  ############################################### END: PLOTTING PART ###############################################################
  
  # OUTPUT
  return(list(evalCrit = evalCrit, 
              minevalCrit = minevalCrit, 
              indexbest = indexbest, 
              lambdaopt = lambda, 
              alphaopt = alpha))
}