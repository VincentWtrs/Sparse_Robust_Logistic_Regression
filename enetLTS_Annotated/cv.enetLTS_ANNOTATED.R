
# IN THIS FILE I WILL ANNOTATE THE enetLTS FUNCTION:

# This function handles the cross validation scheme for the enetLTS fitting procedure:
cv.enetLTS_ANNOTATED <- function (index = NULL, xx, yy, family, h, alphas, lambdas, nfold, 
                        repl, ncores, plot = TRUE) 
{
  ## Setting all loss functions to NULL:
  # MNLL = Mean Negative Loglikelihood
  # TMNLL ?
  # RTMSPE = 
  # RMSPE = Root Mean Squared Prediction Error
  MNLL <- TMNLL <- RTMSPE <- RMSPE <- NULL
  
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
  
  ## Initiating evalCrit matrix over lambdas/alphas
  # Filling with NA; row = lambda values / col = alpha values
  evalCrit <- matrix(NA, nrow = length(lambdas), ncol = length(alphas))
  
  # Giving appropriate headers (i.e. each row and column has its OWN name!)
  dimnames(evalCrit) <- list(paste("lambdas", lambdas), paste("alpha", alphas))
  
  # Making grid which has POSITIONS of the values e.g. 1-1, 2-1, 3-1, ... 1-2, 2-2, ...
  combis_ind <- expand.grid(1:length(lambdas), 1:length(alphas)) # Somekind of table of indices!
  
  # Creating sequence of 1, 2, ..., AMOUNT OF COMBINATIONS
  indcombi <- 1:nrow(combis_ind)
  
  # calc_evalCrit() function (I HAVE SPLIT THIS OFF!!!!)
  
  # Multicore applying over the calc_evalCrit
  temp_result <- mclapply(1:nrow(combis_ind), FUN = calc_evalCrit, 
                          combis_ind = combis_ind, alphas = alphas, lambdas = lambdas, 
                          index = index, xx = xx, yy = yy, nfold = nfold, repl = repl, 
                          mc.cores = ncores, mc.allow.recursive = FALSE)
  temp_result2 <- matrix(unlist(temp_result), ncol = repl + 2, byrow = TRUE)
  for (k in 1:nrow(temp_result2)) {
    i <- temp_result2[k, 1]
    j <- temp_result2[k, 2]
    evalCrit[i, j] <- mean(temp_result2[k, 3:(repl + 2)])
  }
  
  optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), 
                  arr.ind = TRUE)[1, ]
  minevalCrit <- evalCrit[optind[1], optind[2]]
  indexbest <- index[, optind[1], optind[2]]
  alphas <- round(alphas, 4)
  alpha <- alphas[optind[2]]
  lambdas <- round(lambdas, 4)
  lambda <- lambdas[optind[1]]
  
  ### PLOTTING PERFORMANCE GRID
  if (plot == TRUE) {
    print(paste("optimal model: lambda =", lambda, "alpha =", 
                alpha))
    
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
    pushViewport(viewport(layout = grid.layout(1, 1)))
    print(mspeplot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  }
  
  # Return values
  return(list(evalCrit = evalCrit, 
              minevalCrit = minevalCrit, 
              indexbest = indexbest, 
              lambdaopt = lambda, 
              alphaopt = alpha))
}