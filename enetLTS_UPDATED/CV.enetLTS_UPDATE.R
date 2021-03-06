cv.enetLTS_UPDATE <- function (index = NULL, xx, yy, family, h, alphas, lambdas, nfold, 
          repl, ncores, plot = TRUE, ic_type = NULL) 
{
  MNLL <- TMNLL <- RTMSPE <- RMSPE <- NULL
  n <- nrow(xx)
  p <- ncol(xx)
  wh <- (alphas < 0 | alphas > 1)
  if (sum(wh) > 0) 
    stop("alphas can take the values only between 0 and 1")
  if (missing(alphas)) 
    stop("provide an alphas sequence")
  if (missing(lambdas)) 
    stop("provide an lambdas sequence")
  evalCrit <- matrix(NA, nrow = length(lambdas), ncol = length(alphas))
  dimnames(evalCrit) <- list(paste("lambdas", lambdas), paste("alpha", 
                                                              alphas))
  combis_ind <- expand.grid(1:length(lambdas), 1:length(alphas))
  indcombi <- 1:nrow(combis_ind)
  
  ## REMOVED: DEFINITION OF calc_evalCrit(), THIS HAS BEEN SPLIT OFF TO SEPARATE FILE
  
  # Running for all combinations of alpha, lambda
  ### NEW: UPDATED TO CALL calc_evalcrit_UPDATE instead of the original one ALSO: added argument family = family and ic_type to be passed on as well
  temp_result <- parallel:::mclapply(1:nrow(combis_ind), FUN = calc_evalCrit_UPDATE, 
                          combis_ind = combis_ind, 
                          alphas = alphas, 
                          lambdas = lambdas, 
                          index = index, 
                          xx = xx, 
                          yy = yy, 
                          nfold = nfold, 
                          repl = repl, 
                          family = family, # NEW 
                          ic_type = ic_type, # NEW
                          mc.cores = ncores, 
                          mc.allow.recursive = FALSE)
  # Restructuring output
  temp_result2 <- matrix(unlist(temp_result), 
                         ncol = repl + 2, # + 2
                         byrow = TRUE)
  # # NEW:
  # if(percentile_method == TRUE){
  #   for (k in 1:nrow(temp_result2)) {
  #     i <- temp_result2[k, 1]
  #     j <- temp_result2[k, 2]
  #     evalCrit[i, j] <- quantile(temp_result2[k, 3:(repl + 2)])
  #   }
  # }
  
  for (k in 1:nrow(temp_result2)) {
    i <- temp_result2[k, 1]
    j <- temp_result2[k, 2]
    evalCrit[i, j] <- mean(temp_result2[k, 3:(repl + 2)])
  }
  optind <- which(evalCrit == min(evalCrit, na.rm = TRUE), 
                  arr.ind = TRUE)[1, ]
  minevalCrit <- evalCrit[optind[1], optind[2]]
  indexbest <- index[, optind[1], optind[2]]
  alphas <- round(alphas, 10)
  alpha <- alphas[optind[2]] 
  lambdas <- round(lambdas, 10) # NEW: less rounding
  lambda <- lambdas[optind[1]] # NEW: less rounding
  
  ### PLOTTING (if plot = TRUE)
  if (plot == TRUE) {
    print(paste("optimal model: lambda =", lambda, "alpha =", 
                alpha))
    lenCol <- length(alphas) * length(lambdas)
    mycol.b <- colorRampPalette(c("black", "blue2", "purple", 
                                  "orange", "yellow"))(lenCol)
    ggmspe <- evalCrit
    rownames(ggmspe) <- lambdas
    colnames(ggmspe) <- alphas
    ggmspe <- reshape:::melt(ggmspe) # NEW: Trying to force getting it from here
    if (is.null(index)) {
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "TMNLL")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = TMNLL)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("TMNLL (minimum at lambda=", 
                                              lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RTMSPE")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = RTMSPE)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("RTMSPE (minimum at lambda=", 
                                              lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
    }
    else {
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "MNLL")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = MNLL)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("MNLL (minimum at lambda=", 
                                              lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RMSPE")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = RMSPE)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("RMSPE (minimum at lambda=", 
                                              lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
    }
    mspeplot <- mspeplot + xlab("lambda") + ylab("alpha")
    grid:::grid.newpage() # NEW: added grid::: to force it getting the right package when changing the function in namespace
    grid:::pushViewport(grid:::viewport(layout = grid:::grid.layout(1, 1)))
    print(mspeplot, vp = grid:::viewport(layout.pos.row = 1, layout.pos.col = 1))
  }
  return(list(evalCrit = evalCrit, minevalCrit = minevalCrit, 
              indexbest = indexbest, lambdaopt = lambda, alphaopt = alpha))
}