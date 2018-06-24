enetLTS <- function (xx, yy, family = c("gaussian", "binomial"), alphas, 
                     lambdas, lambdaw, hsize = 0.75, intercept = TRUE, nsamp = 500, 
                     s1 = 10, nCsteps = 20, nfold = 5, seed = NULL, plot = TRUE, 
                     repl = 5, para = TRUE, ncores = 1, del = 0.0125, tol = -1e+06, 
                     scal = TRUE, type = c("response", "class")) 
{
  matchedCall <- match.call() # match.call great if there are a lot of optional arguments
  matchedCall[[1]] <- as.name("enetLTS")
  
  # Assigning the arguments to internal variables (?)
  family <- match.arg(family)
  type <- match.arg(type)
  xx <- addColnames(as.matrix(xx))  #?
  nc = dim(yy) #?
  if (is.null(nc)) {
    yy <- as.matrix(yy)
  }
  n <- nrow(xx)
  p <- ncol(xx)
  h <- floor((n + 1) * hsize)
  
  # Some error throwing
  if (repl <= 0) 
    stop("repl has to be a positive number")
  if (nCsteps <= 0) 
    stop("nCsteps has to be a positive number")
  if (type == "class" & family == "gaussian") 
    stop("class type is not available for gaussian family")
  ncores <- rep(ncores, length.out = 1)
  if (is.na(ncores)) 
    ncores <- detectCores()
  if (!is.numeric(ncores) || is.infinite(ncores) || ncores < 
      1) {
    ncores <- 1
    warning("invalid value of 'ncores'; using default value")
  }
  else ncores <- as.integer(ncores)
  
  # Binomial-specific preparation:
  if (family == "binomial") {
    rnames <- rownames(yy)
    rownames(yy) <- 1:nrow(yy)
    y_actual = as.factor(yy)
    ntab = table(y_actual)
    minclass = min(ntab) # why minclass?
    classnames = names(ntab)
  }
  
  # Missing alpha sequence: 41 evenly spaced values
  if (missing(alphas)) {
    alphas <- seq(0, 1, length = 41)
  }
  
  ## User specified: check:
  # Sorting to get nice order
  alphas <- sort(alphas)
  
  # Checking if values outside [0, 1] are given
  wh <- (alphas < 0 | alphas > 1)
  if (sum(wh) > 0)  # Error throwing
    stop("alphas can take the values only between 0 and 1")
  
  ## Missing lambda sequence
  # Gaussian
  if (missing(lambdas) & family == "gaussian") {
    # Bivariate Winsorization method (lambda0()):
    l0 <- robustHD::lambda0(xx, yy, normalize = scal, intercept = intercept)
    lambdas <- seq(l0, 0, by = -0.025 * l0)
  }
  # Binomial
  else if (missing(lambdas) & family == "binomial") {
    # Upper limit for lambda penalty (point-biserial)
    l00 <- lambda00(xx, yy, normalize = scal, intercept = intercept)
    lambdas <- seq(l00, 0, by = -0.025 * l00)
  }
  
  # Data Preparation
  sc <- prepara(xx, yy, family, robu = 1) # enetLTS:::prepara
  x <- sc$xnor
  y <- sc$ycen
  
  # Calling enetLTS:::warmCsteps
  WarmCstepresults <- warmCsteps(x, y, h, n, p, family, alphas, 
                                 lambdas, hsize, nsamp, s1, nCsteps, nfold, para, ncores, 
                                 tol, scal, seed)
  #?
  indexall <- WarmCstepresults$indexall #?
  
  # Results for single alpha, lambda combination
  if ((length(alphas) == 1) & (length(lambdas) == 1)) {
    if (plot == TRUE) 
      warning("There is no meaning to see plot for a single combination of lambda and alpha")
    indexbest <- drop(indexall)
    alphabest <- alphas
    lambdabest <- lambdas
  }
  
  # Results for multiple alpha, lambda combinations
  else {
    CVresults <- cv.enetLTS(indexall, x, y, family, h, alphas, 
                            lambdas, nfold, repl, ncores, plot)
    indexbest <- CVresults$indexbest
    alphabest <- CVresults$alphaopt
    lambdabest <- CVresults$lambdaopt
    evalCritCV <- CVresults$evalCrit
  }
  
  # If scaling == TRUE:
  if (scal) {
    scl <- prepara(xx, yy, family, indexbest, robu = 0)
    xs <- scl$xnor
    ys <- scl$ycen
    fit <- glmnet(xs[indexbest, ], ys[indexbest, ], family, 
                  alpha = alphabest, lambda = lambdabest, standardize = FALSE, 
                  intercept = FALSE)
    
    # Binomial 
    if (family == "binomial") {
      # If no intercept:
      a00 <- if (intercept == F) 
        0
      # If yes intercept:
      else drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% (scl$mux/scl$sigx)) # Something with results from prepara()
      
      raw.coefficients <- drop(as.matrix(fit$beta)/scl$sigx)
      raw.residuals <- -(ys * xs %*% as.matrix(fit$beta)) + log(1 + exp(xs %*% as.matrix(fit$beta)))
      raw.wt <- weight.binomial(xx, yy, c(a00, raw.coefficients), 
                                intercept = intercept, del)
      sclw <- prepara(xx, yy, family, which(raw.wt == 1), robu = 0)
      xss <- sclw$xnor
      yss <- sclw$ycen
      
      # Missing lambda
      if (missing(lambdaw)) {
        lambdaw <- cv.glmnet(xss[which(raw.wt == 1), 
                                 ], yss[which(raw.wt == 1)], family = family, 
                             nfolds = 5, alpha = alphabest, standardize = FALSE, 
                             intercept = FALSE, type.measure = "mse")$lambda.min
      }
      # 1 User-supplied lambda
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      # Multiple user-supplied lambda
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        lambdaw <- cv.glmnet(xss[which(raw.wt == 1), 
                                 ], yss[which(raw.wt == 1)], family = family, 
                             lambda = lambdaw, nfolds = 5, alpha = alphabest, 
                             standardize = FALSE, intercept = FALSE, type.measure = "mse")$lambda.min
      }
      fitw <- glmnet(xss[which(raw.wt == 1), ], yss[which(raw.wt == 1)], family, alpha = alphabest, lambda = lambdaw, 
                     standardize = FALSE, intercept = FALSE)
      a0 <- if (intercept == F) 
        0
      else drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% 
                  (sclw$mux/sclw$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      wgt <- weight.binomial(xx, yy, c(a0, coefficients), 
                             intercept, del)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, 
                                                        coefficients)) + log(1 + exp(cbind(1, xx) %*% 
                                                                                       c(a0, coefficients)))
    }
    else if (family == "gaussian") {
      a00 <- if (intercept == F) 
        0
      else drop(scl$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% 
                  (scl$mux/scl$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/scl$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(raw.residuals, indexbest, 
                                del)$we
      sclw <- prepara(xx, yy, family, which(raw.wt == 1), 
                      robu = 0)
      xss <- sclw$xnor
      yss <- sclw$ycen
      if ((missing(lambdaw))) {
        lambdaw <- cv.glmnet(xss[which(raw.wt == 1), 
                                 ], yss[which(raw.wt == 1)], family = family, 
                             nfolds = 5, alpha = alphabest, standardize = FALSE, 
                             intercept = FALSE, type.measure = "mse")$lambda.min
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        lambdaw <- cv.glmnet(xss[which(raw.wt == 1), 
                                 ], yss[which(raw.wt == 1)], family = family, 
                             lambda = lambdaw, nfolds = 5, alpha = alphabest, 
                             standardize = FALSE, intercept = FALSE, type.measure = "mse")$lambda.min
      }
      fitw <- glmnet(xss[which(raw.wt == 1), ], yss[which(raw.wt == 
                                                            1)], family, alpha = alphabest, lambda = lambdaw, 
                     standardize = FALSE, intercept = FALSE)
      a0 <- if (intercept == F) 
        0
      else drop(sclw$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% 
                  (sclw$mux/sclw$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sclw$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, 
                                                      coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(reweighted.residuals, raw.wt == 
                               1, del)$we
    }
  }
  else {
    fit <- glmnet(x[indexbest, ], y[indexbest, ], family, 
                  alpha = alphabest, lambda = lambdabest, standardize = FALSE, 
                  intercept = FALSE)
    if (family == "binomial") {
      a00 <- if (intercept == F) 
        0
      else drop(fit$a0 - as.vector(as.matrix(fit$beta)) %*% 
                  (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- -(y * x %*% as.matrix(fit$beta)) + 
        log(1 + exp(x %*% as.matrix(fit$beta)))
      raw.wt <- weight.binomial(xx, yy, c(a00, raw.coefficients), 
                                intercept, del)
      if (missing(lambdaw)) {
        lambdaw <- cv.glmnet(x[which(raw.wt == 1), ], 
                             y[which(raw.wt == 1)], family = family, nfolds = 5, 
                             alpha = alphabest, standardize = FALSE, intercept = FALSE, 
                             type.measure = "mse")$lambda.min
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        lambdaw <- cv.glmnet(x[which(raw.wt == 1), ], 
                             y[which(raw.wt == 1)], family = family, lambda = lambdaw, 
                             nfolds = 5, alpha = alphabest, standardize = FALSE, 
                             intercept = FALSE, type.measure = "mse")$lambda.min
      }
      fitw <- glmnet(x[which(raw.wt == 1), ], y[which(raw.wt == 
                                                        1)], family, alpha = alphabest, lambda = lambdaw, 
                     standardize = FALSE, intercept = FALSE)
      a0 <- if (intercept == F) 
        0
      else drop(fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% 
                  (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      wgt <- weight.binomial(xx, yy, c(a0, coefficients), 
                             intercept, del)
      reweighted.residuals <- -(yy * cbind(1, xx) %*% c(a0, 
                                                        coefficients)) + log(1 + exp(cbind(1, xx) %*% 
                                                                                       c(a0, coefficients)))
    }
    else if (family == "gaussian") {
      a00 <- if (intercept == F) 
        0
      else drop(sc$muy + fit$a0 - as.vector(as.matrix(fit$beta)) %*% 
                  (sc$mux/sc$sigx))
      raw.coefficients <- drop(as.matrix(fit$beta)/sc$sigx)
      raw.residuals <- yy - cbind(1, xx) %*% c(a00, raw.coefficients)
      raw.rmse <- sqrt(mean(raw.residuals^2))
      raw.wt <- weight.gaussian(raw.residuals, indexbest, 
                                del)$we
      if (missing(lambdaw)) {
        lambdaw <- cv.glmnet(x[which(raw.wt == 1), ], 
                             y[which(raw.wt == 1)], family = family, nfolds = 5, 
                             alpha = alphabest, standardize = FALSE, intercept = FALSE, 
                             type.measure = "mse")$lambda.min
      }
      else if (!missing(lambdaw) & length(lambdaw) == 1) {
        lambdaw <- lambdaw
      }
      else if (!missing(lambdaw) & length(lambdaw) > 1) {
        lambdaw <- cv.glmnet(x[which(raw.wt == 1), ], 
                             y[which(raw.wt == 1)], family = family, lambda = lambdaw, 
                             nfolds = 5, alpha = alphabest, standardize = FALSE, 
                             intercept = FALSE, type.measure = "mse")$lambda.min
      }
      fitw <- glmnet(x[which(raw.wt == 1), ], y[which(raw.wt == 
                                                        1)], family, alpha = alphabest, lambda = lambdaw, 
                     standardize = FALSE, intercept = FALSE)
      a0 <- if (intercept == F) 
        0
      else drop(sc$muy + fitw$a0 - as.vector(as.matrix(fitw$beta)) %*% 
                  (sc$mux/sc$sigx))
      coefficients <- drop(as.matrix(fitw$beta)/sc$sigx)
      reweighted.residuals <- yy - cbind(1, xx) %*% c(a0, 
                                                      coefficients)
      reweighted.rmse <- sqrt(mean(reweighted.residuals^2))
      wgt <- weight.gaussian(reweighted.residuals, raw.wt == 
                               1, del)$we
    }
  }
  num.nonzerocoef <- sum(coefficients != 0)
  intercept <- isTRUE(intercept)
  if (intercept) 
    xx <- addIntercept(xx)
  if (intercept) {
    coefficients <- c(a0, coefficients)
    raw.coefficients <- c(a00, raw.coefficients)
  }
  else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  if (family == "binomial") {
    u <- xx %*% raw.coefficients
    raw.fitted.values <- if (type == "class") {
      ifelse(u <= 0.5, 0, 1)
    }
    else if (type == "response") {
      1/(1 + exp(-u))
    }
    uu <- xx %*% coefficients
    fitted.values <- if (type == "class") {
      ifelse(uu <= 0.5, 0, 1)
    }
    else if (type == "response") {
      1/(1 + exp(-uu))
    }
  }
  else if (family == "gaussian") {
    raw.fitted.values <- xx %*% raw.coefficients
    fitted.values <- xx %*% coefficients
  }
  if (family == "binomial") {
    objective <- h * (mean((-yy[indexbest] * (xx[indexbest, 
                                                 ] %*% coefficients)) + log(1 + exp(xx[indexbest, 
                                                                                       ] %*% coefficients))) + lambdabest * sum(1/2 * (1 - 
                                                                                                                                         alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  }
  else if (family == "gaussian") {
    objective <- h * ((1/2) * mean((yy[indexbest] - xx[indexbest, 
                                                       ] %*% coefficients)^2) + lambdabest * sum(1/2 * (1 - 
                                                                                                          alphabest) * coefficients^2 + alphabest * abs(coefficients)))
  }
  if (intercept) {
    coefficients <- coefficients[-1]
    raw.coefficients <- raw.coefficients[-1]
  }
  else {
    coefficients <- coefficients
    raw.coefficients <- raw.coefficients
  }
  inputs <- list(xx = xx, yy = yy, family = family, alphas = alphas, 
                 lambdas = lambdas, lambdaw = lambdaw, hsize = hsize, 
                 h = h, nsamp = nsamp, s1 = s1, nCsteps = nCsteps, nfold = nfold, 
                 intercept = intercept, repl = repl, para = para, ncores = ncores, 
                 del = del, scal = scal)
  
  # Binomial
  if (family == "binomial") {
    output <- list(objective = objective, best = sort(indexbest), 
                   raw.wt = raw.wt, wt = wgt, a00 = a00, raw.coefficients = raw.coefficients, 
                   a0 = a0, coefficients = coefficients, alpha = alphabest, 
                   lambda = lambdabest, lambdaw = lambdaw, num.nonzerocoef = num.nonzerocoef, 
                   h = h, raw.residuals = drop(raw.residuals), residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), raw.fitted.values = drop(raw.fitted.values), 
                   classnames = classnames, classsize = ntab, inputs = inputs, 
                   indexall = indexall, call = sys.call())
  }
  # Gaussian
  else if (family == "gaussian") {
    output <- list(objective = objective, best = sort(indexbest), 
                   raw.wt = raw.wt, wt = wgt, a00 = a00, raw.coefficients = raw.coefficients, 
                   a0 = a0, coefficients = coefficients, alpha = alphabest, 
                   lambda = lambdabest, lambdaw = lambdaw, num.nonzerocoef = num.nonzerocoef, 
                   h = h, raw.residuals = drop(raw.residuals), residuals = drop(reweighted.residuals), 
                   fitted.values = drop(fitted.values), raw.fitted.values = drop(raw.fitted.values), 
                   raw.rmse = raw.rmse, rmse = reweighted.rmse, inputs = inputs, 
                   indexall = indexall, call = sys.call())
  }
  class(output) <- "enetLTS"
  output$call <- matchedCall
  output
}