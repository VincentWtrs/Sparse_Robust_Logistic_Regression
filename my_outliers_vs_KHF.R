### Type of outliers example ###

## PART 1: Single fit using Weighted Bianco-Yohai and MLE ##

# EXPLANTION: in their paper Kurnatz et al. (2017) propose the following outlier sytem, I think it does not make sense...
#  hence I show an univariate simple example (p = 1) where their outliers do not really cause any trouble for the current 
#  methods but mine do. First I generate some data according to my idea of outliers and their idea of outliers using the ...
#  binary_dgp() function, then I run the ordinary MLE using glm() to show that KHF-type outliers are no problem at all, 
#  while mine are, then using glmrob(, method = "WBY") which is weighted Bianco-Yohai robust estimation, I run the thing again, ...
#  this robust method can indeed handle both situations, as a check. Some plotting was also done

# Seed
# set.seed(1234)

# CHOOSE TO SET BETA POS OR NEG
pos <- TRUE # set to FALSE for negative beta

# Setting beta on user setting
beta <- ifelse(pos == TRUE,
               yes = 1.5,
               no = -1.5)
# beta = 1.5 works good, see later to see that beta = 1 give some trouble sometimes

## Generating data
# My idea of vertical outliers  
train_VW <- binary_reg_dgp2(n = 100, 
                          beta = beta, 
                          beta0 = 0.2, 
                          dirty = 0.1, 
                          v_outlier = "VW",
                          type = "bernoulli",
                          outlier_mean = 10) 
# KHF (2017) idea of outliers
train_KHF <- binary_reg_dgp2(n = 100, 
                          beta = beta, 
                          beta0 = 0.2,
                          dirty = 0.1,
                          v_outlier = "KHF",
                          type = "bernoulli",
                          outlier_mean = 10)

## Fitting logistic regression models with MLE (ordinary)
# My of outliers as data
logit_VW <- glm(y ~ X, data = train_VW, family = "binomial")

# KHF idea of outliers as data
logit_KHF <- glm(y ~ X, data = train_KHF, family = "binomial")

# Checking coefficients
coef(logit_VW) # Sign has turned versus true
coef(logit_KHF) # Same sign, almost exactly true values

# Checking loss (Negative loglik)
-1 * (logLik(logit_VW)) # 57.8694
-1 * (logLik(logit_KHF)) # 48.2235 (BEST, minimizing risk/loss) -> Hence probably not real outliers

## Plotting logistic fit
# Generating predicted probabilities along a fine grid of x values
x_grid <- seq(-100, 100, by = 0.01)

# Getting predictions (on the grid)
logit_VW_preds <- predict(logit_VW, 
                        newdata = data.frame(X = x_grid),
                        type = "response")
logit_KHF_preds <- predict(logit_KHF,
                        newdata = data.frame(X = x_grid),
                        type = "response")

# Plotting (Zoomed-out)
xlim <- c(-80, 80)
par(mfrow = c(1, 2))
plot(train_VW$X, train_VW$y, 
     xlim = xlim,
     main = "VW outliers (Non-robust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_VW_preds)

plot(train_KHF$X, train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (Non-robust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_KHF_preds)
par(mfrow = c(1, 1))

# Plotting (Zoomed-in)
xlim <- c(-25, 25)
par(mfrow = c(1, 2))
plot(train_VW$X, train_VW$y,
     xlim = xlim,
     main = "VW outliers (Non-robust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_VW_preds)

plot(train_KHF$X, train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (Non-robust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_KHF_preds)
par(mfrow = c(1, 1))


## Fitting Weighted Bianco-Yohai (glmrob, WBY)
logit_rob_VW <- glmrob(y ~ X, data = train_VW, family = "binomial", method = "WBY")
logit_rob_KHF <- glmrob(y ~ X, data = train_KHF, family = "binomial", method = "WBY")
# Ignore warnings

# Checking coefficients
coef(logit_rob_VW) 
coef(logit_rob_KHF)

## Plotting results
# Predictions
logit_rob_VW_preds <- predict(logit_rob_VW, newdata = data.frame(X = x_grid), type = "response") # These don't want to work for some reason:
logit_rob_KHF_preds <- predict(logit_rob_KHF, newdata = data.frame(X = x_grid), type = "response") # These don't want to work for some reason:

# MANUALLY getting predictions
logit_rob_VW_preds <- exp(coef(logit_rob_VW)[1] + coef(logit_rob_VW)[2] * x_grid)/(1 + exp(coef(logit_rob_VW)[1] + coef(logit_rob_VW)[2] * x_grid))
logit_rob_KHF_preds <- exp(coef(logit_rob_KHF)[1] + coef(logit_rob_KHF)[2] * x_grid)/(1 + exp(coef(logit_rob_KHF)[1] + coef(logit_rob_KHF)[2] * x_grid))


# Plotting (Robust, Zoomed-out)
xlim <- c(-80, 80)
par(mfrow = c(1, 2))
plot(x = train_VW$X, 
     y = train_VW$y,
     xlim = xlim,
     main = "VW outliers (WBY fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_VW_preds)

plot(x = train_KHF$X, 
     y = train_KHF$y,
     xlim = xlim,
     main = "KHF Outliers (WBY fit)",
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_KHF_preds)
par(mfrow = c(1, 1))

# Plotting (Robust, Zoomed-in)
xlim <- c(-25, 25)
par(mfrow = c(1, 2))
plot(x = train_VW$X, 
     y = train_VW$y, 
     xlim = xlim,
     main = "VW outliers (WBY fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_VW_preds)

plot(x = train_KHF$X, 
     y = train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (WBY fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_KHF_preds)
par(mfrow = c(1, 1))


## TEMPORARY CONCLUSION: It can be seen that the my (VW) type of outliers cause a lot more trouble to the estimators...
# to the estimators, even the WBY estimator sometimes fails. To get an idea about the frequency of this failutre we.
# Only for the case where beta < 0, does the KHF method seemingly (unknowingly?) give the correct type of problem-causing...
# outliers. The strange beta case for positive and negative values is investigated through a simulation study running ..
# multiple runs for multiple betas

#############################################################################################################################
#############################################################################################################################


## PART 2: Simulation study using multiple values of beta ##


# Setting simulation settings
n <- 100 # Sample size
runs <- 49 # Amount of simulation runs
beta0 <- 0.2 # intercept
beta <- c(seq(from = 0.1, to = 1.8, by = 0.2)) # Beta range
ylim <- c(min(beta) - 0.5, max(beta) + 2) # ylim range for plotting based on the betas

# Creating data for each run (i.e. the sampling variability)
data <- lapply(beta, FUN = function(x) lapply(1:runs, FUN = function(z) binary_reg_dgp2(n = 100, 
                                                                                        beta = x, #!
                                                                                        beta0 = beta0, 
                                                                                        dirty = 0.1, 
                                                                                        v_outlier = "VW",
                                                                                        type = "bernoulli")))
# Double lapply, for each beta, for and multiple runs to get sampling variability


# Fitting logistic by MLE
logit_mle <- lapply(1:length(beta), FUN = function(x) lapply(1:runs, FUN = function(z) glm(y ~ X,
                                                                                           family = "binomial",
                                                                                           data = data[[x]][[z]])))

# Gathering coefficients MLE
coefs_mle <- vector("list", length = length(beta))
for(i in 1:length(beta)){
  coefs_mle[[i]] <- matrix(NA, nrow = runs, ncol = 2)
  colnames(coefs_mle[[i]]) <- c("Intercept", "beta1")
  for(j in 1:runs){
    coefs_mle[[i]][j, ] <- coef(logit_mle[[i]][[j]])
  }
}


## Plotting MLE
# Boxplot
for(i in 1:length(beta)){
  boxplot(coefs_mle[[i]], 
          ylim = ylim,
          main = "Box plot of simulation using MLE")
  points(x = c(1, 2),
         y = c(beta0, beta[i]),
         col = "red",
         pch = 20, # Red filled dot
         cex = 2) # size = 2
}



# Fitting WBY
logit_wby <- lapply(1:length(beta), FUN = function(x) lapply(1:runs, FUN = function(z) glmrob(y ~ X,
                                                                                           family = "binomial",
                                                                                           method = "WBY",
                                                                                           data = data[[x]][[z]])))

# Gathering coefficients WBY
coefs_wby <- vector("list", length = length(beta))
for(i in 1:length(beta)){
  coefs_wby[[i]] <- matrix(NA, nrow = runs, ncol = 2)
  colnames(coefs_wby[[i]]) <- c("Intercept", "beta1")
  for(j in 1:runs){
    coefs_wby[[i]][j, ] <- coef(logit_wby[[i]][[j]])
  }
}


## Plotting WBY
# Boxplot
for(i in 1:length(beta)){
  boxplot(coefs_wby[[i]], 
          ylim = ylim,
          main = "Box plot of simulation using WBY")
  points(x = c(1, 2),
         y = c(beta0, beta[i]),
         col = "red",
         pch = 20, # Red filled dot
         cex = 2) # size = 2
}

## CONCLUSION: The MLE completely fails in all contamination settings. The WBY does better but only for quite...
# large true betas. If below 1.2-ish they seem to get severly skewed towards 0, even though 1.2 is a quite big...
# value, certainly in terms of odds ratio (exp(beta)). 

## SPECULATION: This worrying behaviour of the WBY might have to do with the still M-estimator style? type. ...
# Might also have to do with the properties of the outliers.