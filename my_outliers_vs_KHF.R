### Type of outliers example ###

# EXPLANTION: in their paper Kurnatz et al. (2017) propose the following outlier sytem, I think it does not make sense...
#  hence I show an univariate simple example (p = 1) where their outliers do not really cause any trouble for the current 
#  methods but mine do. First I generate some data according to my idea of outliers and their idea of outliers using the ...
#  binary_dgp() function, then I run the ordinary MLE using glm() to show that KHF-type outliers are no problem at all, 
#  while mine are, then using glmrob(, method = "WBY") which is weighted Bianco-Yohai robust estimation, I run the thing again, ...
#  this robust method can indeed handle both situations, as a check. Some plotting was also done



## Generating data
# My idea of vertical outliers  
train_VW <- binary_reg_dgp2(n = 100, 
                          beta = 1, 
                          beta0 = 1, 
                          dirty = 0.1, 
                          v_outlier = "VW") 
# KHF (2017) idea of outliers
train_KHF <- binary_reg_dgp2(n = 100, 
                          beta = 1, 
                          beta0 = 1,
                          dirty = 0.1,
                          v_outlier = "KHF")

## Fitting logistic regression models with MLE (ordinary)
# My of outliers as data
logit_VW <- glm(y ~ X, data = train_VW, family = "binomial")

# KHF idea of outliers as data
logit_KHF <- glm(y ~ X, data = train_KHF, family = "binomial")

# Checking coefficients
coef(logit_VW) # Sign has turned versus true
coef(logit_KHF) # Same sign

# Checking loss (Negative loglik)
-logLik(logit_VW) # 60.56415
-logLik(logit_KHF) # 37.76741 (BEST, minimizing risk/loss) -> Hence probably not real outliers

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
     main = "My (VW) outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_VW_preds)

plot(train_KHF$X, train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_KHF_preds)
par(mfrow = c(1, 1))

# Plotting (Zoomed-in)
xlim <- c(-25, 25)
par(mfrow = c(1, 2))
plot(train_VW$X, train_VW$y,
     xlim = xlim,
     main = "My (VW) outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_VW_preds)

plot(train_KHF$X, train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit_KHF_preds)
par(mfrow = c(1, 1))


## Trying glmRob
logit_rob_VW <- glmrob(y ~ X, data = train_VW, family = "binomial", method = "WBY")
logit_rob_KHF <- glmrob(y ~ X, data = train_KHF, family = "binomial", method = "WBY")
# Ignore warnings

coef(logit_rob_VW)
coef(logit_rob_KHF)

# Predictions
logit_rob_VW_preds <- predict(logit_rob_VW, newdata = data.frame(X = x_grid), type = "response") # These don't want to work for some reason:
logit_rob_KHF_preds <- predict(logit_rob_KHF, newdata = data.frame(X = x_grid), type = "response") # These don't want to work for some reason:

# MANUALLY getting predictions
logit_rob_VW_preds <- exp(coef(logit_rob_VW)[1] + coef(logit_rob_VW)[2] * x_grid)/(1 + exp(coef(logit_rob_VW)[1] + coef(logit_rob_VW)[2] * x_grid))
plot(x = x_grid, 
     y = logit_rob_VW_preds,
     type = "l")

logit_rob_KHF_preds <- exp(coef(logit_rob_KHF)[1] + coef(logit_rob_KHF)[2] * x_grid)/(1 + exp(coef(logit_rob_KHF)[1] + coef(logit_rob_KHF)[2] * x_grid))
plot(x = x_grid, 
     y = logit_rob_KHF_preds,
     type = "l")

# Plotting (Robust, Zoomed-out)
xlim <- c(-80, 80)
par(mfrow = c(1, 2))
plot(x = train_VW$X, 
     y = train_VW$y,
     xlim = xlim,
     main = "My (VW) outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_VW_preds)

plot(x = train_KHF$X, 
     y = train_KHF$y,
     xlim = xlim,
     main = "KHF Outliers (ROBUST)",
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
     main = "My (VW) outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_VW_preds)

plot(x = train_KHF$X, 
     y = train_KHF$y, 
     xlim = xlim,
     main = "KHF Outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob_KHF_preds)
par(mfrow = c(1, 1))
