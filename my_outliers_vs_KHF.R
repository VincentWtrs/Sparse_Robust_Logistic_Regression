# Generating data
train1 <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0.1, v_outlier = 0) # My idea of vertical outliers
train2 <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0.1, v_outlier = 1) # KHF idea of vertical outliers

## Fitting logits (ordinary)
logit1 <- glm(y ~ X, data = train1, family = "binomial")
logit2 <- glm(y ~ X, data = train2, family = "binomial")

coef(logit1)
coef(logit2)

# Checking loss
-logLik(logit1)/nrow(train1) # 0.602
-logLik(logit2)/nrow(train2) # 0.364 (BEST, minimizing risk/loss) -> Hence probably not real outliers

## Plotting logistic fit
# Generating predicted probabilities along a fine grid
x_grid <- seq(-100, 100, by = 0.01)

# Getting predictions
logit1_preds <- predict(logit1, 
                        newdata = data.frame(X = x_grid),
                        type = "response")
logit2_preds <- predict(logit2,
                        newdata = data.frame(X = x_grid),
                        type = "response")

# Plotting (Zoomed-out)
xlim <- c(-80, 80)
par(mfrow = c(1, 2))
plot(train1$X, train1$y, 
     xlim = xlim,
     main = "My (VW) outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit1_preds)

plot(train2$X, train2$y, 
     xlim = xlim,
     main = "KHF Outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit2_preds)
par(mfrow = c(1, 1))

# Plotting (Zoomed-in)
xlim <- c(-25, 25)
par(mfrow = c(1, 2))
plot(train1$X, train1$y,
     xlim = xlim,
     main = "My (VW) outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit1_preds)

plot(train2$X, train2$y, 
     xlim = xlim,
     main = "KHF Outliers (NONrobust fit)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, y = logit2_preds)
par(mfrow = c(1, 1))


## Trying glmRob
logit_rob1 <- glmrob(y ~ X, data = train1, family = "binomial", method = "WBY")
logit_rob2 <- glmrob(y ~ X, data = train2, family = "binomial", method = "WBY")

coef(logit_rob1)
coef(logit_rob2)

# Predictions
# These don't want to work for some reason:
logit_rob1_preds <- predict(logit_rob1, newdata = data.frame(X = x_grid), type = "response")
logit_rob2_preds <- predict(logit_rob2, newdata = data.frame(X = x_grid), type = "response")

# MANUALLY getting predictions
logit_rob1_preds <- exp(coef(logit_rob1)[1] + coef(logit_rob1)[2] * x_grid)/(1 + exp(coef(logit_rob1)[1] + coef(logit_rob1)[2] * x_grid))
plot(x_grid, logit_rob1_preds)

logit_rob2_preds <- exp(coef(logit_rob2)[1] + coef(logit_rob2)[2] * x_grid)/(1 + exp(coef(logit_rob2)[1] + coef(logit_rob2)[2] * x_grid))
plot(x_grid, logit_rob2_preds)

# Plotting (Robust, Zoomed-out)
xlim <- c(-80, 80)
par(mfrow = c(1, 2))
plot(x = train1$X, 
     y = train1$y,
     xlim = xlim,
     main = "My (VW) outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob1_preds)

plot(x = train2$X, 
     y = train2$y,
     xlim = xlim,
     main = "KHF Outliers (ROBUST)",
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob2_preds)
par(mfrow = c(1, 1))

# Plotting (Robust, Zoomed-in)
xlim <- c(-25, 25)
par(mfrow = c(1, 2))
plot(x = train1$X, 
     y = train1$y, 
     xlim = xlim,
     main = "My (VW) outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob1_preds)

plot(x = train2$X, 
     y = train2$y, 
     xlim = xlim,
     main = "KHF Outliers (ROBUST)", 
     xlab = "X",
     ylab = "Probability")
lines(x = x_grid, 
      y = logit_rob2_preds)
par(mfrow = c(1, 1))