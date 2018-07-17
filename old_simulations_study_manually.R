
# Clean
reps <- 100
data <- vector("list", length = reps)
test <- vector("list", length = reps)
models <- vector("list", length = reps)
predictions <- matrix(NA, nrow = reps, ncol = 100)
loss <- vector("numeric", length = reps)
for(i in 1:reps){
  data[[i]] <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0, v_outlier = 0) # My idea of vertical outliers
  test[[i]] <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0, v_outlier = 0)
  data_now <- data[[i]]
  test_now <- test[[i]]
  y_test_now <- test_now[, 1]
  X_test_now <- test_now[, -1, drop = FALSE] # otherwise it makes it a vector :(
  
  models[[i]] <- glm(y ~ X, data = data_now, family = "binomial")
  predictions[i, ] <- predict(models[[i]], newdata = X_test_now, type = "response")
  loss[i] <- logit_loss(y_true = y_test_now, prob_hat = predictions[i, ])
  
}

coefs_clean <- t(sapply(models, coef))
mean_loss <- mean(loss)


# Nonrobust
reps <- 100
data <- vector("list", length = reps)
models <- vector("list", length = reps)
for(i in 1:reps){
  data[[i]] <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0.1, v_outlier = 0) # My idea of vertical outliers
  data_now <- data[[i]]
  models[[i]] <- glm(y ~ X, data = data_now, family = "binomial")
}

# Getting coefs of all models
coefs <- t(sapply(models, coef))


# ROBUST
reps <- 100
data <- vector("list", length = reps)
models <- vector("list", length = reps)
for(i in 1:reps){
  data[[i]] <- binary_reg_dgp1(100, 1, 1, beta = 1, dirty = 0.1, v_outlier = 0) # My idea of vertical outliers
  data_now <- data[[i]]
  models[[i]] <- glmrob(y ~ X, data = data_now, family = "binomial", method =  "WBY")
  # Needs (W)BY otherwise it screws up
}

coefs_rob <- t(sapply(models, coef))

# Plotting
par(mfrow = c(1, 3))
boxplot(x = coefs_clean,
        ylim = c(-3, 3),
        main = "glm fits CLEAN")
boxplot(x = coefs, 
        ylim = c(-3, 3),
        main = "glm fits DIRTY")
boxplot(x = coefs_rob, 
        ylim = c(-3, 3),
        main = "glmrob WBY DIRTY")

## Predictions
# Clean
clean_preds <- predict()


#Hence it can be seen, the coefs for dirty glm() are completely off, glmrob WBY is doing well in comparison with clean glm!
  
###High dimensional case

reps <- 100
data <- vector("list", length = reps)
models <- vector("list", length = reps)

# Trying 1 run:
temp <- binary_reg_dgp1(n = 100, p = 60, p_a = 6, beta = c(rep(1, 6), rep(0, 54)), dirty = 0.05, v_outlier = 0)
temp_test <- binary_reg_dgp1(n = 100, p = 60, p_a = 6, beta = c(rep(1, 6), rep(0, 54)), dirty = 0, v_outlier = 0)

temp_y_test <- temp_test[, 1]
temp_X_test <- as.matrix(temp_test[, -1]) # as.matrix!

temp_y <- temp[, 1]
temp_X <- as.matrix(temp[, -1])

temp_enetlts1 <- enetLTS(xx = temp_X,
                         yy = temp_y,
                         family = "binomial",
                         nfold = 5,
                         repl = 5)
temp_enetlts1_preds <- unname(unlist(predict(temp_enetlts1, newX = temp_X_test, type = "response")))
logit_loss(y_true = temp_y_test, prob_hat = temp_enetlts1_preds)


temp_glmrob1 <- glmrob(y ~ ., data = temp, family = "binomial", method = "BY")

temp_glmrob1_preds <- predict(temp_glmrob1, type = "response")
logit_loss(y_true = temp_y, prob_hat = temp_glmrob1_preds)



for(i in 1:reps){
  data[[i]] <- binary_reg_dgp1(n = 100, p = 60, p_a = 6, beta = c(rep(1, 6), rep(0, 54)), dirty = 0.1, v_outlier = 0) # My idea of vertical outliers
  data_now <- data[[i]]
  models[[i]] <- enetLTS(y ~ X, data = data_now, family = "binomial", method =  "BY")
  # Needs WBY otherwise it screws up
}

temp_y * log(temp_enetlts1_preds) + (1 - temp_y) * log(1 - temp_enetlts1_preds)

#Ofcourse logit loss on a test set WITH OUTLIERS is very high! because in that case loss would be infinity!!!!!!!
  