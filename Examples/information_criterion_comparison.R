# ### EXAMPLE SCRIPT: INFORMATION CRITERION VISUALIZED ###
# 
# # Generating data (load binary_reg_dgp2() if needed)
# data <- binary_reg_dgp2(n = 100,
#                         p = 50,
#                         p_a = 5,
#                         beta = c(rep(1, 5), rep(0, 45)),
#                         dirty = 0,
#                         type = "bernoulli")
# 
# data_y <- data.matrix(data[, 1])
# data_X <- data.matrix(data[, -1]) # 100 x 50 matrix (intercept is not included)
#  
# # Fitting LASSO models (over path of lambdas)
# models <- glmnet(x = data_X,
#                  y = data_y,
#                  family = "binomial",
#                  alpha = 1)
# 
# # Getting predictions (on training data)
# preds <- predict(object = models,
#                  newx = data_X,
#                  type = "response")
# 
# # Getting lambdas used in the fitting path
# lambdas <- models$lambda
# 
# # Checking how many models were fitted
# L <- length(lambdas)
# L # 92
# 
# # Getting degrees of freedom for each model (= number of nonzero coefs)
# dfs <- models$df
# 
# # Getting fitted coefs
# betas <- as.matrix(models$beta) # as.matrix() to get rid of sparse structure
# dim(betas) # 50 x 92: Each column contains the beta vector for the respective lambda value
# 
# # Getting average negative loglikelihood loss
# loss <- apply(preds, MARGIN = 2, FUN = function(x) logit_loss(y_true = data_y,
#                                                               prob_hat = x)$avg_loss)
# 
# # Selecting some information criteria to use
# ic <- c("AIC", "AIC_C", "EBIC", "EBIC2", "GIC", "ERIC")
# 
# for(i in 1:length(ic)){
#   for(j in 1:L){
#     
#   }
# }
# 
# penalty <- data.frame(matrix(NA, nrow = L, ncol = length(ic)))
# names(penalty) <- ic
# for(i in 1:L){
#   lambda_now <- lambdas[i]
#   model_now <- glmnet(x = data_X,
#                       y = data_y,
#                       family = "binomial",
#                       alpha = 1,
#                       lambda = lambda_now)
#   for(j in 1:length(ic)){
#     penalty[i, j] <- ic_penalty(type = ic[j],
#                           model = model_now,
#                           X = data_X,
#                           alpha = 1,
#                           EBIC_sigma = 1,
#                           intercept = FALSE) # put to false for testing atm
#   }
# }
# 
# 
# aic <- loss + penalty$AIC
# aicc <- loss + penalty$AIC_C
# ebic <- loss + penalty$EBIC
# ebic2 <- loss + penalty$EBIC2
# gic <- loss + penalty$GIC
# eric <- loss + penalty$ERIC
# 
# overview <- data.frame(aic, aicc, ebic, ebic2, gic, eric)
# 
# plot(x = dfs, 
#      y = overview$aic,
#      col = 1,
#      type = "l",
#      lwd = 2,
#      ylim = c(0, 2),
#      xlim = c(0, 10))
# for(i in 2:ncol(overview)){
#   lines(x = dfs,
#         y = overview[, i],
#         col = i,
#         lwd = 2)
# }

# SEEMS TO BE SOMETHING WRONG WITH THIS



models <- glmnet(x = train1_X,
                 y = train1_y,
                 family = "binomial",
                 alpha = 1)

# Getting predictions (on training data)
preds <- predict(object = models,
                 newx = data_X,
                 type = "response")

# Getting lambdas used in the fitting path
lambdas <- models$lambda

# Checking how many models were fitted
L <- length(lambdas)
L # 94

# Getting degrees of freedom for each model (= number of nonzero coefs)
dfs <- models$df

# Getting fitted coefs
betas <- as.matrix(models$beta) # as.matrix() to get rid of sparse structure
dim(betas) # 50 x 92: Each column contains the beta vector for the respective lambda value

# Getting average negative loglikelihood loss
loss <- apply(preds, MARGIN = 2, FUN = function(x) logit_loss(y_true = data_y,
                                                              prob_hat = x)$avg_loss)



## EBIC
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "EBIC",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE,
                           EBIC_sigma = 0)
}

ebic <- loss + penalty
plot(x = dfs, y = ebic, type = "l", main = "Amount of covariates versus EBIC", ylim = c(0, 1.5), xlim = c(0, 20), xaxt = "none")
axis(side = 1, at = seq(1:20))
abline(v = 5, lwd = 1, lty = 2)
dfs[which.min(ebic)]
betas[, which.min(ebic)]


## BIC
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "BIC",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

bic <- loss + penalty
lines(x = dfs, y = bic, type = "l", lty = 1, lwd = 2, col = "red")
dfs[which.min(bic)]
betas[, which.min(bic)]

## EBIC2
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "EBIC2",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

ebic2 <- loss + penalty
lines(x = dfs, y = ebic2, type = "l", lty = 1, lwd = 2, col = "blue")
dfs[which.min(ebic2)]
betas[, which.min(ebic2)]

## ERIC
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "ERIC",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

eric <- loss + penalty
lines(x = dfs, y = eric, type = "l", lty = 1, lwd = 2, col = "green")
dfs[which.min(eric)]
betas[, which.min(eric)]

## AIC
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "AIC",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

aic <- loss + penalty
lines(x = dfs, y = aic,  type = "l", lty = 1, lwd = 2, col = "purple")
dfs[which.min(aic)]
betas[, which.min(aic)]


## BIC-HD
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "BIC_HD",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

bichd <- loss + penalty
lines(x = dfs, y = bichd,  type = "l", lty = 1, lwd = 2, col = "yellow")
dfs[which.min(bichd)]
betas[, which.min(bichd)]

## GIC
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "GIC",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

gic <- loss + penalty
lines(x = dfs, y = gic,  type = "l", lty = 1, lwd = 2, col = "pink")
dfs[which.min(gic)]
betas[, which.min(gic)]

# BIC_WLL
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "BIC_WLL",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

bicwll <- loss + penalty
lines(x = dfs, y = bicwll,  type = "l", lty = 1, lwd = 2, col = "gray")
dfs[which.min(bicwll)]
betas[, which.min(bicwll)]


## CV
cverror <- cv.glmnet(x = train1_X,
                     y = train1_y,
                     family = "binomial",
                     lambda = lambdas,
                     alpha = 1,
                     type.measure = "deviance",
                     nfolds = 20)$cvm

lines(x = dfs, y = cverror, type = "l", lty = 1, lwd = 2, col = "brown")

dfs[which.min(cverror)]
betas[, which.min(cverror)]

# AIC_C
penalty <- numeric(length = L)
for(i in 1:L){
  model_now <- glmnet(x = train1_X,
                      y = train1_y,
                      family = "binomial",
                      alpha = 1,
                      lambda = lambdas[i])
  penalty[i] <- ic_penalty(type = "AIC_C",
                           model = model_now,
                           X = train1_X,
                           alpha = 1,
                           intercept = FALSE)
}

aicc <- loss + penalty
lines(x = dfs, y = aicc,  type = "l", lty = 1, lwd = 2, col = "orange")
dfs[which.min(aicc)]
betas[, which.min(aicc)]

