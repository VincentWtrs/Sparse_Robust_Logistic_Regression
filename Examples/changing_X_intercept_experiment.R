### Effect of changing X measurements (Scaling, standardizing) on the coefficients ###

# Seed for reproducing
set.seed(1234)

## CASE 1: Linear model: No scaling, no centering
# Setting sample size:
N <- 1000

# Generating covariate data
X <- seq(from = -N/2, to = (N/2)-1, by = 1) # -1 because otherwise it's a length = 1001 vector
sd(X) # 288.8194

# Setting true beta
beta0 <- 1
beta1 <- 2

# Specifying DGP for y
y <- beta0 + beta1 * X + rnorm(n = N, mean = 0, sd = 1)

# Plotting
plot(x = X,
     y = y)

# Putting in data frame for ease of use
data <- data.frame(y = y, X = X)

# Fitting linear model
lm1 <- lm(y ~ X, data = data)

# Checking results
lm1

## CASE 2: Linear model: Scaling only
# Scaling data
data$X_scaled <- data$X/sd(data$X)
sd(data$X_scaled) # OK, sd = 1

# Fitting linear model
lm2 <- lm(y ~ X_scaled, data = data)

# Checking results
lm2

## CONCLUSION: If one fits on the standardized data, the coefficient increases as follows: beta_scaled = beta_unscaled * sd(X)
# intercept is still exactly the same!

## CASE 3: Logistic regression, no centering, no scaling
# Setting seed
set.seed(1234)

# Creating covariate data
X <- rnorm(n = N, mean = 0, sd = 5)
sd(X) # 4.986689

# Setting true beta
beta0 <- 1
beta1 <- 2

# Making linear predictor eta (just sum)
eta <- beta0 + beta1 * X

# Creating parameter p of Bernouilli using logistic function
p <- exp(eta)/(1 + exp(eta))

# Creating outcome
y <- rbinom(n = 1000,
            size = 1,
            p = p)
table(y)/N # Quite balanced: 464 0-outcomes

# Plotting
plot(x = X,
     y = y)

# Putting in daa frame
data <- data.frame(y = y, X = X)

# Fitting logistic regression (link = "logit")
logit1 <- glm(y ~ X,
              family = "binomial",
              data = data) # warning about 0/1 is ok
logit1
# might want to check if convergence: logit1$converged == TRUE?

# What about fitting without intercept
logit2 <- glm(y ~ -1 + X,
              family = "binomial",
              data = data)
logit2

## CONCLUSION: Yes the beta1_hat becomes more biased, but not too much. What if the classes are imbalanced a bit more

# Removing some of the y = 0 outcomes
data0 <- data[data$y == 0, ][1:200, ] # So we remove 464-200 = 264 0-outcomes
data1 <- data[data$y == 1, ]
data_unbalanced <- data.frame(rbind(data0, data1))

# Checking proportions
table(data$y) # 0: 200 ; 1 : 536, OK: more unbalanced now

# Fitting again with intercept
logit3 <- glm(y ~ X,
              family = "binomial",
              data = data_unbalanced)

logit3

# Fitting without intercept
logit4 <- glm(y ~ -1 + X,
              family = "binomial",
              data = data_unbalanced)
logit4

## CONCLUSION Seems like even more bias

## CASE 4: Logistic regression, scaling data, no centering
# Checking original data sd
sd(data$X) #  4.986689

# Rescaling
data$X_scaled <- data$X/sd(data$X)
sd(data$X_scaled) # OK = 1

# Fitting
logit5 <- glm(y ~ X_scaled,
              family = "binomial",
              data = data)

logit5

# Compared to original
coef(logit1); coef(logit5)
coef(logit1)["X"] * sd(data$X) ; coef(logit5)["X_scaled"]

## CONCLUSION: Yes the idea for the linear model holds also for the logistic regression model, also the intercept is unchanged

## CASE 5: Logistic regression, no scaling, with centering
# Checking mean X
mean(data$X) # -0.132986

# Centering
data$X_centered <- data$X - mean(data$X)

# Checking sd
sd(data$X_centered); sd(data$X) # ofcourse still exactly the same

# Fitting
logit6 <- glm(y ~ X_centered,
              family = "binomial",
              data = data)
logit6 # Doesn't do much

# What about recentering with + 5
data$X_shifted <- data$X + 5
mean(data$X_shifted) # 4.867014

# Checking sd
sd(data$X); sd(data$X_shifted) # The same ofcourse

# Fitting on shifted data
logit7 <- glm(y ~ X_shifted,
              family = "binomial",
              data = data)
logit7 # But intercept -9.075 now!

# Checking log-likelihood:
logLik(logit6); logLik(logit7) # -127.8722 same, so the models are equivalent

## CASE 7: Standardized data
# Standardizing
data$X_norm <- scale(data$X, center = TRUE, scale = TRUE)

# Fitting
logit8 <- glm(y ~ X_norm,
              family = "binomial",
              data = data)
logit8

coef(logit8)[1] - mean(data$X)/sd(data$X) *  coef(logit8)[2]
coef(logit1)[1]

## Conclusion: yes the formula still holds for the logistic regression model




