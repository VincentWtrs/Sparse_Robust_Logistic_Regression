### Showing overfitting in logistic regression using polynomials ###

# Setting seed
set.seed(1234)

### Generating data
## Class 1 (y = 0)
# x1_1: X (x1) coordinate (y = 0)
x1_1 <- rnorm(n = 25,
             mean = 1,
             sd = 0.5)
# x2_2: y (x2) coordinate (y = 0)
x2_1 <- rnorm(n = 25,
             mean = 1,
             sd = 1)

## Class 2 (y = 1)
# x1_2: X (x1) coordinate (y = 1)
x1_2 <- rnorm(n = 25,
             mean = 2,
             sd = 0.5)
# x2_2: y (x2) coordinate (y  = 1)
x2_2 <- rnorm(n = 25,
             mean = 2,
             sd = 1)

## Combining into dataset
x1 <- c(x1_1, x1_2)
x2 <- c(x2_1, x2_2)
y <- c(rep(0, 25), rep(1, 25))
data <- data.frame(y = y, x1 = x1, x2 = x2)


## Plotting
plot(x = x1_1, 
     y = x2_1, 
     xlim = c(0, 5), 
     ylim = c(0, 5), 
     pch = 1, 
     xlab = "x1", 
     ylab = "x2", 
     cex = 1.5)
points(x = x1_2, 
       y = x2_2, 
       pch = 4, 
       cex = 1.5)
legend(x = "topright",
       legend = c("y = 0",
                  "y = 1"),
       pch = c(1, 4))


### Fitting
## First-order logit
logit1 <- glm(y ~ x1 + x2, data = data, family = "binomial")
logit1

# Getting coefs
coef1 <- coef(logit1)

# Grid for plotting
x1_grid <- seq(-5, 5, by = 0.001)

# Calculating x2 based on x1 (transformed from 1 = b0 + b1*x1 + b2*x2)
x2_grid <- (1 - coef1[1] - coef1[2] * x1_grid)/coef1[2] 

# Adding decision boundary
lines(x = x1_grid, y = x2_grid)

## Fitting 2nd order model
logit2 <- glm(y ~ x1 + x2 + I(x1^2) + I(x2^2), family = "binomial", data = data)
logit2


# Other way of working necessary for higher order models: making grid and plotting contour
x1_grid <- seq(from = -2, to = 6, by = 0.01)
x2_grid <- seq(from = -2, to = 6, by = 0.01)
grid <- data.frame(expand.grid(x1 = x1_grid, 
                               x2 = x2_grid))

# Predictions
preds2 <- predict(object = logit2,
                      newdata = grid,
                      type = "response")

# Adding contour to existing plot
contour(x = x1_grid, 
        y = x2_grid, 
        z = matrix(preds2, nrow = length(x1_grid)), 
        levels = 0.5, 
        col = "blue",
        drawlabels = FALSE,
        add = TRUE)


## 7-th order logit
logit3 <- glm(y ~ x1 + x2 + I(x1^2) + I(x2^2) + I(x1^3) + I(x2^3) + I(x1^4) + I(x2^4) + I(x1^5) + I(x2^5) + I(x1^6) + I(x2^6) + I(x1^7) + I(x2^7), 
              family = "binomial",
              control = glm.control(maxit = 200),
              data = data)

# Predictions 
preds3 <- predict(object = logit3,
                 newdata = grid,
                 type = "response")

# Contour
contour(x = x1_grid, 
        y = x2_grid, 
        z = matrix(preds3, nrow = length(x1_grid)), 
        levels = 0.5, 
        col = "red",
        drawlabels = FALSE,
        add = TRUE)