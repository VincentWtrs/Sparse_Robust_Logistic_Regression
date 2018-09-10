### Small example showing the negative loglikelihood loss for logistic regression ###

# Defining simple logistic loss function
logloss <- function(y, pi){
  loss <- -(y * log(pi) + (1 - y) * log(1 - pi))/length(y)
}

# Creating a fine grid of fictional predictions pi within unit interval
pi_grid <- seq(from = 0.01,
               to = 0.99,
               by = 0.001)

# Creating vector of outcomes (0 and 1 respectively)
y0 <- rep(0, length = length(pi_grid))
y1 <- rep(1, length = length(pi_grid))

# Calculating loss based on the predictions and true outcomes
loss0 <- logloss(y = y0,
                 pi = pi_grid)

loss1 <- logloss(y = y1,
                 pi = pi_grid)

## Plotting
# For 0-outcomes (BLACK)
plot(x = pi_grid,
     y = loss0,
     type = "l",
     lty = 1,
     lwd = 1,
     col = "black",
     xlab = "Predicted Probability",
     ylab = "Negative log-likelihood loss")
# For 1-outcomes (RED)
lines(x = pi_grid,
      y = loss1,
      lty = 1,
      lwd = 1,
      col = "red")
# Legend
legend(x = "topright",
       legend = c("y = 0",
                  "y = 1"),
       col = c("black", "red"),
       lty = 1)

## IN TERMS OF LINEAR PREDICTOR
eta_grid <- log(pi_grid/(1 - pi_grid))

## Plotting
# For 0-outcomes (BLACK)
plot(x = eta_grid,
     y = loss0,
     type = "l",
     lty = 1,
     lwd = 1,
     col = "black",
     xlab = "Linear predictor value",
     ylab = "Negative log-likelihood loss")
# For 1-outcomes (RED)
lines(x = eta_grid,
      y = loss1,
      lty = 1,
      lwd = 1,
      col = "red")
# Legend
legend(x = "topright",
       legend = c("y = 0",
                  "y = 1"),
       col = c("black", "red"),
       lty = 1)