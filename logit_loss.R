logit_loss <- function(y_true, prob_hat){
  # Function that returns the loss from logistic regression (cross entropy loss)
  # Note: the loss is normalized (divided by the sample size) making it comparable for differing sizes
  
  # Renaming arguments for clearer code
  y <- y_true
  p <- prob_hat
  
  loss <- - (y * log(p) + (1 - y) * log(1 - p))
  
  # Making correction if the prediction is spot on (0 if 0, 1 if 1) to avoid NaN
  loss[which(p == y)] <- 0
  
  # Calculating average
  avg_loss <- mean(loss)
  
  # Output
  return(list(loss = loss, avg_loss = avg_loss))
}