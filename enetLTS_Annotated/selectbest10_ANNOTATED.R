selectbest10_ANNOTATED <- function (x, y, family, h, hsize, alpha, lambda, nsamp, s1, para, ncores, scal, seed) 
{
  # Initialization
  obj <- NULL # ?
  all_subsets <- enetLTS:::InitialSubset(x = x, 
                                         y = y, 
                                         family = family, 
                                         h = h, 
                                         hsize = hsize, 
                                         alpha = alpha, 
                                         lambda = lambda, 
                                         nsamp = nsamp, 
                                         para = para, 
                                         ncores = ncores, 
                                         scal = scal, 
                                         seed = seed)
  # all_subsets: list of length 2: $subsets and $index.subsets
  # all_subsets$subsets: list of length 500 (nsamp)
  # all_subsets$index.subsets: list of 2000 (nsamp x 4 (4: size of elemental subset?))
    
  subsets <- all_subsets$subsets
  index.subsets <- all_subsets$index.subsets 
  # subsets <- all_subsets$subsets each element contains 2 elements: $obj: value of objective function (Bianco-Yohai)
  # and $indx: the indices of h (75) elements (TO DO ???)
  # index.subsets <- all_subsets$index.subsets is a matrix of 4 x nsamp (TO DO ???)
  
  ## Gathering values of objective function (just extracting it basiacally) for all nsamp (500) subsets
  # Parallel case
  if (para) {
    obj <- unlist(parallel:::mclapply(1:nsamp, function(ob, sub) {
      ob_val <- subsets[[ob]]$obj
    }, subsets, mc.cores = ncores, mc.allow.recursive = FALSE))
    
     # Non-parallel case
  } else {
    for (i in 1:nsamp) { # For each of the nsamp (500) elemental subsets
      obj <- c(obj, subsets[[i]]$obj) # Just keeping adding the subsets objective function
      # Basically a nasty way of doing obj[i] <- subsets[[i]]$obj
    }
  } # At the end of the loop obj is a vector of length = nsamp (500) with values for the obj function (Bianco-Yohai)
  
  ## Sorting the nsamp (500) objective function values
  # Binomial Case
  if (family == "binomial") {
    obj_sorted <- sort(obj, decreasing = TRUE, index.return = TRUE) # Higher is better!!
    
  # Gaussian Case
  } else if (family == "gaussian") {
    obj_sorted <- sort(obj, decreasing = FALSE, index.return = TRUE) # MSE: lower is better!!
  } # obj°sorted is still vector of length = nsamp (500) with $x: objective, $ix: indices
  
  ## Gathering the s1 (10) highest objective functions
  obj <- obj_sorted$x[1:s1] # Getting 10 highest (may include inf)
  s1_new <- length(obj[!is.infinite(obj)]) # s1_new: s1 excluding the infs (will be equal or smaller than s1!)
  idx <- obj_sorted$ix[1:s1_new] # Gathering the INDICES (ids) OF THE BEST 10 subsets
  
  # Some error throwing if all s1 models had inf, hence s1_new == 0:
  if (s1_new == 0) {
    stop(paste("Model is not suitable for alpha", alpha, 
               "lambda", lambda, "for this data set. Choose another lambda."))
  }
  
  ## Gathering the h!!! observation indices for the best s1_new (10, often) best subsets
  # Parallel case
  if (para) {
    bestindex <- parallel:::mclapply(1:s1_new, function(c, idx, subsets) {
      indx <- subsets[[idx[c]]]$indx
    }, idx, subsets, mc.cores = ncores)
    
  # Non-parallel case
  } else {
    bestindex <- lapply(1:s1_new, function(c, idx, subsets) {
      indx <- subsets[[idx[c]]]$indx
    }, idx, subsets)
  } # bestindex is a list of length s1_new (10) with each element having index of h (75) observations
  return(list(idxbest = bestindex, s1_new = s1_new, subsets = subsets, index.subsets = index.subsets))
}