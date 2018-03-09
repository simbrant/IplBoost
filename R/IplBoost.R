##############################################################################################
## This file contains the main functionionality for the package IplBoost, both the function ##
## IplBoost() for fitting models, and cv.IplBoost for tuning the number of iterations via   ##
## cross-validataion                                                                        ##
##############################################################################################


IplBoost <- function(times, status, mat, lms, w, M, lambda, verbose, standardise=TRUE){
  ## This is the main function of the package, that fits sliding landmark models
  ## by boosting van Houwelingens integrated partial likelihood, following the strategy
  ## of CoxBoost.

  # Check input
  if (length(lambda) == 1){
    lambda = rep(lambda, length(lms))
  } else if (!(length(lms)==length(lambda))) {
    stop("Must supply a vector of lambdas matching the vector of landmarks")
  }
  
  if(standardise){
    sds <- apply(mat, 2, function(x){sqrt(var(x))}) 
    mat <- scale(mat)
  }
  
  # Make sure the observations are in increasing order
  status <- status[order(times)]
  mat <- mat[order(times), ]
  times <- times[order(times)]
  
  # Create list of martrices of estimates and a vector of ipl-values, initialise the first entry
  estimates <- vector("list", M+1)
  ipl.vals <- vector("numeric", M+1)
  estimates[[1]] <- sparseMatrix(i=c(1), j=c(1), x = c(0), dims = c(length(lms), dim(mat)[2]))
  ipl.vals[[1]] <- .compute_ipl(times, status, mat, as.matrix(estimates[[1]]), lms, w, length(lms), length(times), dim(mat)[2])
  
  # Loop to the given number of iterations, update the coefficients
  for (m in 2:(M+1)){
    if(verbose){
      cat(c("m: ", m-1, "\n"))
    }
    estimates[[m]] <- .IPLBOOST.iter(times, status, mat, estimates[[m-1]], lms, w, lambdas)
    ipl.vals[[m]] <- compute_ipl(times, status, mat, as.matrix(estimates[[m]]),
                                 lms, w, length(lms), length(times), dim(mat)[2])
  }
  
  if (standardise) {
    for (m in 2:M){
      estimates[[m]] <- as.matrix(estimates[[m]])
      .scale_columns(estimates[[m]], sds, dim(estimates[[m]])[1], dim(estimates[[m]])[2])
    }
  }
  
  return(list(estimates=estimates, ipl=ipl.vals))
}


cv.IPLBOOST <- function(times, status, mat, lms, w, M, lambda, verbose, folds, standardise=TRUE){
  ## This function performs K-fold cross-valitation for IplBoost, to tune the number of
  ## Iterations
  
  if(standardise){
    mat <- scale(mat)
  }
  
  # Estimate the K times M models
  cv.mods <- lapply(1:max(folds), function(k){IplBoost(times[folds!=k], status[folds!=k], mat[folds!=k, ],
                                                       lms, w, M, lambda, FALSE, FALSE)})
  
  # compute the Cross-validated ipl for each iteration
  ipl.cv <- vector("numeric", M+1)
  for (m in 0:M){
    ipl.curr <- vector("numeric", max(folds))
    for (k in 1:max(folds)){
      ipl.curr[k] <- compute_ipl(t[folds==k], d[folds==k], mat[folds==k, ],
                                 as.matrix(cv.mods[[k]][[1]][[m+1]]), lms, w,
                                 length(lms), length(t[folds==k]), dim(mat)[2])
    }
    ipl.cv[m+1] <- mean(ipl.curr)
  }
  
  return(list(ipl.cv=ipl.cv, opt.m = which(ipl.cv == max(ipl.cv))))
}
