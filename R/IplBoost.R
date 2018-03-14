##############################################################################################
## This file contains the main functionionality for the package IplBoost, both the function ##
## IplBoost() for fitting models, and cv.IplBoost for tuning the number of iterations via   ##
## cross-validataion                                                                        ##
##############################################################################################


IplBoost <- function(times, status, mat, lms, w, M, lambda, verbose=FALSE, standardise=TRUE,
                     compute.ipl=TRUE){
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
  
  if (compute.ipl){
    ipl.vals[[1]] <- .compute_ipl(times, status, mat, as.matrix(estimates[[1]]), lms, w, length(lms), length(times), dim(mat)[2])
  }
  
  # Loop to the given number of iterations, update the coefficients
  for (m in 2:(M+1)){
    if(verbose){
      cat(c("m: ", m-1, "\n"))
    }
    
    # Update estimates
    estimates[[m]] <- .IplBoost.iter(times, status, mat, estimates[[m-1]], lms, w, lambda)

    # Compute ipl
    if (compute.ipl){
      ipl.vals[[m]] <- .compute_ipl(times, status, mat, as.matrix(estimates[[m]]),
                                    lms, w, length(lms), length(times), dim(mat)[2])
    }
  }
  
  if (standardise) {
    for (m in 2:M){
      estimates[[m]] <- as.matrix(estimates[[m]])
      #.scale_columns(estimates[[m]], sds, dim(estimates[[m]])[1], dim(estimates[[m]])[2])
      for (j in 1:length(sds)){
        estimates[[m]][, j] <-  estimates[[m]][, j]/sds[j]
      }
    }
  }
  
  return(list(estimates=estimates, ipl=ipl.vals))
}


cv.IplBoost <- function(times, status, mat, lms, w, M, lambda, folds, verbose=FALSE,
                        standardise=TRUE, parallel=FALSE){

  ## This function performs K-fold cross-valitation for IplBoost,
  ## to tune the number of iterations
  
  # Make sure the observations are in increasing order
  status <- status[order(times)]
  mat <- mat[order(times), ]
  times <- times[order(times)]
  
  if(standardise){
    mat <- scale(mat)
  }
  
  if(parallel){
    if (sfParallel()) {
      cat("IplBoost Running in parallel mode on", sfCpus(), "nodes.\n")
      sfExportAll()
    }
    else{
      cat("IPLBOOST Running in sequential mode.\n")
    }
  }
  
  # Estimate the K times M models
  if (parallel){
    cv.mods <- sfLapply(1:max(folds), function(k){print(k);IplBoost(times=times[folds!=k],
                                                                    status=status[folds!=k],
                                                                    mat=mat[folds!=k, ], lms=lms, w=w,
                                                                    M=M, lambda=lambda, standardise=FALSE,
                                                                    compute.ipl=FALSE)})
  } else {
    cv.mods <- lapply(1:max(folds), function(k){print(k);IplBoost(times=times[folds!=k],
                                                                  status=status[folds!=k],
                                                                  mat=mat[folds!=k, ], lms=lms, w=w,
                                                                  M=M, lambda=lambda, standardise=FALSE,
                                                                  compute.ipl=FALSE)})
  }
  
  
  # Compute the cross validated integrated partial likelihood
  
  # Define two functions to perform the double for-loop over the models from each
  # iteration for each fold
  cv.ipl.k <- function(betas, times, status, mat, lms, w){
    .compute_ipl(times=times, status=status, mat=mat,
                 betas=betas, lms=lms, w=w, S=length(lms),
                 n=length(times), p=dim(mat)[2])
  }
  
  cv.ipl.m <- function(m, folds, times, status, mat, mods, lms, w){
    cvs <- lapply(1:max(folds), function(k) cv.ipl.k(betas=as.matrix(mods[[k]]$estimates[[m+1]]),
                                                     times=times[folds==k], status=status[folds==k],
                                                     mat=mat[folds==k, ], lms=lms, w=w))
    return(mean(as.numeric(cvs)))
  }

  # Compute the ipl itself by nested lapply calls
  ipl.cv <- as.numeric(lapply(0:M, cv.ipl.m, folds=folds, times=times, status=status, mat=mat,
                              mods=cv.mods, lms=lms, w=w))

   return(list(ipl.cv=ipl.cv, opt.m = which(ipl.cv == max(ipl.cv)) - 1))
}
