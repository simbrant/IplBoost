##############################################################################################
## This file contains the main functionionality for the package IplBoost, both the function ##
## IplBoost() for fitting models, and cv.IplBoost for tuning the number of iterations via   ##
## cross-validataion                                                                        ##
##############################################################################################

IplBoost <- function(times, status, mat, landmarks, w, M, lambda, ...) UseMethod("IplBoost")


IplBoost.default <- function(times, status, mat, landmarks, w, M, lambda, verbose=FALSE,
                             standardise=TRUE, compute.ipl=TRUE){
  ##' IplBoost
  ##' @description This is the main function of the package, that fits sliding landmark models
  ##' by boosting van Houwelingens integrated partial likelihood, following the strategy
  ##' of \link[CoxBoost]{CoxBoost}.
  ##' @param times A n-dimensional vector of survival times
  ##' @param status A n-dimensional vector of censoring indicators
  ##' @param mat A n x p matrix of covariate values
  ##' @param landmarks A S-dimensinal vector of landmark points to produce dynamic 
  ##' predictions from
  ##' @param w A number. The "landmark interval width" or how far ahead the 
  ##' survival predictions will be made
  ##' @param M A number. The number of boosting iterations to perform
  ##' @param lambda A number or an S-dimensional vector. The regularisation
  ##'  parameter for the boosting algorithm.
  ##' @param verbose Boolean. Indicates whether the iteration number should
  ##'  be printed to the console.
  ##' @param standardise Boolean. Indicates if covariates should be standardised prior to
  ##' fitting the model. Coefficient estimates are transformed back to the scale of the data if true.
  ##' @param compute.ipl Boolean. Indicated if the ipl should be computed for each step.
  ##' @return A list of (M + 1) elements containing the landmark coefficients as a (S x p)
  ##' matrix for each iteration.
  ##' @return A vector of the ipl computed for each step.
  ##' @examples
  ##' # Tune the number of iterations via cross validation (see \link{cv.Iplboost})
  ##' cv.mod <- cv.IplBoost(times, status, design, landmarks=seq(0, 10, 0.1),
  ##'                       w=5, M=100, lambda=100,
  ##'                       folds=Kfold(length(times), 10))
  ##' # Fit the model using the tuned number of iterations
  ##' mod <- IplBoost(times, status, design, landmarks=seq(0, 10, 0.1),
  ##'                 w=5, M=cv.mod$opt.m, lambda=100)
  ##' estimates <- mod$estimates[[cv.mod$opt.m + 1]]

  
  # Check input
  if (length(lambda)==1){
    lambda <- rep(lambda, length(landmarks))
  } else if (!(length(landmarks)==length(lambda))) {
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
  
  # Add random noise to break ties if nescessary
  if (length(times) != length(unique(times))){
    times <- times + abs(rnorm(length(times), 0, 10**(-8)))
  }
  
  # Create list of martrices of estimates and a vector of ipl-values, initialise
  # the first entry
  estimates <- vector("list", M+1)
  ipl.vals <- vector("numeric", M+1)
  estimates[[1]] <- sparseMatrix(i=c(1), j=c(1), x=c(0), dims=c(length(landmarks),
                                                                dim(mat)[2]))
  
  if (compute.ipl){
    ipl.vals[[1]] <- .compute_ipl(times, status, mat, as.matrix(estimates[[1]]),
                                  landmarks, w, length(landmarks), length(times),
                                  dim(mat)[2])
  }
  
  # Loop to the given number of iterations, update the coefficients
  for (m in 2:(M+1)){
    if(verbose){
      cat(c("m: ", m-1, "\n"))
    }
    
    # Update estimates
    estimates[[m]] <- .IplBoost.iter(times, status, mat, estimates[[m-1]], landmarks,
                                     w, lambda)

    # Compute ipl
    if (compute.ipl){
      ipl.vals[[m]] <- .compute_ipl(times, status, mat, as.matrix(estimates[[m]]),
                                    landmarks, w, length(landmarks), length(times),
                                    dim(mat)[2])
    }
  }
  
  if (standardise) {
    for (m in 2:(M + 1)){
      estimates[[m]] <- as.matrix(estimates[[m]])
      .scale_columns(estimates[[m]], sds, dim(estimates[[m]])[1],
                     dim(estimates[[m]])[2])
    }
  }

  result <- list(estimates=estimates, ipl=ipl.vals, call=match.call())
  class(result) <- "IplBoost"
  
  return(result)
  
}


cv.IplBoost <- function(times, status, mat, landmarks, w, M, lambda, folds, ...) UseMethod("cv.IplBoost")

cv.IplBoost.default <- function(times, status, mat, landmarks, w, M, lambda, folds,
                                verbose=FALSE, standardise=TRUE, parallel=FALSE){
  ##' cv.IplBoost
  ##' @description This function performs K-fold cross-valitation for \link{IplBoost},
  ##' to tune the number of iterations.
  ##' @param times A n-dimensional vector of survival times
  ##' @param status A n-dimensional vector of censoring indicators
  ##' @param mat A n x p matrix of covariate values
  ##' @param landmarks A S-dimensinal vector of landmark points to produce dynamic
  ##' predictions from
  ##' @param w A number. The "landmark interval width" or how far ahead the survival
  ##'  predictions will be made
  ##' @param M A number. The number maximum of boosting iterations to perform
  ##' @param lambda A number or an S-dimensional vector. The regularisation parameter 
  ##' for the boosting algorithm.
  ##' @param folds A n-dimensional vector that assigns a number between 1 and K to 
  ##' each observations for cross-validataion.
  ##' @param verbose Boolean. Indicates whether the iteration number should be 
  ##' printed to the console.
  ##' @param standardise Boolean. Indicates if covariates should be standardised.
  ##' @param parallel Boolean. Indicates if the cross validation should be performed
  ##' in parallel. Relies on the package snowfall. Cluster must be initialised before
  ##' calling cv.IplBoost in the case of this being true, see the
  ##' \link[snowfall]{snowfall} package and \link[snowfall]{sfInit}.
  ##' @return A list containg a (M+1)-dimensional vector of the cross validated ipl
  ##' as the first element, and a number indicating the optimal number of iterations
  ##' as the second.
  ##' @examples 
  ##' # Tune the number of iterations via cross validation
  ##' cv.mod <- cv.IplBoost(times, status, design, lamdmarks=seq(0, 10, 0.1),
  ##'                       w=5, M=100, lambda=100,
  ##'                       folds=Kfold(length(times), 10))
  ##'                       
  ##' # Plot the cross-validated likelihood
  ##' plot(cv.mod)
  ##' 
  ##' # To tune in parallel on 2 cores:
  ##' library(snowfall)
  ##' sfInit(parallel=TRUE, cpus=2)
  ##' cv.mod <- cv.IplBoost(times, status, design, landmarks=seq(0, 10, 0.1),
  ##'                       w=5, M=100, lambda=100, folds=Kfold(length(times), 10),
  ##'                       parallel=TRUE)
  ##' sfStop()
  ##' 

  # Make sure the observations are in increasing order
  status <- status[order(times)]
  mat <- mat[order(times), ]
  times <- times[order(times)]
  
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
                                                                    mat=mat[folds!=k, ],
                                                                    landmarks=landmarks,
                                                                    w=w, M=M, lambda=lambda,
                                                                    standardise=standardise,
                                                                    compute.ipl=FALSE)})
  } else {
    cv.mods <- lapply(1:max(folds), function(k){print(k);IplBoost(times=times[folds!=k],
                                                                  status=status[folds!=k],
                                                                  mat=mat[folds!=k, ],
                                                                  landmarks=landmarks,
                                                                  w=w, M=M, lambda=lambda,
                                                                  standardise=standardise,
                                                                  compute.ipl=FALSE)})
  }

  # Compute the cross validated integrated partial likelihood
  
  # Define two functions to perform the double for-loop over the models from each
  # iteration for each fold
  cv.ipl.k <- function(betas, times, status, mat, landmarks, w){
    .compute_ipl(times=times, status=status, mat=mat,
                 betas=betas, lms=landmarks, w=w, S=length(landmarks),
                 n=length(times), p=dim(mat)[2])
  }
  
  cv.ipl.m <- function(m, folds, times, status, mat, mods, landmarks, w){
    cvs <- lapply(1:max(folds), function(k) cv.ipl.k(betas=as.matrix(mods[[k]]$estimates[[m+1]]),
                                                     times=times[folds==k], status=status[folds==k],
                                                     mat=mat[folds==k, ], landmarks=landmarks, w=w))
    return(mean(as.numeric(cvs)))
  }

  # Compute the ipl itself by nested lapply calls
  ipl.cv <- as.numeric(lapply(0:M, cv.ipl.m, folds=folds, times=times, status=status,
                              mat=mat, mods=cv.mods, landmarks=landmarks, w=w))

  
  result <- list(ipl.cv=ipl.cv, opt.m = (which(ipl.cv == max(ipl.cv)) - 1)[1],
                 call=match.call())
  class(result) <- "cv.IplBoost"
  
  return(result)
}

plot.cv.IplBoost <- function(x, ...){
  plot(0:(length(x$ipl.cv) - 1), x$ipl.cv, lwd=2, "l", main="Cross-validated ipl",
       ylab="Integrated partial likelihood", xlab="Number of iterations")
  abline(v=x$opt.m, lty=3)
  legend(x="bottomright", legend = c("ipl", "Opt. m"), col = c(1,1), lty = c(1, 3),
         lwd=c(2, 1))
}
