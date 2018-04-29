#########################################################
# This file contains a function that given a number of  #
# observations n, and a number of folds K returns       #
# randomly assigned foldids to fold k=1, 2, ..., K      #
# for n ids.                                            #
#########################################################

Kfold <- function(n, K){
  ## Function that assigns n observations to K
  ## folds of roughly equal size.
  foldids <- c()
  
  if (n > K) {
    
    for (k in 1:(ceiling(n/K) - 1)){
      foldids <- c(foldids, 1:K)
    }
    foldids <- c(foldids, sample(x=1:K,
                                 size=n - (ceiling(n/K)*K - K),
                                 replace=FALSE))
  } else if (n==K){
    foldids <- 1:n
  }
  return (sample(x=foldids, size=n, replace=FALSE))  
  
}
