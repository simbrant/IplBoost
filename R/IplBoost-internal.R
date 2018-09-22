
.IplBoost.iter <- function(times, status, mat, betas, lms, w, lambda) {
  ## Internal function that performs one iteration of the 
  ## IplBoost algorithm
  
  # Matrix of risk functions for the coefficients for each landmark point
  risk.s <- exp(mat %*% t(betas))

  # Call Cpp function to compute S0 for each landmark
  S0 <- .compute_S0(as.matrix(risk.s), length(times), length(lms))
  
  # Call C++ functions sequentially to compute S1.j and S2.j for each
  # landmark for each covariate j (loops over j)
  S1 <- lapply(1:dim(mat)[2], .compute_S1_j, risk=as.matrix(risk.s),
               mat=mat, n=length(times), S=length(lms))
  S2 <- lapply(1:dim(mat)[2], .compute_S2_j, risk=as.matrix(risk.s),
               mat=mat, n=length(times), S=length(lms))
  
  # Call C++ functions to sequentially compute the first derivative and the
  # negative of the second derivative for each landmark, for each covariate j
  first.der <- lapply(1:dim(mat)[2],
                      function(j){.compute_u_j(j, status, mat, times, S0,
                                               S1[[j]], length(times),
                                               length(lms), lms, w)})
  
  neg.second.der <- lapply(1:dim(mat)[2],
                           function(j){.compute_negI_j(j, status, times, S0,
                                                       S1[[j]], S2[[j]],
                                                       length(times),
                                                       length(lms), lms, w,
                                                       lambda)})
  
  # Compute scores (proportional to second order Taylor expansion of the ipl)
  score.vars <- as.numeric(lapply(1:dim(mat)[2],
                                  function(j){sum(first.der[[j]]**2/neg.second.der[[j]])}))
  
  # Choose the variable that maximises the approximation
  j.star <- which(score.vars == max(score.vars))

  # Update coefficients
  betas <- betas
  betas[, j.star] <- betas[, j.star] + first.der[[j.star]]/neg.second.der[[j.star]]

  return(betas)
}
