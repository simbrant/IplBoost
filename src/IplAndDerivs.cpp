/*
 This file contains internal C++ functions to compute the first and second derivatives of van Houwelingens
 ipl for one covariate index, for all landmarks, as well as evaluate the ipl itself.
*/

#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_u_j(int j, NumericVector status, NumericMatrix mat, NumericVector times,
                          NumericMatrix S0, NumericMatrix S1j, int n, int S, NumericVector lms,
                          double w){
  /*
   Computes u_j = d/dbeta_{js} ipl(beta_s), s = 1,..., S
  
  Args:
  - j is an integer indicating the column to compute the derivatives for
  - times is an n-dimensional vector of lifetimes
  - status is an n-dimensional vector of censoring indicators
  - mat is the design matrix (n x p)
  - n i the length of times (number of observations)
  - S is the number of rows of risk (number of landmarks)
  - lms is an S-dimensional vector of landmark points
  - w is the width of each landmark interval
  - S0, and S1j are (n x S) matrices as computed by compute_S0 and compute S1 
  returns a (n x S) matrix
  */
  
  // Empty vector to fill with the first derivatives
  NumericVector u(S, 0.0);
  
  int lmsind = 0;
  for (int s=0; s<S; s++){
    // Assume times are ordered
    int i = lmsind;
    while (times[i] < lms[s] & i < n){
      i += 1;
    }
    lmsind = i;
    while (times[i] <= lms[s] + w & i < n){
      u[s] += status[i]*(mat(i, j-1) - S1j(i, s)/S0(i, s));
      i += 1;
    } 
  }   
  return u;
}


// [[Rcpp::export]]
NumericVector compute_negI_j(int j, NumericVector status, NumericVector times, NumericMatrix S0,
                             NumericMatrix S1j, NumericMatrix S2j, int n, int S, NumericVector lms,
                             double w, NumericVector lambda){
  /*
   Computes -I_j = d**2/dbeta_{js}**2 ipl(beta_s), s = 1,..., S
  
  Args:
  - j is an integer indicating the column to compute the derivatives for
  - times is an n-dimensional vector of lifetimes
  - status is an n-dimensional vector of censoring indicators
  - n i the length of times (number of observations)
  - S is the number of rows of risk (number of landmarks)
  - lms is an S-dimensional vector of landmark points
  - w is the width of each landmark interval
  - S0, S1j and S2j are (n x S) matrices as computed by compute_S0, compute S1_j and compute_S2_j
  returns a (n x S) matrix
  */
  
  NumericVector negI(S, 0.0);
  
  int lmsind = 0;
  for (int s=0; s<S; s++){
    // Assume times are ordered
    int i = lmsind;
    while (times[i] < lms[s] & i < n){
      i += 1;
    }
    lmsind = i;
    while (times[i] <= lms[s] + w & i < n){
      negI[s] += status[i]*((S2j(i, s)*S0(i, s)- S1j(i, s)*S1j(i, s))/(S0(i, s)*S0(i, s)));
      i += 1;
    } 
  }
  for (int s = 0; s<S; s++){
    negI[s] += lambda[s];
  }
  return I;
}


//[[Rcpp::export]]
double compute_ipl(NumericVector times, NumericVector status, NumericMatrix mat, NumericMatrix betas,
                   NumericVector lms, double w, int S, int n, int p) {
  double ipl = 0.0;
  NumericMatrix prognostic(n, S);
  NumericMatrix s0(n, S);
  
  // Compute matrix of prognostic indexes for each landmark
  for (int s=0; s < S; s++){
    for (int i=0; i<n; i++){
      for (int j=0; j<p; j++){
        prognostic(i, s) = prognostic(i, s) + mat(i, j)*betas(s, j);
      }
    }
  }
  
  // Compute matrix containing S0 (one column for each landmark)
  for (int s=0; s<S; s++){
    
    for (int i=0; i<n; i++){
      s0(i, s) = exp(prognostic(i, s));
    }
    
    for (int i = n-2; i >= 0; i--){
      s0(i, s) = s0(i+1, s) + s0(i, s);
    }
  }
  
  // Compute the ipl
  for (int s=0; s<S; s++){
    int i = 0;
    
    while(times[i] < lms[s] & i < n){
      i += 1;
    }
    
    while (times[i] <= lms[s] + w & i < n){
      ipl = ipl + status[i]*(prognostic(i, s) - log(s0(i, s)));
      i += 1;
    }
  }
  
  return ipl;
}


