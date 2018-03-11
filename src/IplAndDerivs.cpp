/*
 This file contains internal C++ functions to compute the first and second derivatives of van Houwelingens
 ipl for one covariate index, for all landmarks, as well as evaluate the ipl itself.
*/

#include <Rcpp.h>
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
  NumericVector u(S);
  
  int lmsind = 0;
  for (int s=0; s<S; s++){
    // Assume times are ordered
    int i = lmsind;
    while (times[i] < lms[s]){
      i += 1;
    }
    lmsind = i;
    while (times[i] <= lms[s] + w){
      u[s] += status[i]*(mat(i, j-1) - S1j(i, s)/S0(i, s));
      i += 1;
    } 
  }   
  return u;
}


// [[Rcpp::export]]
NumericVector compute_minI_j(int j, NumericVector status, NumericVector times, NumericMatrix S0,
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
  
  NumericVector I(S);
  
  int lmsind = 0;
  for (int s=0; s<S; s++){
    // Assume times are ordered
    int i = lmsind;
    while (times[i] < lms[s]){
      i += 1;
    }
    lmsind = i;
    while (times[i] <= lms[s] + w){
      I[s] += status[i]*((S2j(i, s)*S0(i, s)- S1j(i, s)*S1j(i, s))/(S0(i, s)*S0(i, s)));
      i += 1;
    } 
  }
  for (int s = 0; s<S; s++){
    I[s] += lambda[s];
  }
  return I;
}


//[[Rcpp::export]]
double compute_ipl(NumericVector times, NumericVector status, NumericMatrix mat, NumericMatrix betas,
                   NumericVector lms, double w, int S, int n, int p) {
  double ipl = 0;
  
  int lms_first_index = 0;
  for (int s=0; s<S; s++){
    // compute linear predictor / prognostic index
    NumericVector pi_s(n);
    NumericVector risk_s(n);
    
    //  for (int i=0; i<n; i++){
    //  pi_s[i]  = sum(mat(i, _)*betas(s, _));
    //}
    for (int i=0; i<n; i++){
      double this_pi = 0;
      for (int j=0; j<p; j++){
        this_pi += mat(i, j)*betas(s, j);
      }
      pi_s[i] = this_pi;
    }
    
    risk_s = exp(pi_s); // Rcpp sugar
    // Backwards cumsum of risk
    for (int i=n-2; i>=0; i--){
      risk_s[i] = risk_s[i] + risk_s[i+1];
    }
    NumericVector logcumrisksum = log(risk_s);
    
    int i = lms_first_index;
    while (times[i] < lms[s]){
      i += 1;
    }
    lms_first_index = i;
    while (times[i] <= lms[s] + w){
      ipl += status[i]*(pi_s[i] - logcumrisksum[i]);
      i += 1;
    } 
  }
  return ipl;
}


