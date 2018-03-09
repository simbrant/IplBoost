/*
  This file contains internal C++ functions to compute S0, S1, S2
  */

# include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix compute_S0(NumericMatrix risk, NumericVector times, int n,
                         int S) {
  /*
  Computes S0 = sum_{l in R_i} exp(beta(LM_s)^Tx_l), i = 1,..n, s = 1,..S
      
      Args:
       - risk is a (n x S) matrix (exp(...))
       - times is an n-dimensional vector
       - n i the length of times (number of observations)
       - S is the number of rows of risk (number of landmarks)
   
   returns a (n x S) matrix
  */
    
    NumericMatrix S0 = clone(risk);
    for (int s=0; s<S; s++){
      NumericVector tmp = S0(_, s);
      IntegerVector idx = seq_along(tmp) - 1;
      std::sort(idx.begin(), idx.end(),[&](int i, int j){return times[i] < times[i];});
      
      tmp = tmp[idx];
      for (int i=0; i<n; i++){
        S0(i, s) = tmp[i];
      }
      
      for (int i=n-2; i >= 0; i--){
        S0(i, s) = S0(i, s) + S0(i+1, s);
      }
    }
    return S0;
}

// [[Rcpp::export]]
NumericMatrix compute_S1_j(int j, NumericMatrix risk, NumericVector times, NumericMatrix mat,
                           int n, int S) {
  /*
   Computes S1.j = sum_{l in R_i} x_{lj}exp(beta(LM_s)^Tx_l), i = 1,..n, s = 1,..S
  
  Args:
    - j is the covariate index to compute the matrix for
    - risk is a (n x S) matrix (exp(...))
    - times is an n-dimensional vector
    - n i the length of times (number of observations)
    - S is the number of rows of risk (number of landmarks)
    - mat is the design matrix (n x p)
   
   returns a (n x S) matrix
  */
  
    NumericMatrix S1 = clone(risk);
    NumericVector cov_vec = mat(_, j-1);
    
    for (int s=0; s<S; s++){
      NumericVector tmp = S1(_, s);
      IntegerVector idx = seq_along(tmp) - 1;
      std::sort(idx.begin(), idx.end(), [&](int i, int j){return times[i] < times[i];});
      
      tmp = tmp[idx];
      cov_vec = cov_vec[idx];
      for (int i=0; i<n; i++){
        S1(i, s) = cov_vec[i]*tmp[i];
      }
      
      // cumsum   
      for (int i=n-2; i >= 0; i--){
        S1(i, s) = S1(i, s) + S1(i+1, s);
      }
    }
    return S1;
}


// [[Rcpp::export]]
NumericMatrix compute_S2_j(int j, NumericMatrix risk, NumericVector times, NumericMatrix mat, int n,
                           int S) {
  /*
   Computes S2.j = sum_{l in R_i} x_{lj}**2exp(beta(LM_s)^Tx_l), i = 1,..n, s = 1,..S
   
   Args:
   - j is the covariate index to compute the matrix for
   - risk is a (n x S) matrix (exp(...))
   - times is an n-dimensional vector
   - n i the length of times (number of observations)
   - S is the number of rows of risk (number of landmarks)
   - mat is the design matrix (n x p)
   
   returns a (n x S) matrix
   */
    NumericMatrix S2 = clone(risk);
    NumericVector cov_vec = mat(_, j-1);
    
    for (int s=0; s<S; s++){
      NumericVector tmp = S2(_, s);
      IntegerVector idx = seq_along(tmp) - 1;
      std::sort(idx.begin(), idx.end(), [&](int i, int j){return times[i] < times[i];});
      
      tmp = tmp[idx];
      cov_vec = cov_vec[idx];
      for (int i=0; i<n; i++){
        S2(i, s) = cov_vec[i]*cov_vec[i]*tmp[i];
      }
      
      // backwards cumsum 
      for (int i=n-2; i >= 0; i--){
        S2(i, s) = S2(i, s) + S2(i+1, s);
      }
    }
    return S2;
}