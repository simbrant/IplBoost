/*
  This file contains a simple function that scales the columns of a
  matrix by individual weights for each column.
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void scale_columns(NumericMatrix mat, NumericVector vec, int m, int n){
  /*
    mat is a (m x n) matrix
    vec is a n-dimensional vector
    m, and n are the number of rows and columns of mat
  */
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      mat(i, j) = mat(i, j)/vec[j];
    }
  }
}

