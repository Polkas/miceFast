#include <RcppArmadillo.h>
#include "miceFast.h"

// Return a 0/1 indicator vector: 1 if the row of x has no NaN, 0 otherwise
arma::uvec complete_cases_mat(arma::mat &x) {
  unsigned int N = x.n_rows;
  arma::uvec out(N, arma::fill::ones);
  for (unsigned int i = 0; i < N; i++) {
    if (x.row(i).has_nan()) {
      out(i) = 0;
    }
  }
  return out;
}

// Return a 0/1 indicator vector: 1 if y(i) is finite, 0 otherwise
arma::uvec complete_cases_vec(arma::colvec &y) {
  unsigned int N = y.n_elem;
  arma::uvec out(N, arma::fill::ones);
  for (unsigned int i = 0; i < N; i++) {
    if (y.row(i).has_nan()) {
      out(i) = 0;
    }
  }
  return out;
}

// Check that posit_y is not among posit_x
bool different_y_and_x(int posit_y, arma::uvec posit_x) {
  arma::uvec matches = arma::find(posit_x == posit_y);
  return matches.n_elem == 0;
}

// Check that all elements of posit_x are unique
bool different_x(arma::uvec posit_x) {
  arma::uvec uni = arma::unique(posit_x);
  return uni.n_elem == posit_x.n_elem;
}

// Force a matrix to be symmetric: (x + x') / 2
arma::mat sym(arma::mat x) {
  return (x + x.t()) / 2;
}
