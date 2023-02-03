#ifndef trind_H
#define trind_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List trind_generator (int K) ;

int trind (Rcpp::List tri, arma::uvec i_deriv) ;

Rcpp::NumericMatrix  list2derivs (Rcpp::List f, int deriv_order) ;

#endif
