#ifndef transform_H
#define transform_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::NumericMatrix  transform (arma::mat x, Rcpp::String type, arma::vec par, int deriv_order) ;

Rcpp::NumericMatrix  derivs_transform (Rcpp::NumericMatrix f, Rcpp::String type, arma::vec par, Rcpp::List tri, int deriv_order) ;

#endif
