#ifndef copulas_H
#define copulas_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::NumericMatrix  dcop_unrot_cpp (arma::vec u, arma::vec v, arma::vec p, Rcpp::String distr_cop, int deriv_order, Rcpp::List tri, bool logp) ;

Rcpp::NumericMatrix  dcop_cpp (arma::vec u, arma::vec v, arma::vec p, Rcpp::String distr_cop, int rot, int deriv_order, Rcpp::List tri, bool logp) ;

#endif
