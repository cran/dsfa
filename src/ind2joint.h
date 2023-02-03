#ifndef ind2joint_H
#define ind2joint_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::NumericMatrix ind2joint_bi (Rcpp::NumericMatrix f, Rcpp:: NumericMatrix g, Rcpp::List tri_f, Rcpp::List tri_g, Rcpp::List tri_h, int deriv_order) ;

Rcpp::NumericMatrix ind2joint (Rcpp::List f_list, Rcpp::List tri_f_list, Rcpp::List tri_h_list, int deriv_order) ;

#endif
