#ifndef derivs_rules_H
#define derivs_rules_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::NumericMatrix chainrule_bi (Rcpp::NumericMatrix f, Rcpp::NumericMatrix g, Rcpp::List tri, int deriv_order) ;
  
Rcpp::NumericMatrix chainrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix sumrule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix sumrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix differencerule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix differencerule (Rcpp::List f_list, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix productrule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix productrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix quotientrule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order) ;

Rcpp::NumericMatrix quotientrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order) ;


#endif
