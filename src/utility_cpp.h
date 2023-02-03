#ifndef utility_cpp_H
#define utility_cpp_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat erfcinv(arma::mat x) ;

arma::vec OwenT (arma::vec h, arma::vec a) ;

#endif
