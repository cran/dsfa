#include <RcppArmadillo.h>
#include <boost/math/special_functions/owens_t.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat erfcinv (arma::mat x){
  //int ncol = x.n_cols ;
  //int n = x.n_rows ;
  arma::vec x2_vec = vectorise(x/2) ;
  NumericVector qnorm_x = Rcpp::qnorm(NumericVector(x2_vec.begin(),x2_vec.end())) ; 
  arma::vec out_vec = -as<arma::vec>(qnorm_x)/sqrt(2) ;
  arma::mat out = reshape(out_vec, size(x)) ;
  
  return out;
}

/*** R
erfcinv(matrix(c(1:9)/10,ncol=3))
*/

// [[Rcpp::export]]
arma::vec OwenT (arma::vec h, arma::vec a){
  int n = h.size() ;
  arma::vec out (n) ;
  for (int i = 0; i < n; i++) {
    out(i) = boost::math::owens_t(h(i),a(i));
  }
  
  return out;
}
