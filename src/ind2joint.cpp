#include <RcppArmadillo.h>
#include "trind.h"
#include "derivs_rules.h"


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix ind2joint_bi (Rcpp::NumericMatrix f, Rcpp:: NumericMatrix g, Rcpp::List tri_f, Rcpp::List tri_g, Rcpp::List tri_h, int deriv_order){
  
  arma::mat f0 = as<arma::mat>(f) ;
  int n = f0.n_rows ; //Number of observations
  arma::mat f1(n, max(vectorise(as<arma::mat>(tri_f[0]))) ) ; 
  arma::mat f2(n, max(vectorise(as<arma::mat>(tri_f[1]))) ) ; 
  arma::mat f3(n, max(vectorise(as<arma::mat>(tri_f[2]))) ) ;
  arma::mat f4(n, max(vectorise(as<arma::mat>(tri_f[3]))) ) ;
  
  arma::mat g0 = as<arma::mat>(g) ;
  arma::mat g1(n, max(vectorise(as<arma::mat>(tri_g[0]))) ) ; 
  arma::mat g2(n, max(vectorise(as<arma::mat>(tri_g[1]))) ) ; 
  arma::mat g3(n, max(vectorise(as<arma::mat>(tri_g[2]))) ) ;
  arma::mat g4(n, max(vectorise(as<arma::mat>(tri_g[3]))) ) ;
  
  int nparf = f1.n_cols ; //Number of parameters of f
  int nparg = g1.n_cols ; //Number of parameters of g
  int nparh = nparf+nparg ;
  //List tri_h = trind_generator(nparh) ;
  
  arma::mat h0 = join_rows(f0,g0) ;
  arma::mat h1(n, max(vectorise(as<arma::mat>(tri_h[0]))) ) ; 
  arma::mat h2(n, max(vectorise(as<arma::mat>(tri_h[1]))) ) ; 
  arma::mat h3(n, max(vectorise(as<arma::mat>(tri_h[2]))) ) ;
  arma::mat h4(n, max(vectorise(as<arma::mat>(tri_h[3]))) ) ;
  
  if(deriv_order>0){
    g1 += as<arma::mat>(g.attr("d1")) ;
    g2 += as<arma::mat>(g.attr("d2")) ;
    f1 += as<arma::mat>(f.attr("d1")) ;
    f2 += as<arma::mat>(f.attr("d2")) ;
    if(deriv_order>2){ 
      g3 += as<arma::mat>(g.attr("d3")) ;
      g4 += as<arma::mat>(g.attr("d4")) ;
      f3 += as<arma::mat>(f.attr("d3")) ;
      f4 += as<arma::mat>(f.attr("d4")) ;
    } 
  }
  
  if(deriv_order>0){
    for (int i = 0; i < nparh; i++) {
      //First derivative
      arma::vec val(n) ;
     
      if(i<nparf){
        val += f1.col(trind(tri_f, {(unsigned int)i})) ;
      } 
      
      if(i>=nparf){
        val += g1.col(trind(tri_g, {(unsigned int)(i-nparf)})) ;
      }
      
      h1.col(trind(tri_h, {(unsigned int)i})) += val ;
      
      for (int j = i; j < nparh; j++) {
        //Second derivative
        arma::vec val(n) ;
        
        if((i<nparf) & (j<nparf)){
          val += f2.col(trind(tri_f, {(unsigned int)i,(unsigned int)(j)})) ;
        } 
        
        if((i>=nparf) & (j>=nparf)){
          val += g2.col(trind(tri_g, {(unsigned int)(i-nparf),(unsigned int)(j-nparf)})) ;
        }
        
        h2.col(trind(tri_h, {(unsigned int)i, (unsigned int)j})) += val ;
        
        if(deriv_order>2){
          for (int k = j; k < nparh; k++) {
            //Third derivative
            arma::vec val(n) ;
            
            if((i<nparf) & (j<nparf) & (k<nparf)){
              val += f3.col(trind(tri_f, {(unsigned int)i,(unsigned int)(j),(unsigned int)(k)})) ;
            } 
            
            if((i>=nparf) & (j>=nparf) & (k>=nparf)){
              val += g3.col(trind(tri_g, {(unsigned int)(i-nparf),(unsigned int)(j-nparf),(unsigned int)(k-nparf)})) ;
            }
            
            h3.col(trind(tri_h, {(unsigned int)i, (unsigned int)j, (unsigned int)k})) += val ;
            for (int l = k; l < nparh; l++) {
              //Fourth derivative
              arma::vec val(n) ;
              
              if((i<nparf) & (j<nparf) & (k<nparf) & (l<nparf)){
                val += f4.col(trind(tri_f, {(unsigned int)i,(unsigned int)(j),(unsigned int)(k),(unsigned int)(l)})) ;
              } 
              
              if((i>=nparf) & (j>=nparf) & (k>=nparf) & (l>=nparf)){
                val += g4.col(trind(tri_g, {(unsigned int)(i-nparf),(unsigned int)(j-nparf),(unsigned int)(k-nparf),(unsigned int)(l-nparf)})) ;
              }
              
              h4.col(trind(tri_h, {(unsigned int)i, (unsigned int)j, (unsigned int)k, (unsigned int)l})) += val ;
            }  
          }
        }
      }
    }
  }
  
  NumericMatrix out = list2derivs(List::create(h0, h1, h2, h3, h4), 4) ;
  
  
  return out ; 
}

//' Independent to joint function
//'
//' Combines multiple derivs objects into a single derivs object.
//' 
//' @return Returns a derivs object.
//'
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = (f_1(x_{n1}), f_2(x_{n2}), ... ,f_M(x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the independent functions are combined.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @param tri_f_list list of length \eqn{K} trind_generator objects, the \eqn{kth} element corresponds to \eqn{kth} derivs object.
//' @param tri_h_list list of length \eqn{K} trind_generator objects, the \eqn{kth} element corresponds to a derivs object with \eqn{k \cdot (k+1)/2} parameters.
//' @inheritParams chainrule
//' 
//' @examples
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' B_derivs<-transform(A, type="exp", par=0, deriv_order=4)
//' ind2joint (list(A_derivs,B_derivs),
//'            list(trind_generator(1),trind_generator(1)),
//'            list(trind_generator(1),trind_generator(1+1)), 4)
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ind2joint (Rcpp::List f_list, Rcpp::List tri_f_list, Rcpp::List tri_h_list, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[0]; //first element of list
  
  for (int i = 0; i < n-1; i++) {
    out = ind2joint_bi (out, f_list[i+1], tri_h_list[i], tri_f_list[i+1], tri_h_list[i+1], deriv_order) ;
  }
  
  return out ;
}
