#include <RcppArmadillo.h>
#include "trind.h"
#include "derivs_rules.h"
#include "utility_cpp.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' transform
//'
//' Transforms a matrix via the specified function.
//' 
//' @details Takes the numeric matrix x as an input for the function specified by \code{type} and evaluates it together with the derivatives.
//'  
//' @return Returns an object of class \code{derivs}.
//'
//' @examples
//' A<-matrix(c(1:9)/10, ncol=3)
//' A_mat<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' transform(x=transform(x = A, type="exp", par=0, deriv_order=4), type="log", deriv_order=4, par = 0)
//'
//' @inheritParams trind
//' @inheritParams list2derivs
//' 
//' @param x numeric matrix to be transformed.
//' @param type string, specifies the transformation function. Available are:
//' \enumerate{
//' \item  `identity`: \eqn{f(x)=x}.
//' \item  `exp`: \eqn{f(x)=\exp\{x\}}.
//' \item  `log`: \eqn{f(x)=\log\{x\}}.
//' \item  `glogit`: \eqn{f(x)=\log\{(-x + min)/(x - max)}, where \code{par=c(min, max)}.
//' \item  `glogitinv`: \eqn{f(x)=\exp\{x\} \cdot (max + min)/(1 + \exp\{x\}) }, where \code{par=c(min, max)}.
//' \item  `inv`: \eqn{f(x)=\frac{1}{x}}.
//' \item  `pnorm`: \eqn{f(x)=\Phi(x)}.
//' \item  `qnorm`: \eqn{f(x)=\Phi^{-1}(x)}.
//' \item  `mexp`: \eqn{f(x)=-\exp\{x\}}.
//' \item  `zeta`: \eqn{f(x)=\log\{2 \cdot \Phi(x)\}}.
//' \item  `constant`: \eqn{f(x)=c}.
//' \item  `chainrule_utility`: \eqn{f(x)=f'(x)=f''(x)=f'''(x)=f''''(x)}.
//' }
//' @param par numeric vector, additional parameters, e.g. min and max for \code{glogit}.
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix  transform (arma::mat x, Rcpp::String type, arma::vec par, int deriv_order){
  
  arma::mat ret(x.n_rows, x.n_cols ) ;
  arma::mat d1(x.n_rows, x.n_cols ) ;
  arma::mat d2(x.n_rows, x.n_cols ) ;
  arma::mat d3(x.n_rows, x.n_cols ) ;
  arma::mat d4(x.n_rows, x.n_cols ) ;
  
  Rcpp::StringVector types = {"identity","exp","log","glogit","glogitinv","inv","pnorm","qnorm","mexp","zeta","constant","chainrule_utility"} ;
  
  //identity function
  if(type==types(0)){ 
    ret += x ;
    if(deriv_order>0){
      d1 += 1 ;
      d2 += 0 ;
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //exp function
  if(type==types(1)){ 
    ret += exp(x) ;
    if(deriv_order>0){
      d1 += ret ;
      d2 += ret ;
      if(deriv_order>2){
        d3 += ret ;
        d4 += ret ;
      }
    }
  }
  
  //log function
  if(type==types(2)){ 
    ret += log(x) ;
    if(deriv_order>0){
      d1 += 1/x ;
      d2 += -pow(x,-2) ;
      if(deriv_order>2){
        d3 += 2/pow(x,3) ;
        d4 += 6/pow(x,4) ;
      }
    }
  }
  
  //glogit function
  if(type==types(3)){
    int min = par(0) ;
    int max = par(1) ;
    
    ret += log((-x + min)/(x - max)) ;
    if(deriv_order>0){
      d1 += (-max + min)/((x - max)%(x - min)) ;
      d2 += pow(x - max,-2) - pow(x - min,-2) ;
      if(deriv_order>2){
        d3 += 2*(pow(-x + max,-3) + pow(x - min,-3)) ;
        d4 += 6/pow(x - max,4) - 6/pow(x - min,4) ;
      }
    }
  }
  
  //glogitinv function
  if(type==types(4)){
    int min = par(0) ;
    int max = par(1) ;
    
    ret += (exp(x)*max + min)/(1 + exp(x)) ;
    if(deriv_order>0){
      d1 += (exp(x)*(max - min))/pow(1 + exp(x),2) ;
      d2 += -((exp(x)%(-1 + exp(x))*(max - min))/pow(1 + exp(x),3)) ;
      if(deriv_order>2){
        d3 += (exp(x)%(1 - 4*exp(x) + exp(2*x))*(max - min))/pow(1 + exp(x),4) ;
        d4 += -((exp(x)%(-1 + 11*exp(x) - 11*exp(2*x) + exp(3*x))*(max - min))/pow(1 + exp(x),5)) ;
      }
    }
  }
  
  //inverse function
  if(type==types(5)){ 
    ret += 1/x ;
    if(deriv_order>0){
      d1 += -pow(x,-2) ;
      d2 += 2/pow(x,3) ;
      if(deriv_order>2){
        d3 += -6/pow(x,4) ;
        d4 += 24/pow(x,5) ;
      }
    }
  }
  
  //pnorm function
  if(type==types(6)){ 
    ret += (1 + erf(x/sqrt(2)))/2. ;
    if(deriv_order>0){
      d1 += 1/(exp(pow(x,2)/2.)*sqrt(2*datum::pi)) ;
      d2 += -(x/(exp(pow(x,2)/2.)*sqrt(2*datum::pi))) ;
      if(deriv_order>2){
        d3 += (-1 + pow(x,2))/(exp(pow(x,2)/2.)*sqrt(2*datum::pi)) ;
        d4 += -((x%(-3 + pow(x,2)))/(exp(pow(x,2)/2.)*sqrt(2*datum::pi))) ;
      }
    }
  }
  
  //qnorm function
  if(type==types(7)){ 
    ret += -(sqrt(2)*erfcinv(2*x)) ;
    if(deriv_order>0){
      d1 += exp(pow(erfcinv(2*x),2))*sqrt(2*datum::pi) ;
      d2 += -2*sqrt(2)*exp(2*pow(erfcinv(2*x),2))*datum::pi%erfcinv(2*x) ;
      if(deriv_order>2){
        d3 += 2*sqrt(2)*exp(3*pow(erfcinv(2*x),2))*pow(datum::pi,1.5)%(1 + 4*pow(erfcinv(2*x),2)) ;
        d4 += -4*sqrt(2)*exp(4*pow(erfcinv(2*x),2))*pow(datum::pi,2)%erfcinv(2*x)%(7 + 12*pow(erfcinv(2*x),2)) ;
      }
    }
  }
  
  //mexp function
  if(type==types(8)){ 
    ret += -exp(x) ;
    if(deriv_order>0){
      d1 += ret ;
      d2 += ret ;
      if(deriv_order>2){
        d3 += ret ;
        d4 += ret ;
      }
    }
  }
  
  //zeta function
  if(type==types(9)){
    if(x.has_nonfinite()){
      x.clamp(-100000,100000) ;
      x.elem( find_nan(x) ).zeros() ;
    }
    arma::mat x2 = pow(x,2) ;
    
    arma::mat log_normcdf = log(normcdf(x)) ;
    ret += log_normcdf + log(2) ;
    
    if(deriv_order>0){
      d1 += -x/(1 - 1/(x2 + 2) + 1/((x2 + 2) % (x2 + 4)) - 5/((x2 + 2) % (x2 + 4) % (x2 + 6)) + 9/((x2 + 2) % (x2 + 4) % (x2 + 6) % (x2 + 8)) - 129/((x2 + 2) % (x2 + 4) % (x2 + 6) % (x2 + 8) % (x2 + 10))) ;
      arma::umat xg50 = find(x > (-50)) ;
      d1.elem( xg50 ) = exp(log_normpdf(x.elem( xg50 )) - log_normcdf.elem( xg50 ));
      d2 += -d1 % (x + d1) ;
      
      if(deriv_order>2){
        d3 += (- d2 % (x + d1) - d1 % (1 + d2)) ;
        d4 += (-d3 % (x + 2 * d1) - 2 * d2 % (1 + d2)) ;
      }
    }
  }
  
  //Constant function
  if(type==types(10)){ 
    ret += x ;
    if(deriv_order>0){
      d1 += 0 ;
      d2 += 0 ;
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Chainrule_utility function
  if(type==types(11)){ 
    ret += x ;
    if(deriv_order>0){
      d1 += x ;
      d2 += x ;
      if(deriv_order>2){
        d3 += x ;
        d4 += x ;
      }
    }
  }
  
  Rcpp::NumericMatrix out = list2derivs(List::create(ret, d1, d2, d3, d4), deriv_order) ;
  return out;
}

//' derivs_transform
//'
//' Transforms a derivs object via the specified function and applies the chainrule.
//' 
//' @details Takes the derivs object \code{f} as an input for the function specified by \code{type} and evaluates it together with the derivatives utilizing the chainrule.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @return Returns an object of class \code{derivs}
//'
//' @examples
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_mat<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' derivs_transform(f =derivs_transform(f = A, type="exp", par=0,
//'                                      tri=trind_generator(1), deriv_order=4),
//'                    type="log", par=0, tri=trind_generator(1), deriv_order=4)
//' 
//' @inheritParams transform
//' @inheritParams trind
//' @inheritParams list2derivs
//' @param f derivs object.
//' 
//' @family derivs
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix  derivs_transform (Rcpp::NumericMatrix f, Rcpp::String type, arma::vec par, Rcpp::List tri, int deriv_order){
  Rcpp::NumericMatrix out = chainrule (List::create(transform( as<arma::mat>(f), type,  par, deriv_order), f), tri, deriv_order) ; 
  
  return out;
}

