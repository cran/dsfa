#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' trind function
//'
//' Provides the column index of the required derivative for the specified order of a \code{derivs} object. 
//'
//' @param tri list; created by the function [trind_generator()].
//' @param part_deriv_var integer vector; specifies \eqn{\frac{\partial^J f(\cdot)}{\partial x_{ni_1} ... \partial x_{ni_J}}}.
//' The length of the vector is denoted as \eqn{J} and determines the order of the partial derivatives with maximum four. The element \eqn{i_j \in \{0,...,K-1\}}  
//' specifies the variable with respect to which the derivative is taken, where \eqn{j \in \{1,...,J\}}, The order corresponds to the order of derivatives.
//' For example \code{c(0,0,1,2)} is equal to \eqn{\frac{\partial^4 f(\cdot)}{\partial x_{n1} \partial x_{n1} \partial x_{n2} \partial x_{n3}}}.
//' See details for more information.
//' 
//' @details  Let \eqn{f:\mathbb{R}^K -> \mathbb{R}^L, (x_{n1},x_{n2},...,x_{nK}) -> f(x_{n1},x_{n2},...,x_{nK})} be differentiable up to order four w.r.t all parameters \eqn{x_{nk}}, where \eqn{k \in \{1,...,K\}} and \eqn{n \in \{1,...,N\}}.
//' Then a \code{derivs} class object is a numeric matrix with \eqn{N} rows and \eqn{L} columns. \eqn{N} is the length  of the input vectors. Further, it has the following attributes:
//' \enumerate{
//' \item  `d1`: a numeric matrix of the first derivatives w.r.t all parameters,
//'  where the \eqn{nth} row  corresponds to: \eqn{(\frac{\partial   f(\cdot)}{\partial x_{n1}}, \frac{\partial f(\cdot)}{\partial x_{n1}},...,\frac{\partial f(\cdot)}{\partial x_{nK}})}
//' \item  `d2`: a numeric matrix of the second derivatives w.r.t all parameters,
//'  where the \eqn{nth} row  corresponds to: \eqn{(\frac{\partial^2 f(\cdot)}{\partial x_{n1} \partial x_{n1}}, \frac{\partial^2 f(\cdot)}{\partial x_{n1} \partial x_{n2}},...,\frac{\partial^2 f(\cdot)}{\partial x_{nK} \partial x_{nK}})}
//' \item  `d3`: a numeric matrix of the third derivatives w.r.t all parameters,
//'  where the \eqn{nth} row  corresponds to: \eqn{(\frac{\partial^3 f(\cdot)}{\partial x_{n1} \partial x_{n1} \partial x_{n1}}, \frac{\partial^3 f(\cdot)}{\partial x_{n1} \partial x_{n1} \partial x_{n2}},...,\frac{\partial^3 f(\cdot)}{\partial x_{nK} \partial x_{nK} \partial x_{nK}})}
//' \item  `d4`: a numeric matrix of the fourth derivatives w.r.t all parameters,
//'  where the \eqn{nth} row  corresponds to: \eqn{(\frac{\partial^4 f(\cdot)}{\partial x_{n1} \partial x_{n1} \partial x_{n1} \partial x_{n1}}, \frac{\partial^4 f(\cdot)}{\partial x_{n1} \partial x_{n1} \partial x_{n1} \partial x_{n2}},...,\frac{\partial^4 f(\cdot)}{\partial x_{nK} \partial x_{nK} \partial x_{nK} \partial x_{nK}})}
//' }
//' The function \code{trind()} provides the index for the corresponding derivatives. The \code{derivs} class object allows for a modular system which can be easily extended and is faster than numerical derivatives.
//' The advantage compared to analytical derivatives provided by 'mathematica' or \code{\link[stats:deriv]{deriv()}} is that asymptotics and approximations can be used for individual parts.
//' Handwritten derivatives can be tedious at times and may be prone to errors. Thus, the \code{derivs} class object can be used by lazy users.
//' Mainly intended for internal use.
//' 
//' @return Integer, the index for a derivs object.
//'
//' @examples
//' tri=trind_generator(3)
//' trind(tri, c(2,1))
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
int trind (Rcpp::List tri, arma::uvec part_deriv_var){
  int deriv_order = part_deriv_var.size() ;
  
  //Initialize out
  int out = 0;
  
  if(deriv_order==1){
    arma::umat i1 = tri["i1"] ;
    out += i1(part_deriv_var(0), 0) ;
  }
  
  if(deriv_order==2){
    arma::umat i2 = tri["i2"] ;
    out += i2(part_deriv_var(0), part_deriv_var(1)) ;
  }
  
  if(deriv_order==3){
    arma::umat i3 = tri["i3"] ;
    uvec i_start = tri["i_start"] ;
    uvec i_end = tri["i_end"] ;
    out += i3.submat(span(i_start(part_deriv_var(0)), i_end(part_deriv_var(0))), span::all)(part_deriv_var(1),part_deriv_var(2)) ;
  }
  
  if(deriv_order==4){
    arma::umat i4 = tri["i4"] ;
    uvec i_start = tri["i_start"] ;
    uvec i_end = tri["i_end"] ;
    out += i4.submat(span(i_start(part_deriv_var(0)), i_end(part_deriv_var(0))), span(i_start(part_deriv_var(1)), i_end(part_deriv_var(1))))(part_deriv_var(2),part_deriv_var(3)) ;
  }
  
  //Correct for C++ indexing start at 0
  out+= -1 ;
  return out ;
}

//' Trind_generator function
//'
//' Generates index matrices for upper triangular storage up to order four. 
//'
//' @param K integer; determines the number of parameters.
//'
//' @details Useful when working with higher order derivatives, which generate symmetric arrays. Mainly intended for internal use. Similar to 'mgcv::trind.generator'. Mostly internal function.
//' 
//' @return Returns a list with index matrices for the first to fourth derivative, which can be accessed via the function [trind()].
//' The numerical vectors \code{i_start} and \code{i_end} hold the starting and ending indexes, which are required by [trind()] for derivatives greater than two.
//'
//' @examples
//' tri<-trind_generator(3)
//' tri_mgcv<-mgcv::trind.generator(3)
//' 
//' for(i in 1:3){
//'   print(i==trind(tri, part_deriv_var=c(i)-1)+1)
//'   for(j in i:3){
//'     print(tri_mgcv$i2[i,j]==trind(tri, part_deriv_var=c(i,j)-1)+1)
//'     for(k in j:3){
//'       print(tri_mgcv$i3[i,j,k]==trind(tri, part_deriv_var=c(i,j,k)-1)+1)
//'       for(l in k:3){
//'         print(tri_mgcv$i4[i,j,k,l]==trind(tri, part_deriv_var=c(i,j,k,l)-1)+1)
//'       } 
//'     } 
//'   }
//'}
//'
//' @family derivs
//'
//' @export
// [[Rcpp::export]]
Rcpp::List trind_generator (int K){
  //Initialize indicator integer matrices and i_start and i_end integer vectors
  arma::umat i1(K,1) ;
  arma::umat i2(K,K) ;
  arma::umat i3(K*K,K) ;
  arma::umat i4(K*K,K*K) ;
  
  arma::uvec i_start(K) ;
  arma::uvec i_end(K) ;
  
  //Fill i_start and i_end with starting and ending indexes
  int i_start_index = 0 ;
  
  for (int i = 0; i < K; i++) {
    i_start(i) = i_start_index ;
    i_start_index +=K ;
    i_end(i) = i_start_index-1 ;
  }
  
  //Initiliaze indexes
  int m1 = 1;
  int m2 = 1;
  int m3 = 1;
  int m4 = 1;
  
  //Fill i1, i2, i3, i4 with indexes
  for (int i = 0; i < K; i++) {
    i1(i,0)=m1 ;
    m1+=1 ;
    for (int j = i; j < K; j++) {
      i2(i,j)=m2 ;
      i2(j,i)=m2 ;
      m2+=1 ;
      for (int k = j; k < K; k++) {
        i3.submat(span(i_start(i), i_end(i)), span::all)(j,k)= m3 ;
        i3.submat(span(i_start(i), i_end(i)), span::all)(k,j)= m3 ;
        
        i3.submat(span(i_start(j), i_end(j)), span::all)(i,k)= m3 ;
        i3.submat(span(i_start(j), i_end(j)), span::all)(k,i)= m3 ;
        
        i3.submat(span(i_start(k), i_end(k)), span::all)(i,j)= m3 ;
        i3.submat(span(i_start(k), i_end(k)), span::all)(j,i)= m3 ;
        m3+=1 ;
        for (int l = k; l < K; l++) {
          i4.submat(span(i_start(i), i_end(i)), span(i_start(j), i_end(j)))(k,l)= m4 ;
          i4.submat(span(i_start(i), i_end(i)), span(i_start(j), i_end(j)))(l,k)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(i), i_end(i)))(k,l)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(i), i_end(i)))(l,k)= m4 ;
          
          i4.submat(span(i_start(i), i_end(i)), span(i_start(k), i_end(k)))(j,l)= m4 ;
          i4.submat(span(i_start(i), i_end(i)), span(i_start(k), i_end(k)))(l,j)= m4 ;
          i4.submat(span(i_start(k), i_end(k)), span(i_start(i), i_end(i)))(j,l)= m4 ;
          i4.submat(span(i_start(k), i_end(k)), span(i_start(i), i_end(i)))(l,j)= m4 ;
          
          i4.submat(span(i_start(i), i_end(i)), span(i_start(l), i_end(l)))(j,k)= m4 ;
          i4.submat(span(i_start(i), i_end(i)), span(i_start(l), i_end(l)))(k,j)= m4 ;
          i4.submat(span(i_start(l), i_end(l)), span(i_start(i), i_end(i)))(j,k)= m4 ;
          i4.submat(span(i_start(l), i_end(l)), span(i_start(i), i_end(i)))(l,k)= m4 ;
          
          i4.submat(span(i_start(k), i_end(k)), span(i_start(j), i_end(j)))(i,l)= m4 ;
          i4.submat(span(i_start(k), i_end(k)), span(i_start(j), i_end(j)))(l,i)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(k), i_end(k)))(i,l)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(k), i_end(k)))(l,i)= m4 ;
          
          i4.submat(span(i_start(l), i_end(l)), span(i_start(j), i_end(j)))(i,k)= m4 ;
          i4.submat(span(i_start(l), i_end(l)), span(i_start(j), i_end(j)))(k,i)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(l), i_end(l)))(i,k)= m4 ;
          i4.submat(span(i_start(j), i_end(j)), span(i_start(l), i_end(l)))(k,i)= m4 ;
          
          i4.submat(span(i_start(k), i_end(k)), span(i_start(l), i_end(l)))(i,j)= m4 ;
          i4.submat(span(i_start(k), i_end(k)), span(i_start(l), i_end(l)))(j,i)= m4 ;
          i4.submat(span(i_start(l), i_end(l)), span(i_start(k), i_end(k)))(i,j)= m4 ;
          i4.submat(span(i_start(l), i_end(l)), span(i_start(k), i_end(k)))(j,i)= m4 ;
          m4+=1 ;
        }
      }
    }
  }
  
  //Initialize list output 
  Rcpp::List out = Rcpp::List::create(Named("i1")=i1, Named("i2")=i2, Named("i3")=i3, Named("i4")=i4, Named("i_start")=i_start, Named("i_end")=i_end) ;
  
  return out ;
}

//' list2derivs
//'
//' Transforms a list of matrices d0, d1, d2, d3, d4 to a \code{derivs} object.
//'
//' @param f list of matrices; d0, d1, d2, d3, d4 
//' @param deriv_order integer; maximum order of derivative. Available are \code{0},\code{2} and \code{4}.
//' 
//' @return Mostly internal function. Returns an object of class \code{derivs}
//' For more details see [trind()] and [trind_generator()].
//' 
//' @examples
//' A<-matrix(c(1:9)/10, ncol=3)
//' list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix  list2derivs (Rcpp::List f, int deriv_order){
  //Transforms a list of matrices to derivs object
  Rcpp::NumericMatrix out = f[0] ;
  
  if(deriv_order>0){
    out.attr("d1") = f[1] ;
    out.attr("d2") = f[2] ;
    if(deriv_order>2){
      out.attr("d3") = f[3] ;
      out.attr("d4") = f[4] ;
    }
  }
  
  return out ;
}
