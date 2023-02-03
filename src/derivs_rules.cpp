#include <RcppArmadillo.h>
#include "trind.h"
#include "transform.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix chainrule_bi (Rcpp::NumericMatrix f, Rcpp::NumericMatrix g, Rcpp::List tri, int deriv_order){
  
  arma::mat g0 = as<arma::mat>(g) ;
  int n = g0.n_rows ; //Number of observations
  
  arma::mat g1(n, max(vectorise(as<arma::mat>(tri[0]))) ) ; 
  arma::mat g2(n, max(vectorise(as<arma::mat>(tri[1]))) ) ;
  arma::mat g3(n, max(vectorise(as<arma::mat>(tri[2]))) ) ;
  arma::mat g4(n, max(vectorise(as<arma::mat>(tri[3]))) ) ;
  
  arma::mat f0 = as<arma::mat>(f) ;
  arma::mat f1(n, g1.n_cols ) ; 
  arma::mat f2(n, g2.n_cols ) ; 
  arma::mat f3(n, g3.n_cols ) ;
  arma::mat f4(n, g4.n_cols ) ;
  
  arma::mat h0 = f0 ;
  arma::mat h1(n, g1.n_cols ) ; 
  arma::mat h2(n, g2.n_cols ) ; 
  arma::mat h3(n, g3.n_cols ) ;
  arma::mat h4(n, g4.n_cols ) ;
  
  if(deriv_order>0){
    g1 += as<arma::mat>(g.attr("d1")) ;
    g2 += as<arma::mat>(g.attr("d2")) ;
    f1 += repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
    f2 += repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
    
    if(deriv_order>2){ 
      g3 += as<arma::mat>(g.attr("d3")) ;
      g4 += as<arma::mat>(g.attr("d4")) ;
      f3 += repmat(as<arma::mat>(f.attr("d3")), 1, g3.n_cols/as<arma::mat>(f.attr("d3")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;
      f4 += repmat(as<arma::mat>(f.attr("d4")), 1, g4.n_cols/as<arma::mat>(f.attr("d4")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;f3.reshape( n, g3.n_cols ) ;
    } 
  }
  
  if(deriv_order>0){
    int nparg = g1.n_cols ; //Number of parameters of g
    
    arma::uvec i_deriv = uvec(4) ;
    
    for (int i = 0; i < nparg; i++) {
      //First derivative
      int i_index = trind(tri, {(unsigned int)i}) ;
      
      h1.col(i_index) += f1.col(i_index) % g1.col(i_index) ;
      
      for (int j = i; j < nparg; j++) {
        //Second derivative
        int j_index = trind(tri, {(unsigned int)j}) ;
        
        int ij_index = trind(tri, {(unsigned int)i, (unsigned int)j}) ;
        
        h2.col(ij_index) += f2.col(ij_index) % g1.col(j_index) % g1.col(i_index) + f1.col(i_index) % g2.col(ij_index) ;
        
        if(deriv_order>2){
          for (int k = j; k < nparg; k++) {
            //Third derivative
            int k_index = trind(tri, {(unsigned int)k}) ;
            
            int ik_index = trind(tri, {(unsigned int)i, (unsigned int)k}) ;
            int jk_index = trind(tri, {(unsigned int)j, (unsigned int)k}) ;
            
            int ijk_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)k}) ;
            
            h3.col(ijk_index) += f3.col(ijk_index) % g1.col(k_index) % g1.col(j_index) % g1.col(i_index)+
              f2.col(ij_index) % g2.col(jk_index) % g1.col(i_index)+
              f2.col(ij_index) % g1.col(j_index) % g2.col(ik_index)+
              f2.col(ik_index) % g1.col(k_index) % g2.col(ij_index)+
              f1.col(i_index) % g3.col(ijk_index) ;
            
            for (int l = k; l < nparg; l++) {
              //Fourth derivative
              int l_index = trind(tri, {(unsigned int)l}) ;
              
              int il_index = trind(tri, {(unsigned int)i, (unsigned int)l}) ;
              int jl_index = trind(tri, {(unsigned int)j, (unsigned int)l}) ;
              int kl_index = trind(tri, {(unsigned int)k, (unsigned int)l}) ;
              
              int ijl_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)l}) ;
              int ikl_index = trind(tri, {(unsigned int)i, (unsigned int)k, (unsigned int)l}) ;
              int jkl_index = trind(tri, {(unsigned int)j, (unsigned int)k, (unsigned int)l}) ;
              
              int ijkl_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)k, (unsigned int)l}) ;
              
              h4.col(ijkl_index) += f4.col(ijkl_index)%g1.col(l_index)%g1.col(k_index)%g1.col(j_index)%g1.col(i_index)+
                f3.col(ijk_index)%g2.col(kl_index)%g1.col(j_index)%g1.col(i_index)+
                f3.col(ijk_index)%g1.col(k_index)%g2.col(jl_index)%g1.col(i_index)+
                f3.col(ijk_index)%g1.col(k_index)%g1.col(j_index)%g2.col(il_index)+
                //#next line
                f3.col(ijl_index)%g1.col(l_index)%g2.col(jk_index)%g1.col(i_index)+
                f2.col(ij_index)%g3.col(jkl_index)%g1.col(i_index)+
                f2.col(ij_index)%g2.col(jk_index)%g2.col(il_index)+
                //#next line
                f3.col(ijl_index)%g1.col(l_index)%g1.col(j_index)%g2.col(ik_index)+
                f2.col(ij_index)%g2.col(jl_index)%g2.col(ik_index)+
                f2.col(ij_index)%g1.col(j_index)%g3.col(ikl_index)+
                //#next line
                f3.col(ikl_index)%g1.col(l_index)%g1.col(k_index)%g2.col(ij_index)+
                f2.col(ik_index)%g2.col(kl_index)%g2.col(ij_index)+
                f2.col(ik_index)%g1.col(k_index)%g3.col(ijl_index)+
                //#next line
                f2.col(il_index)%g1.col(l_index)%g3.col(ijk_index)+
                f1.col(i_index)%g4.col(ijkl_index) ;
            }  
          }
        }
      }
    }
  }
  
  NumericMatrix out = list2derivs(List::create(h0, h1, h2, h3, h4), deriv_order) ;
  
  return out ;
}

//' Chainrule
//'
//' Chainrule for derivs objects.
//' 
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = f_1(\cdot) \circ f_2(\cdot) ... \circ f_M(x_{n1},x_{n2},...,x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the chainrule is applied.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @return Returns an object of class \code{derivs} for the function \eqn{h(\cdot)}.
//'
//' @inheritParams trind
//' @inheritParams list2derivs
//' @param f_list, list of \code{derivs} objects of length \eqn{M}, e.g. \eqn{list(f_1(\cdot), f_2(\cdot),...,f_M(\cdot))}
//' 
//' @examples 
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' B_derivs<-transform(A, type="exp", par=0, deriv_order=4)
//' C_derivs<-transform(B_derivs, type="log", par=0, deriv_order=4)
//' chainrule(list(C_derivs, B_derivs), trind_generator(1), deriv_order=4) #equal to A_derivs
//'
//' @family derivs
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix chainrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[n-1]; //last element of list
  
  for (int i = 0; i < n-1; i++) {
    out = chainrule_bi(f_list[(n-1)-(i+1)], out, tri, deriv_order) ;
  }
  
  return out ;
} 



// [[Rcpp::export]]
Rcpp::NumericMatrix sumrule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order){
  
  arma::mat g0 = as<arma::mat>(g) ;
  int n = g0.n_rows ; //Number of observations
  
  arma::mat g1(n, max(vectorise(as<arma::mat>(tri[0]))) ) ; //g1 = as<arma::mat>(g.attr("d1")) ;
  arma::mat g2(n, max(vectorise(as<arma::mat>(tri[1]))) ) ; //g2 = as<arma::mat>(g.attr("d2")) ;
  arma::mat g3(n, max(vectorise(as<arma::mat>(tri[2]))) ) ;
  arma::mat g4(n, max(vectorise(as<arma::mat>(tri[3]))) ) ;
  
  arma::mat f0 = as<arma::mat>(f) ;
  arma::mat f1(n, g1.n_cols ) ; //f1 = repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
  arma::mat f2(n, g2.n_cols ) ; //f2 = repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
  arma::mat f3(n, g3.n_cols ) ;
  arma::mat f4(n, g4.n_cols ) ;
  
  arma::mat h0(n, g0.n_cols ) ;
  h0 += f0 + g0 ;
  arma::mat h1(n, g1.n_cols ) ; //f1 = repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
  arma::mat h2(n, g2.n_cols ) ; //f2 = repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
  arma::mat h3(n, g3.n_cols ) ;
  arma::mat h4(n, g4.n_cols ) ;
  
  if(deriv_order>0){
    g1 += as<arma::mat>(g.attr("d1")) ;
    g2 += as<arma::mat>(g.attr("d2")) ;
    f1 += repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
    f2 += repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
    h1 += f1 + g1 ;
    h2 += f2 + g2 ;
    if(deriv_order>2){ 
      g3 += as<arma::mat>(g.attr("d3")) ;
      g4 += as<arma::mat>(g.attr("d4")) ;
      f3 += repmat(as<arma::mat>(f.attr("d3")), 1, g3.n_cols/as<arma::mat>(f.attr("d3")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;
      f4 += repmat(as<arma::mat>(f.attr("d4")), 1, g4.n_cols/as<arma::mat>(f.attr("d4")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;f3.reshape( n, g3.n_cols ) ;
      h3 += f3 + g3 ;
      h4 += f4 + g4 ;
    } 
  }
  
  NumericMatrix out = list2derivs(List::create(h0, h1, h2, h3, h4), deriv_order) ;
  
  return out ; 
}

//' Sumrule
//'
//' Sumrule for derivs objects.
//' 
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = f_1(\cdot) + f_2(\cdot) ... + f_M(x_{n1},x_{n2},...,x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the sumrule is applied.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @return Returns an object of class \code{derivs} for the function \eqn{h(\cdot)}.
//' 
//' @inheritParams chainrule
//' 
//' @examples 
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' sumrule(list(A_derivs, A_derivs), trind_generator(1), deriv_order=4) #equal to 2*A_derivs
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sumrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[0]; //last element of list
  
  for (int i = 0; i < n-1; i++) {
    out = sumrule_bi(f_list[i+1], out, tri, deriv_order) ;
  }
  
  return out ;
} 

// [[Rcpp::export]]
Rcpp::NumericMatrix differencerule_bi (Rcpp::NumericMatrix f,Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order){
  
  arma::mat g0 = as<arma::mat>(g) ;
  int n = g0.n_rows ; //Number of observations
  
  arma::mat g1(n, max(vectorise(as<arma::mat>(tri[0]))) ) ; //g1 = as<arma::mat>(g.attr("d1")) ;
  arma::mat g2(n, max(vectorise(as<arma::mat>(tri[1]))) ) ; //g2 = as<arma::mat>(g.attr("d2")) ;
  arma::mat g3(n, max(vectorise(as<arma::mat>(tri[2]))) ) ;
  arma::mat g4(n, max(vectorise(as<arma::mat>(tri[3]))) ) ;
  
  arma::mat f0 = as<arma::mat>(f) ;
  arma::mat f1(n, g1.n_cols ) ; //f1 = repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
  arma::mat f2(n, g2.n_cols ) ; //f2 = repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
  arma::mat f3(n, g3.n_cols ) ;
  arma::mat f4(n, g4.n_cols ) ;
  
  arma::mat h0(n, g0.n_cols ) ;
  h0 += f0 - g0 ;
  arma::mat h1(n, g1.n_cols ) ; //f1 = repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
  arma::mat h2(n, g2.n_cols ) ; //f2 = repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
  arma::mat h3(n, g3.n_cols ) ;
  arma::mat h4(n, g4.n_cols ) ;
  
  if(deriv_order>0){
    g1 += as<arma::mat>(g.attr("d1")) ;
    g2 += as<arma::mat>(g.attr("d2")) ;
    f1 += repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
    f2 += repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
    h1 += f1 - g1 ;
    h2 += f2 - g2 ;
    if(deriv_order>2){ 
      g3 += as<arma::mat>(g.attr("d3")) ;
      g4 += as<arma::mat>(g.attr("d4")) ;
      f3 += repmat(as<arma::mat>(f.attr("d3")), 1, g3.n_cols/as<arma::mat>(f.attr("d3")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;
      f4 += repmat(as<arma::mat>(f.attr("d4")), 1, g4.n_cols/as<arma::mat>(f.attr("d4")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;f3.reshape( n, g3.n_cols ) ;
      h3 += f3 - g3 ;
      h4 += f4 - g4 ;
    } 
  }
  
  NumericMatrix out = list2derivs(List::create(h0, h1, h2, h3, h4), deriv_order) ;
  
  return out ; 
}


//' Differencerule
//'
//' Differencerule for derivs objects.
//' 
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = f_1(\cdot) - f_2(\cdot) ... - f_M(x_{n1},x_{n2},...,x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the difference rule is applied.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @return Returns an object of class \code{derivs} for the function \eqn{h(\cdot)}.
//' 
//' @inheritParams chainrule
//' 
//' @examples 
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=4)
//' differencerule(list(A_derivs, A_derivs), trind_generator(1), deriv_order=4) #equal to 0
//'
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix differencerule (Rcpp::List f_list, Rcpp::List tri, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[0]; //first element of list
  
  for (int i = 0; i < n-1; i++) {
    out = differencerule_bi(f_list[i+1], out, tri, deriv_order) ;
  }
  
  return out ;
} 

// [[Rcpp::export]]
Rcpp::NumericMatrix productrule_bi (Rcpp::NumericMatrix f, Rcpp::NumericMatrix g, Rcpp::List tri, int deriv_order){
  
  arma::mat g0 = as<arma::mat>(g) ;
  int n = g0.n_rows ; //Number of observations
  
  arma::mat g1(n, max(vectorise(as<arma::mat>(tri[0]))) ) ; 
  arma::mat g2(n, max(vectorise(as<arma::mat>(tri[1]))) ) ;
  arma::mat g3(n, max(vectorise(as<arma::mat>(tri[2]))) ) ;
  arma::mat g4(n, max(vectorise(as<arma::mat>(tri[3]))) ) ;
  
  arma::mat f0 = as<arma::mat>(f) ;
  arma::mat f1(n, g1.n_cols ) ; 
  arma::mat f2(n, g2.n_cols ) ; 
  arma::mat f3(n, g3.n_cols ) ;
  arma::mat f4(n, g4.n_cols ) ;
  
  arma::mat h0 = f0 % g0 ;
  arma::mat h1(n, g1.n_cols ) ; 
  arma::mat h2(n, g2.n_cols ) ; 
  arma::mat h3(n, g3.n_cols ) ;
  arma::mat h4(n, g4.n_cols ) ;
  
  if(deriv_order>0){
    g1 += as<arma::mat>(g.attr("d1")) ;
    g2 += as<arma::mat>(g.attr("d2")) ;
    f1 += repmat(as<arma::mat>(f.attr("d1")), 1, g1.n_cols/as<arma::mat>(f.attr("d1")).n_cols ) ;
    f2 += repmat(as<arma::mat>(f.attr("d2")), 1, g2.n_cols/as<arma::mat>(f.attr("d2")).n_cols ) ;
    
    if(deriv_order>2){ 
      g3 += as<arma::mat>(g.attr("d3")) ;
      g4 += as<arma::mat>(g.attr("d4")) ;
      f3 += repmat(as<arma::mat>(f.attr("d3")), 1, g3.n_cols/as<arma::mat>(f.attr("d3")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;
      f4 += repmat(as<arma::mat>(f.attr("d4")), 1, g4.n_cols/as<arma::mat>(f.attr("d4")).n_cols ) ;//f1.reshape( n, g1.n_cols ) ;f3.reshape( n, g3.n_cols ) ;
    } 
  }
  
  if(deriv_order>0){
    int nparg = g1.n_cols ; //Number of parameters of g
    
    for (int i = 0; i < nparg; i++) {
      //First derivative
      int i_index = trind(tri, {(unsigned int)i}) ;
      
      h1.col(i_index) += f1.col(i_index) % g0.col(i_index) + f0.col(i_index) % g1.col(i_index) ;
      
      for (int j = i; j < nparg; j++) {
        //Second derivative
        int j_index = trind(tri, {(unsigned int)j}) ;
        
        int ij_index = trind(tri, {(unsigned int)i, (unsigned int)j}) ;
        
        h2.col(ij_index) += f2.col(ij_index) % g0.col(i_index) + 
                            f1.col(i_index) % g1.col(j_index) + 
                            f1.col(j_index) % g1.col(i_index) +
                            f0.col(i_index) % g2.col(ij_index) ;
        
        if(deriv_order>2){
          for (int k = j; k < nparg; k++) {
            //Third derivative
            int k_index = trind(tri, {(unsigned int)k}) ;
            
            int ik_index = trind(tri, {(unsigned int)i, (unsigned int)k}) ;
            int jk_index = trind(tri, {(unsigned int)j, (unsigned int)k}) ;
            
            int ijk_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)k}) ;
            
            h3.col(ijk_index) += f3.col(ijk_index) % g0.col(i_index)+
                                 f2.col(ij_index) % g1.col(k_index)+
                                 f2.col(ik_index) % g1.col(j_index)+
                                 f1.col(i_index) % g2.col(jk_index)+
                                 f2.col(jk_index) % g1.col(i_index)+
                                 f1.col(j_index) % g2.col(ik_index)+
                                 f1.col(k_index) % g2.col(ij_index)+
                                 f0.col(i_index) % g3.col(ijk_index);
            
            for (int l = k; l < nparg; l++) {
              //Fourth derivative
              int l_index = trind(tri, {(unsigned int)l}) ;
              
              int il_index = trind(tri, {(unsigned int)i, (unsigned int)l}) ;
              int jl_index = trind(tri, {(unsigned int)j, (unsigned int)l}) ;
              int kl_index = trind(tri, {(unsigned int)k, (unsigned int)l}) ;
              
              int ijl_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)l}) ;
              int ikl_index = trind(tri, {(unsigned int)i, (unsigned int)k, (unsigned int)l}) ;
              int jkl_index = trind(tri, {(unsigned int)j, (unsigned int)k, (unsigned int)l}) ;
              
              int ijkl_index = trind(tri, {(unsigned int)i, (unsigned int)j, (unsigned int)k, (unsigned int)l}) ;
              
              h4.col(ijkl_index) += f4.col(ijkl_index)%g0.col(i_index)+
                                    f3.col(ijk_index)%g1.col(l_index)+
                                    f3.col(ijl_index)%g1.col(k_index)+
                                    f2.col(ij_index)%g2.col(kl_index)+
                                    //next line
                                    f3.col(ikl_index)%g1.col(j_index)+
                                    f2.col(ik_index)%g2.col(jl_index)+
                                    f2.col(il_index)%g1.col(jk_index)+
                                    f1.col(i)%g3.col(jkl_index)+
                                    //next line
                                    f3.col(jkl_index)%g1.col(i_index)+
                                    f2.col(jk_index)%g2.col(il_index)+
                                    f2.col(jl_index)%g2.col(ik_index)+
                                    f1.col(j_index)%g3.col(ikl_index)+
                                    //next line
                                    f2.col(kl_index)%g2.col(ij_index)+
                                    f1.col(k_index)%g3.col(ijl_index)+
                                    f1.col(l_index)%g3.col(ijk_index)+
                                    f0.col(i_index)%g4.col(ijkl_index);
            }  
          }
        }
      }
    }
  }
  
  NumericMatrix out = list2derivs(List::create(h0, h1, h2, h3, h4), deriv_order) ;
  
  return out ;
}

//' Productrule
//'
//' Productrule for derivs objects.
//' 
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = f_1(\cdot) \cdot f_2(\cdot) ... \cdot f_M(x_{n1},x_{n2},...,x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the productrule is applied.
//' For more details see [trind()] and [trind_generator()].
//' 
//' @return Returns an object of class \code{derivs} for the function \eqn{h(\cdot)}.
//' 
//' @inheritParams chainrule
//' 
//' @examples 
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=2)
//' B_derivs<-derivs_transform(A, type="inv", par=0,  trind_generator(1), deriv_order=2)
//' productrule (list(A_derivs, B_derivs), trind_generator(1), deriv_order=2) #identity
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix productrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[0]; //first element of list
  
  for (int i = 0; i < n-1; i++) {
    out = productrule_bi(f_list[i+1], out, tri, deriv_order) ;
  }
  
  return out ;
} 


// [[Rcpp::export]]
Rcpp::NumericMatrix quotientrule_bi (Rcpp::NumericMatrix f, Rcpp:: NumericMatrix g, Rcpp::List tri, int deriv_order){
  
  NumericMatrix out = derivs_transform(f = differencerule(List::create(derivs_transform(f , "log", {0}, tri, deriv_order),
                                                                derivs_transform(g , "log", {0}, tri, deriv_order)),
                                                                tri, deriv_order), "exp", {0}, tri, deriv_order) ;
  return out ;
}

//' Quotientrule
//'
//' Quotientrule for derivs objects.
//' 
//' @details Let \eqn{f_m} be a function defined in [trind()], where \eqn{m \in {1,...,M}}.
//' Define \eqn{h((x_{n1},x_{n2},...,x_{nK})) = f_1(\cdot) / f_2(\cdot) ... / f_M(x_{n1},x_{n2},...,x_{nK}))}.
//' In order to get the derivatives of \eqn{h(\cdot)} w.r.t all parameters \eqn{x_{nk}}, the quotientrule is applied.
//' For more details see [trind()] and [trind_generator()].The values of the \code{derivs} objects must be positive.
//' Numerically not precise, but included for  reasons of completeness.
//' 
//' @return Returns an object of class \code{derivs} for the function \eqn{h(\cdot)}.
//' 
//' @inheritParams chainrule
//' 
//' @examples 
//' A<-matrix(c(1:9)/10, ncol=1)
//' A_derivs<-list2derivs(list(A, A^0, A^2, A^3, A^4), deriv_order=2)
//' B_derivs<-derivs_transform(A, type="inv", par=0,  trind_generator(1), deriv_order=2)
//' quotientrule (list(A_derivs, B_derivs), trind_generator(1), deriv_order=2) #A/(1/A)=A^2
//' 
//' @family derivs
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix quotientrule (Rcpp::List f_list, Rcpp::List tri, int deriv_order){
  int n = f_list.length() ; //get length of list
  
  NumericMatrix out = f_list[0]; //first element of list
  
  for (int i = 0; i < n-1; i++) {
    out = quotientrule_bi(f_list[i+1], out, tri, deriv_order) ;
  }
  
  return out ;
} 
