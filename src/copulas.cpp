#include <RcppArmadillo.h>
#include "utility_cpp.h"
#include "derivs_rules.h"
#include "ind2joint.h"
#include "trind.h"
#include "transform.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix  dcop_unrot_cpp (arma::vec u, arma::vec v, arma::vec p, Rcpp::String distr_cop, int deriv_order, List tri, bool logp){
  
  // Initialize
  int n = u.size() ;
  
  arma::mat d0 (n, 1)  ;
  arma::mat d1 (n, 3)  ;
  arma::mat d2 (n, 6)  ;
  arma::mat d3 (n, 10)  ;
  arma::mat d4 (n, 15)  ;
  
  Rcpp::StringVector distributions_cop = {"independent","normal","clayton","gumbel","frank","joe","amh"} ;
  
  //Independent copula
  if(distr_cop==distributions_cop(0)){ 
    d0 += 0 ;
    if(deriv_order>0){
      d1 += 0 ;
      d2 += 0 ;
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Normal copula
  if(distr_cop==distributions_cop(1)){ 
    d0 += (p%(p%pow(erfcinv(2*u),2) - 2*erfcinv(2*u)%erfcinv(2*v) + p%pow(erfcinv(2*v),2)))/(-1 + pow(p,2)) - log(1 - pow(p,2))/2. ;
    if(deriv_order>0){
      //First derivatives of za wrt x1,v,p
      d1.col(0) += (-2*exp(pow(erfcinv(2*u),2))%p*sqrt(datum::pi)%(p%erfcinv(2*u) - erfcinv(2*v)))/(-1 + pow(p,2));
      d1.col(1) += (-2*exp(pow(erfcinv(2*v),2))%p*sqrt(datum::pi)%(-erfcinv(2*u) + p%erfcinv(2*v)))/(-1 + pow(p,2));
      d1.col(2) += (p - pow(p,3) + 2*(p%erfcinv(2*u) - erfcinv(2*v))%(-erfcinv(2*u) + p%erfcinv(2*v)))/pow(-1 + pow(p,2),2);
      
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += (2*exp(2*pow(erfcinv(2*u),2))%p*datum::pi%(p + 2*p%pow(erfcinv(2*u),2) - 2*erfcinv(2*u)%erfcinv(2*v)))/(-1 + pow(p,2));
      d2.col(1) += (-2*exp(pow(erfcinv(2*u),2) + pow(erfcinv(2*v),2))%p*datum::pi)/(-1 + pow(p,2));
      d2.col(2) += (-2*exp(pow(erfcinv(2*u),2))*sqrt(datum::pi)%(-2*p%erfcinv(2*u) + (1 + pow(p,2))%erfcinv(2*v)))/pow(-1 + pow(p,2),2);
      d2.col(3) += (2*exp(2*pow(erfcinv(2*v),2))%p*datum::pi%(p - 2*erfcinv(2*u)%erfcinv(2*v) + 2*p%pow(erfcinv(2*v),2)))/(-1 + pow(p,2)) ;
      d2.col(4) += (-2*exp(pow(erfcinv(2*v),2))*sqrt(datum::pi)%((1 + pow(p,2))%erfcinv(2*u) - 2*p%erfcinv(2*v)))/pow(-1 + pow(p,2),2);
      d2.col(5) += (-1 + pow(p,4) + (2 + 6*pow(p,2))%pow(erfcinv(2*u),2) - 4*p%(3 + pow(p,2))%erfcinv(2*u)%erfcinv(2*v) +(2 + 6*pow(p,2))%pow(erfcinv(2*v),2))/pow(-1 + pow(p,2),3);
        
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Clayton copula
  if(distr_cop==distributions_cop(2)){ 
    d0 += log(1 + p) - (1 + p)%log(u%v) + (-2 - 1/p)%log(-1 + pow(u,-p) + pow(v,-p)) ;
    if(deriv_order>0){
      //First derivatives of za wrt u,v,p
      d1.col(0) += -((p%pow(v,p) + (1 + p)%pow(u,p)%(-1 + pow(v,p)))/(u%(-pow(v,p) + pow(u,p)%(-1 + pow(v,p))))) ;
      d1.col(1) +=  (-1 - p + ((1 + 2*p)%pow(u,p))/(pow(v,p) - pow(u,p)%(-1 + pow(v,p))))/v;
      d1.col(2) +=  1/(1 + p) + ((1 + 2*p)%(pow(v,p)%log(u) + pow(u,p)%log(v)))/(p%(pow(v,p) - pow(u,p)%(-1 + pow(v,p)))) - log(u%v) + log(-1 + pow(u,-p) + pow(v,-p))/pow(p,2);
      
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += (-(p%pow(v,2*p)) + (-1 + p + 2*pow(p,2))%pow(u,p)%pow(v,p)%(-1 + pow(v,p)) + (1 + p)%pow(u,2*p)%pow(-1 + pow(v,p),2))/(pow(u,2)%pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2));
      d2.col(1) += (p%(1 + 2*p)%pow(u,-1 + p)%pow(v,-1 + p))/pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2);
      d2.col(2) += (pow(v,2*p) - pow(u,2*p)%pow(-1 + pow(v,p),2) + (1 + 2*p)%pow(u,p)%pow(v,p)%((-1 + pow(v,p))%log(u) + log(v)))/(u%pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2));
      d2.col(3) += (1 + p + (p%(1 + 2*p)%pow(u,2*p))/pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2) + ((1 + p)%(1 + 2*p)%pow(u,p))/(-pow(v,p) + pow(u,p)%(-1 + pow(v,p))))/pow(v,2);
      d2.col(4) += (pow(u,2*p) - pow(-1 + pow(u,p),2)%pow(v,2*p) + (1 + 2*p)%pow(u,p)%pow(v,p)%(log(u) + (-1 + pow(u,p))%log(v)))/(v%pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2));
      d2.col(5) += -pow(1 + p,-2) + (2*(pow(v,p)%log(u) + pow(u,p)%log(v)))/(pow(p,2)%(-pow(v,p) + pow(u,p)%(-1 + pow(v,p)))) + ((1 + 2*p)%pow(pow(v,p)%log(u) + pow(u,p)%log(v),2))/(p%pow(pow(v,p) - pow(u,p)%(-1 + pow(v,p)),2)) + ((1 + 2*p)%(pow(v,p)%pow(log(u),2) + pow(u,p)%pow(log(v),2)))/(-(p%pow(u,p)) + p%(-1 + pow(u,p))%pow(v,p)) - (2*log(-1 + pow(u,-p) + pow(v,-p)))/pow(p,3);
      
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Gumbel copula
  if(distr_cop==distributions_cop(3)){
    arma::vec lu = log(u) ;
    arma::vec lv = log(v) ;
    arma::vec lmlu = log(-lu) ;
    arma::vec lmlv = log(-lv) ;
    arma::vec mlu_p = pow(-lu,p) ;
    arma::vec mlv_p = pow(-lv,p) ;
    
    d0 +=  -pow(mlu_p + mlv_p,1/p) - log(u%v) + log(1 + (-1 + p)/pow(pow(-log(u),p) + mlv_p,1/p)) + ((-1 + p)%(-2*log(mlu_p + mlv_p) + p%log(lu%lv)))/p;
    if(deriv_order>0){
      //First derivatives of za wrt u,v,p
      d1.col(0) += (-1 + (-((mlu_p%((-1 + p)%p + 2*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p)))/(-1 + p + pow(mlu_p + mlv_p,1/p))) + (-1 + p)%mlv_p)/(lu%(mlu_p + mlv_p)))/u;
      d1.col(1) += (-1 + ((-1 + p)%mlu_p - (((-1 + p)%p + 2*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%mlv_p)/(-1 + p + pow(mlu_p + mlv_p,1/p)))/((mlu_p + mlv_p)%lv))/v;
      d1.col(2) += 1/(-1 + p + pow(mlu_p + mlv_p,1/p)) + (-(p%mlu_p%(1 - 3*p + 2*pow(p,2) + 3*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%lmlu) + (1 - p + (-3 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%(mlu_p + mlv_p)%log(mlu_p + mlv_p) - p%(1 - 3*p + 2*pow(p,2) + 3*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%mlv_p%lmlv)/(pow(p,2)%(-1 + p + pow(mlu_p + mlv_p,1/p))%(mlu_p + mlv_p)) + log(lu%lv);
      
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += (1 - p + lu - p%lu + pow(lu,2) + ((-1 + p)%pow(-lu,2*p)%((-1 + p)%p%(-1 + 2*p) + (2 + 5*(-1 + p)%p)%pow(mlu_p + mlv_p,1/p) + 2*(-1 + 2*p)%pow(mlu_p + mlv_p,2/p) + pow(mlu_p + mlv_p,3/p)))/(pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2)) + (mlu_p%(1 - p + lu)%(1 - 3*p + 2*pow(p,2) + 3*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p)))/((-1 + p + pow(mlu_p + mlv_p,1/p))%(mlu_p + mlv_p)))/(pow(u,2)%pow(lu,2));
      d2.col(1) += ((-1 + p)%pow(-lu,-1 + p)%((-1 + p)%p%(-1 + 2*p) + (2 + 5*(-1 + p)%p)%pow(mlu_p + mlv_p,1/p) + 2*(-1 + 2*p)%pow(mlu_p + mlv_p,2/p) + pow(mlu_p + mlv_p,3/p))%pow(-lv,-1 + p))/(u%v%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2));
      d2.col(2) += (pow(p,2)%(-(pow(-lu,2*p)%(pow(-1 + p,2) + (-1 + 2*p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))) - mlu_p%pow(mlu_p + mlv_p,1/p)%mlv_p + pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(-lv,2*p)) + mlu_p%(-(mlu_p%(2 - 3*p + pow(p,2) + 2*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%pow(mlu_p + mlv_p,1/p)%(p%lmlu - log(mlu_p + mlv_p))) + mlv_p%(-(pow(-1 + p,2)%pow(p,2)%(-1 + 2*p)%(lmlu - lmlv)) - 2*(-1 + p)%pow(mlu_p + mlv_p,2/p)%(2*pow(p,2)%lmlu - log(mlu_p + mlv_p) + (1 - 2*p)%p%lmlv) + pow(mlu_p + mlv_p,3/p)%(-(pow(p,2)%lmlu) + log(mlu_p + mlv_p) + (-1 + p)%p%lmlv) - (-1 + p)%pow(mlu_p + mlv_p,1/p)%(pow(p,2)%(-4 + 5*p)%lmlu - (-2 + p)%log(mlu_p + mlv_p) + p%(-2 - 5*(-1 + p)%p)%lmlv))))/(pow(p,2)%u%lu%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2));
      d2.col(3) += (1 + (-((-1 + p)%pow(-lu,2*p)%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%(1 + lv)) - mlu_p%(-1 + p + pow(mlu_p + mlv_p,1/p))%mlv_p%(1 - 3*pow(p,2) + 2*pow(p,3) + pow(mlu_p + mlv_p,2/p)%(-1 + p - lv) + (-1 + p)%pow(mlu_p + mlv_p,1/p)%(-1 + 3*p - lv) + lv - p%lv) +pow(-lv,2*p)%(pow(mlu_p + mlv_p,3/p)%lv + pow(-1 + p,2)%p%(1 + lv) + (-1 + p)%pow(mlu_p + mlv_p,2/p)%(1 + 3*lv) + (-1 + p)%pow(mlu_p + mlv_p,1/p)%(2*p + (-2 + 3*p)%lv)))/(pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2)%pow(lv,2)))/pow(v,2);
      d2.col(4) += (pow(p,2)%(pow(-lu,2*p)%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2) - mlu_p%pow(mlu_p + mlv_p,1/p)%mlv_p - (pow(-1 + p,2) + (-1 + 2*p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%pow(-lv,2*p)) + mlv_p%((2 - 3*p + pow(p,2) + 2*(-1 + p)%pow(mlu_p + mlv_p,1/p) + pow(mlu_p + mlv_p,2/p))%pow(mlu_p + mlv_p,1/p)%mlv_p%(log(mlu_p + mlv_p) - p%lmlv) + mlu_p%(pow(-1 + p,2)%pow(p,2)%(-1 + 2*p)%(lmlu - lmlv) + 2*(-1 + p)%pow(mlu_p + mlv_p,2/p)%(p%(-1 + 2*p)%lmlu + log(mlu_p + mlv_p) - 2*pow(p,2)%lmlv) + pow(mlu_p + mlv_p,3/p)%((-1 + p)%p%lmlu + log(mlu_p + mlv_p) - pow(p,2)%lmlv) + (-1 + p)%pow(mlu_p + mlv_p,1/p)%(p%(2 + 5*(-1 + p)%p)%lmlu + (-2 + p)%log(mlu_p + mlv_p) + (4 - 5*p)%pow(p,2)%lmlv))))/(pow(p,2)%v%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2)%lv );
      d2.col(5) += (4*log(mlu_p + mlv_p))/pow(p,3) - (4*(mlu_p%lmlu + mlv_p%lmlv))/(pow(p,2)%(mlu_p + mlv_p)) + 
        (2*(-1 + p)%pow(mlu_p%lmlu + mlv_p%lmlv,2))/(p%pow(mlu_p + mlv_p,2)) - (pow(mlu_p + mlv_p,-2 + 1/p)%pow(p%mlu_p%lmlu - (mlu_p + mlv_p)%log(mlu_p + mlv_p) + p%mlv_p%lmlv,2))/pow(p,4) + 
        ((-2 + 2/p)%(mlu_p%pow(lmlu,2) + mlv_p%pow(lmlv,2)))/(mlu_p + mlv_p) - pow(mlu_p%(pow(p,2) - (-1 + p)%p%lmlu + (-1 + p)%log(mlu_p + mlv_p)) + mlv_p%(pow(p,2) + (-1 + p)%(log(mlu_p + mlv_p) - p%lmlv)),2)/
          (pow(p,4)%pow(-1 + p + pow(mlu_p + mlv_p,1/p),2)%pow(mlu_p + mlv_p,2)) + (pow(mlu_p + mlv_p,-2 + 1/p)%(-(pow(p,2)%mlu_p%mlv_p%pow(lmlu,2)) - 2*pow(mlu_p + mlv_p,2)%log(mlu_p + mlv_p) + 
            p%mlv_p%lmlv%(2*mlv_p + mlu_p%(2 - p%lmlv)) + 2*p%mlu_p%lmlu%(mlu_p + mlv_p%(1 + p%lmlv))))/pow(p,3) + 
            (2*pow(p,2)%(log(mlu_p + mlv_p) - (p%(mlu_p%lmlu + mlv_p%lmlv))/(mlu_p + mlv_p)) + 
            (-1 + p)%pow(log(mlu_p + mlv_p) - (p%(mlu_p%lmlu + mlv_p%lmlv))/(mlu_p + mlv_p),2) + 
            (-1 + p)%p%(-2*log(mlu_p + mlv_p) + (p%(2*pow(-lu,2*p)%lmlu + 2*pow(-lv,2*p)%lmlv + 
            mlu_p%mlv_p%(-(p%pow(lmlu,2)) + lmlv%(2 - p%lmlv) + 2*lmlu%(1 + p%lmlv))))/pow(mlu_p + mlv_p,2)))/
              (pow(p,4)%(1 + (-1 + p)/pow(mlu_p + mlv_p,1/p))%pow(mlu_p + mlv_p,1/p)) ;
        
        
        
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Frank copula
  if(distr_cop==distributions_cop(4)){ 
    arma::vec monepu = -1 + u ;
    arma::vec monepv = -1 + v ;
    
    d0 += -(p%(u + v)) + log(1 - exp(-p)) - log(pow(exp(-p) - exp(-(p%u)) + exp(-(p%(u + v))),2)) +  log(p) ;//-(p%(u + v)) + log(1 - exp(-p)) - log(pow(exp(-p) - exp(-(p%u)) + exp(-(p%(u + v))),2)) + log(p) ;
    if(deriv_order>0){
      //First derivatives of za wrt u,v,p
      d1.col(0) += -(((-exp(p) + exp(p%(u + v)) + exp(p + p%v))%p)/(exp(p) + exp(p%(u + v)) - exp(p + p%v)));
      d1.col(1) += ((exp(p) - exp(p%(u + v)) + exp(p + p%v))%p)/(exp(p) + exp(p%(u + v)) - exp(p + p%v));
      d1.col(2) += 1/(-1 + exp(p)) + 1/p - u - v + (2*(exp(p%(u + v)) - exp(p + p%v)%u + exp(p)%(u + v)))/(exp(p) + exp(p%(u + v)) - exp(p + p%v));
      
        
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += (2*exp(p%(1 + u + v))%(-1 + exp(p%v))%pow(p,2))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2);
      d2.col(1) += (-2*exp(p%(1 + u + v))%pow(p,2))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2);
      d2.col(2) += (exp(2*p) + exp(2*p%(1 + v)) - 2*exp(p%(2 + v)) - exp(2*p%(u + v)) + 2*exp(p%(1 + u + 2*v))%p%(monepu) - 2*exp(p%(1 + u + v))%p%(monepu + v))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2);
      d2.col(3) += (2*exp(p + p%v)%(exp(p) - exp(p%u))%pow(p,2))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2);
      d2.col(4) += (exp(2*p) - exp(2*p%(1 + v)) - exp(2*p%(u + v)) + 2*exp(p%(1 + u + 2*v)) + 2*exp(p%(2 + v))%p%v - 2*exp(p%(1 + u + v))%p%(monepu + v))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2);
      d2.col(5) += 1/(1 - exp(p)) - pow(-1 + exp(p),-2) - pow(p,-2) + (2*pow(exp(p%(u + v)) - exp(p + p%v)%u + exp(p)%(u + v),2))/pow(exp(p) + exp(p%(u + v)) - exp(p + p%v),2) - (2*(exp(p%(u + v)) - exp(p + p%v)%pow(u,2) + exp(p)%pow(u + v,2)))/(exp(p) + exp(p%(u + v)) - exp(p + p%v));
      
     
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //Joe copula
  if(distr_cop==distributions_cop(5)){ 
    arma::vec pomup = pow(1 - u,p) ;
    arma::vec pomvp = pow(1 - v,p) ;
    arma::vec monepu = -1 + u ;
    arma::vec monepv = -1 + v ;
    
    d0 +=  log(pow(1 - u,-1 + p)%(p - (-1 + pomup)%(-1 + pomvp))%pow(pomup - (-1 + pomup)%pomvp,-2 + 1/p)%pow(1 - v,-1 + p));
    
    if(deriv_order>0){
      //First derivatives of za wrt u,v,p
      d1.col(0) += ((-1 + p)%(-(p%pomup) + (-1 + p + (1 + p)%pomup)%pomvp - (-1 + pomup)%pow(1 - v,2*p)))/((monepu)%(-p + (-1 + pomup)%(-1 + pomvp))%(-pomup + (-1 + pomup)%pomvp));
      d1.col(1) += ((p%(-1 + p + pomup))/(-p + (-1 + pomup)%(-1 + pomvp)) + ((1 - 2*p)%pomup)/(-pomup + (-1 + pomup)%pomvp))/(monepv);
      d1.col(2) += -((p%((-pow(-1 + p,2) + pomup)%pomup + ((-1 + p)%p + (2 + (-1 + p)%p)%pomup - 2*pow(1 - u,2*p))%pomvp + (-1 + pomup)%(-p + pomup)%pow(1 - v,2*p))%log(1 - u) - (-p + (-1 + pomup)%(-1 + pomvp))%(-pomup + (-1 + pomup)%pomvp)%log(pomup - (-1 + pomup)%pomvp) + p%(pow(-1 + pomup,2)%pow(1 - v,2*p)%log(1 - v) + p%pomup%(1 + (-1 + p + pomup)%log(1 - v)) - (-1 + pomup)%pomvp%(p + (-pow(-1 + p,2) + (1 + p)%pomup)%log(1 - v))))/(pow(p,2)%(-p + (-1 + pomup)%(-1 + pomvp))%(pomup - (-1 + pomup)%pomvp)));
      
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += ((-1 + p)%(p%(-1 + p + (1 + p)%pomup)%pow(1 - u,2*p) - (pow(-1 + p,2)%(1 + 2*p) + 2*(-1 + 2*(-1 + p)%p)%pomup + (1 + p)%(1 + 3*p)%pow(1 - u,2*p))%pomup%pomvp + (-pow(-1 + p,2) + (-1 + p)%(-5 + 2*(-2 + p)%p)%pomup + (-7 + p%(-7 + 5*p))%pow(1 - u,2*p) + 3*pow(1 + p,2)%pow(1 - u,3*p))%pow(1 - v,2*p) - (-1 + pomup)%(2 - 2*p + (1 + p)%(-5 + 3*p)%pomup + (1 + p)%(3 + p)%pow(1 - u,2*p))%pow(1 - v,3*p) + pow(-1 + pomup,2)%(-1 + (1 + p)%pomup)%pow(1 - v,4*p)))/(pow(monepu,2)%pow(p - (-1 + pomup)%(-1 + pomvp),2)%pow(pomup - (-1 + pomup)%pomvp,2));
      d2.col(1) += -(((-1 + p)%p%pow(1 - u,-1 + p)%(-2*pow(p,2) - pow(-1 + pomup,2)%pow(-1 + pomvp,2) + p%(-1 + pomup)%(-1 + pomvp)%(3 - pomup + (-1 + pomup)%pomvp))%pow(1 - v,-1 + p))/(pow(p - (-1 + pomup)%(-1 + pomvp),2)% pow(pomup - (-1 + pomup)%pomvp,2)));
      d2.col(2) += (-(((-1 + pomup)%pow(-1 + pomvp,2)%(pomup + (-1 + pomup)%pomvp) + pow(p,2)%(-pomup + (1 + pomup)%pomvp) - 2*p%(-1 + pomup)%(-1 + pomvp)%(-pomup + (1 + pomup)%pomvp))%(-pomup + (-1 + pomup)%pomvp)) + (-1 + p)%pomup%(-1 + pomvp)%(-(p%(-1 + pomvp)%(-pow(1 - u,2*p) + (-1 + pomup)%(3 + pomup)%pomvp)) + 2*pow(p,2)%pomvp + pow(-1 + pomup,2)%pow(-1 + pomvp,2)%pomvp)%log(1 - u) - (-1 + p)%pomup%(-2*pow(p,2) - pow(-1 + pomup,2)%pow(-1 + pomvp,2) + p%(-1 + pomup)%(-1 + pomvp)%(3 - pomup + (-1 + pomup)%pomvp))%pomvp%log(1 - v))/((monepu)%pow(p - (-1 + pomup)%(-1 + pomvp),2)% pow(pomup - (-1 + pomup)%pomvp,2));
      d2.col(3) += ((-1 + p)%(-(pow(-1 + p + pomup,2)%pow(1 - u,2*p)) + (-1 + pomup)%(-1 + p + pomup)%(-1 - p + 2*pow(p,2) + (3 + p)%pomup)%pomup%pomvp - pow(-1 + pomup,2)%(-((-1 + p)%p) + 2*(-1 + (-1 + p)%p)%pomup + (3 + 2*p)%pow(1 - u,2*p))%pow(1 - v,2*p) + (1 + p)%pow(-1 + pomup,3)%(-p + pomup)%pow(1 - v,3*p)))/(pow(p - (-1 + pomup)%(-1 + pomvp),2)%pow(pomup - (-1 + pomup)%pomvp,2)%pow(monepv,2));
      d2.col(4) += (pow(-1 + p + pomup,2)%pow(1 - u,2*p) + (-1 + p)%pomup%(2*pow(p,2) + pow(-1 + pomup,2)%pow(-1 + pomvp,2) - p%(-1 + pomup)%(-1 + pomvp)%(3 - pomup + (-1 + pomup)%pomvp))%pomvp%log(1 - u) + pow(-1 + pomup,3)%pow(1 - v,3*p)%(-1 + 2*p - pomup + (-1 + p)%(-p + pomup)%log(1 - v)) + (-1 + pomup)%pomup%pomvp%((3 - 2*p - 3*pomup)%pomup + (-1 + p)%(-1 + p + pomup)%(-1 + 2*p + pomup)%log(1 - v)) - pow(-1 + pomup,2)%pow(1 - v,2*p)%(pow(-1 + p,2) + pow(1 - u,2*p)%(-3 + 2*(-1 + p)%log(1 - v)) + 2*pomup%(p + pow(-1 + p,2)%log(1 - v))))/(pow(p - (-1 + pomup)%(-1 + pomvp),2)%pow(pomup - (-1 + pomup)%pomvp,2)%(monepv));
      d2.col(5) += (-(pow(p,3)%pow(pomup - (-1 + pomup)%pomvp,2)) + (-1 + p)%pow(p,2)%pomup%(-1 + pomvp)% (-(p%(-1 + pomvp)%(-pow(1 - u,2*p) + (-1 + pomup)%(3 + pomup)%pomvp)) + 2*pow(p,2)%pomvp + pow(-1 + pomup,2)%pow(-1 + pomvp,2)%pomvp)%pow(log(1 - u),2) + 2*pow(p - (-1 + pomup)%(-1 + pomvp),2)% pow(pomup - (-1 + pomup)%pomvp,2)%log(pomup - (-1 + pomup)%pomvp) - 2*p%pomup%log(1 - u)%((-2*p%(-1 + pomup)%(-1 + pomvp) + pow(-1 + pomup,2)%pow(-1 + pomvp,2) + pow(p,2)%(1 + pomup - (-1 + pomup)%pomvp))%(-1 + pomvp)%(-pomup + (-1 + pomup)%pomvp) + (-1 + p)%p%(-2*pow(p,2) - pow(-1 + pomup,2)%pow(-1 + pomvp,2) + p%(-1 + pomup)%(-1 + pomvp)%(3 - pomup + (-1 + pomup)%pomvp))%pomvp%log(1 - v)) + p%(-1 + pomup)%pomvp%log(1 - v)%(-2*pow(-1 + pomup,3)%pow(1 - v,3*p) + pow(-1 + pomup,2)%pow(1 - v,2*p)%(-4 + 2*p%(2 + p) + 6*pomup + (-1 + p)%p%(-p + pomup)%log(1 - v)) + pomup%(2*(pow(-1 + p,2) + (-2 + p%(2 + p))%pomup + pow(1 - u,2*p)) + (-1 + p)%p%(-1 + p + pomup)%(-1 + 2*p + pomup)%log(1 - v)) - 2*(-1 + pomup)%pomvp%(pow(-1 + p,2) + pow(1 - u,2*p)%(3 + (-1 + p)%p%log(1 - v)) + pomup%(-4 + 2*p%(2 + p) + pow(-1 + p,2)%p%log(1 - v)))))/(pow(p,3)%pow(p - (-1 + pomup)%(-1 + pomvp),2)% pow(pomup - (-1 + pomup)%pomvp,2));
      
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  //AMH copula
  if(distr_cop==distributions_cop(6)){ 
    arma::vec monepu = -1 + u ;
    arma::vec monepv = -1 + v ;
    arma::vec p_2 = pow(p,2) ;
    arma::vec uv = u%v ;
    arma::vec monepumonepv = monepu%monepv ;
    
    d0 += log(-((1 + p%(-2 + u + p%monepumonepv + v + uv))/pow(-1 + p%monepumonepv,3))) ;
    if(deriv_order>0){
      //First derivatives of za wrt u,v,p
      d1.col(0) += (-2*p%(-1 + p_2%(monepu)%pow(monepv,2) + 2*v + p%(monepv)*(-2 + u + 2*v + uv)))/
        ((-1 + p%monepumonepv)*(1 + p_2%monepumonepv + p%(-2 + u + v + uv))) ;
      d1.col(1) += (-2*p%(-1 + 2*u + p_2%pow(monepu,2)%(monepv) + p%(monepu)%(-2 + v + u%(2 + v))))/
        ((-1 + p%monepumonepv)%(1 + p_2%monepumonepv + p%(-2 + u + v + uv))) ;
      d1.col(2) += (monepu%(2 - 4*v) - p_2%pow(monepu,2)%pow(monepv,2) + 2*v - 2*p%monepumonepv%(monepu + v + uv))/
        ((-1 + p%monepumonepv)%(1 + p_2%monepumonepv + p%(-2 + u + v + uv))) ;
      
      //Second derivatives of za wrt uu, uv, up, vv, vp, pp
      d2.col(0) += (2*p_2%(1 + pow(p,4)%pow(monepu,2)%pow(monepv,4) - 4*v + pow(v,2) + 
        2*pow(p,3)%(monepu)%pow(monepv,3)*(-2 + u + 2*v + uv) + 
        2*p%(monepv)%(2 - 6*v + pow(v,2) + u%(monepv + 2*pow(v,2))) + 
        p_2%pow(monepv,2)%(6 - 12*v + pow(v,2) + pow(u,2)%pow(1 + v,2) + u%(-6 + 4*v + 4*pow(v,2))))
      )/(pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      d2.col(1) += (4*p%(1 + pow(p,4)%pow(monepu,2)%pow(monepv,2) + p%(2*(-2 + v) + u%(2 + v)) + 
        pow(p,3)%monepumonepv%(2*(-2 + v) + u%(2 + v)) + 
        p_2%(6 - 6*v + pow(v,2) + pow(u,2)%(1 + v + pow(v,2)) + u%(-6 + 2*v + pow(v,2)))))/
          (pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      d2.col(2) += (2*(-1 + 2*v + pow(p,4)%pow(monepu,2)%pow(monepv,3)%(1 + v) + 2*p%(monepv)%(-2 + u + 2*v + uv) + 
        2*pow(p,3)%(monepu)%pow(monepv,2)%(-2 + u + 2*uv) + 
        p_2%(monepv)%(2*(3 - 3*v + pow(v,2)) + 2*u%(-3 - v + pow(v,2)) + 
        pow(u,2)%(1 + 3*v + 2*pow(v,2)))))/
          (pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      d2.col(3) += (2*p_2%(1 - 4*u + pow(u,2) + pow(p,4)%pow(monepu,4)%pow(monepv,2) + 
        2*pow(p,3)%pow(monepu,3)%(monepv)%(-2 + v + u%(2 + v)) + 
        2*p%(monepu)%(2 + u*(-6 + v) - v + pow(u,2)%(1 + 2*v)) + 
        p_2%pow(monepu,2)%(6 - 6*v + pow(v,2) + 2*u%(-6 + 2*v + pow(v,2)) + 
        pow(u,2)%(1 + 4*v + pow(v,2)))))/
          (pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      d2.col(4) += (2*(-1 + 2*u + pow(p,4)%pow(monepu,3)%(1 + u)%pow(monepv,2) + 
        2*pow(p,3)%pow(monepu,2)%(monepv)%(-2 + v + 2*uv) + 2*p%(monepu)%(-2 + v + u%(2 + v)) + 
        p_2%(monepu)%(6 - 6*v + pow(v,2) + 2*pow(u,2)%(1 + v + pow(v,2)) + u%(-6 - 2*v + 3*pow(v,2))))
      )/(pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      d2.col(5) += (1 + pow(p,4)%pow(monepu,4)%pow(monepv,4) - 4*v + 2*pow(v,2) + 
        4*pow(p,3)%pow(monepu,3)%pow(monepv,3)%(monepu + v + uv) + 2*pow(u,2)%(1 - 4*v + pow(v,2)) - 
        4*u%(1 - 4*v + 2*pow(v,2)) + 2*p_2%pow(monepu,2)%pow(monepv,2)%
        (3 - 6*v + pow(v,2) + pow(u,2)%pow(1 + v,2) + 2*u%(-3 + 2*v + pow(v,2))) + 
        4*p%monepumonepv%(-1 + 3*v - pow(v,2) + u%(3 - 7*v + pow(v,2)) + pow(u,2)%(monepv + 2*pow(v,2))))/
          (pow(-1 + p%monepumonepv,2)%pow(1 + p_2%monepumonepv + p%(-2 + u + v + uv),2)) ;
      
      if(deriv_order>2){
        d3 += 0 ;
        d4 += 0 ;
      }
    }
  }
  
  Rcpp::NumericMatrix out = list2derivs(List::create(d0, d1, d2, d3, d4), deriv_order) ;
  
  if(!logp){
    out = derivs_transform(out, "exp", {0}, tri, deriv_order) ;
  }
  
  return out ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix dcop_cpp (arma::vec u, arma::vec v, arma::vec p, Rcpp::String distr_cop, int rot, int deriv_order, List tri, bool logp){
  // Initialize
  int n = u.size() ;
  
  Rcpp::StringVector rotation_vector ; //= {"identity", "identity", "identity"} ;
  
  if(rot!=0){
    if(rot==90){
      rotation_vector = {"onemx", "identity", "identity"} ;
    }
    
    if(rot==180){
      rotation_vector = {"onemx", "onemx", "identity"} ;
    }
    
    if(rot==270){
      rotation_vector = {"identity", "onemx", "identity"} ;
    }
  }
  
  //Initialize output
  NumericMatrix out(n, 1) ;
  Rcpp::List tri_1 = trind_generator(1) ;
  NumericMatrix cop_unrot = dcop_unrot_cpp (u, v, p, distr_cop, deriv_order, tri, logp) ;
  if(rot==0){
    out = cop_unrot ;
  } else {
    out = chainrule(List::create(cop_unrot,ind2joint(List::create(transform(u, rotation_vector(0), {0}, deriv_order),
                                                                  transform(v, rotation_vector(1), {0}, deriv_order),
                                                                  transform(p, rotation_vector(2), {0}, deriv_order)),
                                                     List::create(tri_1,tri_1,tri_1),
                                                     List::create(tri_1,trind_generator(2),tri),
                                                     deriv_order)),
                                                     tri, deriv_order) ;
  }
  
  return out ;
}
