#include <RcppArmadillo.h>
#include "utility_cpp.h"
#include "derivs_rules.h"
#include "trind.h"
#include "transform.h"
#include "copulas.h"
#include "ind2joint.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix dnormhnorm_cpp (arma::vec x, arma::vec m, arma::vec v, arma::vec u, int s, int deriv_order, List tri, bool logp){
  
  // Initialize
  int n = x.size() ;
  
  arma::vec mx = m-x ;
  arma::vec mx2 = arma::pow(mx,2) ;
  
  arma::vec v2 = arma::pow(v,2) ;
  arma::vec v3 = arma::pow(v,3) ;
  arma::vec v4 = arma::pow(v,4) ;
  
  arma::vec u2 = arma::pow(u,2) ;
  arma::vec u3 = arma::pow(u,3) ;
  arma::vec u4 = arma::pow(u,4) ;
  
  arma::vec v2u2 = u2 + v2 ;
  arma::vec v2u2_15 = arma::pow(v2u2,1.5) ;
  arma::vec v2u2_2 = arma::pow(v2u2,2) ;
  arma::vec v2u2_25 = arma::pow(v2u2,2.5) ;
  arma::vec v2u2_3 = arma::pow(v2u2,3) ;
  arma::vec v2u2_35 = arma::pow(v2u2,3.5) ;
  arma::vec v2u2_4 = arma::pow(v2u2,4) ;
  arma::vec v2u2_45 = arma::pow(v2u2,4.5) ;
  arma::vec v2u2_5 = arma::pow(v2u2,5) ;
  
  arma::mat za = -0.9189385332046727 - (0.5*mx2)/(v2u2) - 0.5*log(v2u2) ;
  za.reshape(n,1) ;
  arma::mat za_d1 (n, 3)  ;
  arma::mat za_d2 (n, 6)  ;
  arma::mat za_d3 (n, 10)  ;
  arma::mat za_d4 (n, 15)  ;
  
  arma::mat zb = s*(u%(-m + x))/(v%sqrt(v2u2)) ;
  zb.reshape(n,1) ;
  arma::mat zb_d1 (n, 3)  ;
  arma::mat zb_d2 (n, 6)  ;
  arma::mat zb_d3 (n, 10)  ;
  arma::mat zb_d4 (n, 15)  ;
  
  if(deriv_order>0){
    //First derivatives of za wrt m,v,u
    za_d1.col(0) += ((-m + x))/(v2u2) ;
    za_d1.col(1) += (v%(-(v2u2) + mx2))/v2u2_2 ;
    za_d1.col(2) += (u%(-(v2u2) + mx2))/v2u2_2 ;
    
    //Second derivatives of za wrt mm,mv,mu,vv,vu,uu
    za_d2.col(0) += -1/(v2u2) ;
    za_d2.col(1) += (2*v%(mx))/v2u2_2 ;
    za_d2.col(2) += (2*u%(mx))/v2u2_2 ;
    za_d2.col(3) += (2*v2%(v2u2) - v2u2_2 - 4*v2%mx2 + (v2u2)%mx2)/v2u2_3;
    za_d2.col(4) += (u%v%(2*(v2u2) - 4*mx2))/v2u2_3;
    za_d2.col(5) += (2*u2%(v2u2) - v2u2_2 - 4*u2%mx2 + (v2u2)%mx2)/v2u2_3;
    
    //First derivatives of zb wrt m,v,u
    zb_d1.col(0) += -(u/(v%sqrt(v2u2))) ;
    zb_d1.col(1) += (u%(u2 + 2*v2)%(mx))/(v2%v2u2_15) ;
    zb_d1.col(2) += (v%(-m + x))/v2u2_15 ;
    
    //Second derivatives of zb wrt mm,mv,mu,vv,vu,uu
    zb_d2.col(0) += 0 ;
    zb_d2.col(1) += (u%(u2 + 2*v2))/(v2%v2u2_15) ;
    zb_d2.col(2) += -(v/v2u2_15) ;
    zb_d2.col(3) += -((u%(2*u4 + 5*u2%v2 + 6*v4)%(mx))/(v3%v2u2_25)) ;
    zb_d2.col(4) += -(((u2 - 2*v2)%(mx))/v2u2_25)  ;
    zb_d2.col(5) += (3*u%v%(mx))/v2u2_25  ;
    
    //Multiply with s
    zb_d1 *= s ;
    zb_d2 *=s ;
    
    if(deriv_order>2){
      //Third derivatives of za wrt mmm,mmv,mmu,mvv,mvu,muu,vvv,vvu,vuu,uuu
      za_d3.col(0) += 0 ;
      za_d3.col(1) += (2*v)/v2u2_2 ;
      za_d3.col(2) += (2*u)/v2u2_2 ;
      za_d3.col(3) += (2*(u2 - 3*v2)%(mx))/v2u2_3 ;
      za_d3.col(4) += (-8*u%v%(mx))/v2u2_3 ;
      za_d3.col(5) += (-6*(u2 - 0.3333333333333333*v2)%(mx))/v2u2_3 ;
      za_d3.col(6) += (v%(-8*v2%(v2u2) + 6*v2u2_2 + 24*v2%mx2 - 12*(v2u2)%mx2))/v2u2_4 ;
      za_d3.col(7) += (u%(-8*v2%(v2u2) + 2*v2u2_2 + 24*v2%mx2 - 4*(v2u2)%mx2))/v2u2_4 ;
      za_d3.col(8) += (v%(-8*u2%(v2u2) + 2*v2u2_2 + 24*u2%mx2 - 4*(v2u2)%mx2))/v2u2_4 ;
      za_d3.col(9) += (u%(-8*u2%(v2u2) + 6*v2u2_2 + 24*u2%mx2 - 12*(v2u2)%mx2))/v2u2_4 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,vvuu,vuuu,uuuu
      za_d4.col(0) += 0 ;
      za_d4.col(1) += 0 ;
      za_d4.col(2) += 0 ;
      za_d4.col(3) += (2*u2 - 6*v2)/v2u2_3 ;
      za_d4.col(4) += (-8*u%v)/v2u2_3 ;
      za_d4.col(5) += (-6*u2 + 2*v2)/v2u2_3 ;
      za_d4.col(6) += (-24*(u2%v - v3)%(mx))/v2u2_4 ;
      za_d4.col(7) += (-8*u%(u2 - 5*v2)%(mx))/v2u2_4 ;
      za_d4.col(8) += (-8*v%(-5*v2u2)%(mx))/v2u2_4 ;
      za_d4.col(9) += (24*(u3 - u%v2)%(mx))/v2u2_4 ;
      za_d4.col(10) += (48*v4%(v2u2) - 48*v2%v2u2_2 + 6*v2u2_3 - 192*v4%mx2 + 144*v2%(v2u2)%mx2 - 12*v2u2_2%mx2)/v2u2_5 ;
      za_d4.col(11) += (u%v%(-24*u4 + 24*v4 + pow(m,2)%(72*u2 - 120*v2) + m%(-144*u2 + 240*v2)%x + 72*u2%pow(x,2) - 120*v2%pow(x,2)))/v2u2_5 ;
      za_d4.col(12) += (-6*pow(u,6) - 6*pow(v,6) + pow(m,2)%(20*u4 - 152*u2%v2 + 20*v4) + m%(-40*u4 + 304*u2%v2 - 40*v4)%x + 20*v4%pow(x,2) + u4%(30*v2 + 20*pow(x,2)) + u2%(30*v4 - 152*v2%pow(x,2)))/v2u2_5 ;//(u*v*(24*u4 - 24*v4 + pow(m,2)*(-120*u2 + 72*v2) + m*(240*u2 - 144*v2)*x - 120*u2*pow(x,2) + 72*v2*pow(x,2)))/v2u2_5 ;
      za_d4.col(13) += (u%v%(24*u4 - 24*v4 + pow(m,2)%(-120*u2 + 72*v2) + m%(240*u2 - 144*v2)%x - 120*u2%pow(x,2) + 72*v2%pow(x,2)))/v2u2_5 ;
      za_d4.col(14) += (48*u4%(v2u2) - 48*u2%v2u2_2 + 6*v2u2_3 - 192*u4%mx2 + 144*u2%(v2u2)%mx2 - 12*v2u2_2%mx2)/v2u2_5 ;
      
      //Third derivatives of zb wrt mmm,mmv,mmu,mvv,mvu, (5) muu,vvv,vvu,vuu,uuu
      zb_d3.col(0) += 0 ;
      zb_d3.col(1) += 0 ;
      zb_d3.col(2) += 0 ;
      zb_d3.col(3) += -((u%(2*u4 + 5*u2%v2 + 6*v4))/(v3%v2u2_25)) ;
      zb_d3.col(4) += (-u2 + 2*v2)/v2u2_25 ;
      zb_d3.col(5) += (3*u%v)/v2u2_25 ;
      zb_d3.col(6) += (3*u%(2*pow(u,6) + 7*u4%v2 + 8*u2%v4 + 8*pow(v,6))%(mx))/(v4%v2u2_35) ;
      zb_d3.col(7) += (3*v%(3*u2 - 2*v2)%(mx))/v2u2_35 ;
      zb_d3.col(8) += (3*(u3 - 4*u%v2)%(mx))/v2u2_35 ;
      zb_d3.col(9) += (3*v%(-4*v2u2)%(mx))/v2u2_35 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,(12) vvuu,vuuu,uuuu
      zb_d4.col(0) += 0 ;
      zb_d4.col(1) += 0 ;
      zb_d4.col(2) += 0 ;
      zb_d4.col(3) += 0 ;
      zb_d4.col(4) += 0 ;
      zb_d4.col(5) += 0 ;
      zb_d4.col(6) += (3*u%(2*pow(u,6) + 7*u4%v2 + 8*u2%v4 + 8*pow(v,6)))/(v4%v2u2_35) ;
      zb_d4.col(7) +=  (9*u2%v - 6*v3)/v2u2_35 ;
      zb_d4.col(8) +=  (3*(u3 - 4*u%v2))/v2u2_35 ;
      zb_d4.col(9) +=  (3*v%(-4*v2u2))/v2u2_35 ;
      zb_d4.col(10) +=  (-3*u%(8*pow(u,8) + 36*pow(u,6)%v2 + 63*u4%v4 + 40*u2%pow(v,6) + 40*pow(v,8))%(mx))/(pow(v,5)%v2u2_45) ;
      zb_d4.col(11) +=  (3*(3*u4 - 24*u2%v2 + 8*v4)%(mx))/v2u2_45 ;
      zb_d4.col(12) +=  (15*v%(-3*u3 + 4*u%v2)%(mx))/v2u2_45;//(-3*(4*u4 - 27*u2*v2 + 4*v4)*(mx))/v2u2_45 ;
      zb_d4.col(13) +=  (-3*(4*u4 - 27*u2%v2 + 4*v4)%(mx))/v2u2_45 ;
      zb_d4.col(14) +=  (15*v%(4*u3 - 3*u%v2)%(mx))/v2u2_45 ;
      
      //Multiply with s
      zb_d3 *= s ;
      zb_d4 *= s ;
    }
  }
  
  Rcpp::NumericMatrix out = sumrule(List::create(list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), deriv_order),
                                                 derivs_transform(list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), deriv_order), "zeta", {0}, tri, deriv_order)),
                                           tri, deriv_order) ;      
  
  if(!logp){
    out = derivs_transform(out, "exp", {0}, tri, deriv_order) ;
  }
  
  //List out = List::create(list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), deriv_order),
  //                        derivs_transform(list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), deriv_order), "zeta", {0}, tri, deriv_order)) ;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pnormhnorm_cpp (arma::vec q, arma::vec m, arma::vec v, arma::vec u, int s, int deriv_order, List tri, bool logp){
  
  // Initialize
  int n = q.size() ;
  
  arma::vec mq = m - q ;
  arma::vec mq2 = arma::pow(mq,2) ;
 
  arma::vec v2 = arma::pow(v,2) ;
  arma::vec v3 = arma::pow(v,3) ;
  arma::vec v4 = arma::pow(v,4) ;
  
  arma::vec u2 = arma::pow(u,2) ;
  arma::vec u3 = arma::pow(u,3) ;
  arma::vec u4 = arma::pow(u,4) ;
  
  arma::vec v2u2 = u2 + v2 ;
  arma::vec v2u2_15 = pow(v2u2,1.5) ;
  arma::vec v2u2_2 = pow(v2u2,2) ;
  arma::vec v2u2_25 = pow(v2u2,2.5) ;
  arma::vec v2u2_3 = pow(v2u2,3) ;
  arma::vec v2u2_35 = pow(v2u2,3.5) ;
  arma::vec v2u2_4 = pow(v2u2,4) ;
  arma::vec v2u2_45 = pow(v2u2,4.5) ;
  arma::vec v2u2_5 = pow(v2u2,5) ;
  
  double twopi_05 = sqrt(2*datum::pi) ;
    
  arma::vec za =  (-m + q)/sqrt(v2u2);
  za.reshape(n,1) ;
  arma::mat za_d1 (n, 3)  ;
  arma::mat za_d2 (n, 6)  ;
  arma::mat za_d3 (n, 10)  ;
  arma::mat za_d4 (n, 15)  ;
  
  arma::vec zb = -2*OwenT((-m + q)/sqrt(v2u2), (s*u)/v);
  zb.reshape(n,1) ;
  arma::mat zb_d1 (n, 3)  ;
  arma::mat zb_d2 (n, 6)  ;
  arma::mat zb_d3 (n, 10)  ;
  arma::mat zb_d4 (n, 15)  ;
  
  if(deriv_order>0){
    //First derivatives of za wrt m,v,u
    za_d1.col(0) += -(1/sqrt(v2u2)) ;
    za_d1.col(1) += ((mq)%v)/v2u2_15 ;
    za_d1.col(2) += ((mq)%u)/v2u2_15 ;
    
    //Second derivatives of za wrt mm,mv,mu,vv,vu,uu
    za_d2.col(0) += 0 ;
    za_d2.col(1) += v/v2u2_15 ;
    za_d2.col(2) += u/v2u2_15 ;
    za_d2.col(3) += ((mq)%(u2 - 2*v2))/v2u2_25 ;
    za_d2.col(4) += (3*(-m + q)%u%v)/v2u2_25 ;
    za_d2.col(5) += -(((mq)%(2*u2 - v2))/v2u2_25) ;
    
    
    //First derivatives of zb wrt m,v,u
    zb_d1.col(0) += -(erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2)))/(exp(mq2/(2.*(v2u2)))*twopi_05%sqrt(v2u2)));
    zb_d1.col(1) += (s*u)/(exp((mq2%(1 + (u2)/v2))/(2.*(v2u2)))*datum::pi%(v2u2)) + ((mq)%v%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(exp(mq2/(2.*(v2u2)))*twopi_05%v2u2_15);
    zb_d1.col(2) += -((s*v)/(exp((mq2%(1 + (u2)/v2))/(2.*(v2u2)))*datum::pi%(v2u2))) + ((mq)%u%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(exp(mq2/(2.*(v2u2)))*twopi_05%v2u2_15);
    
    //Second derivatives of zb wrt mm,mv,mu,vv,vu,uu
    zb_d2.col(0) += ((2*s*u%sqrt(v2u2))/exp((mq2*1%u2)/(2.*v2%(v2u2))) + twopi_05*(mq)%v%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(2.*exp(mq2/(2.*(v2u2)))*datum::pi%v%v2u2_15) ;
    zb_d2.col(1) += (-(twopi_05*(mq)*s%u%sqrt(v2u2)%(u2 + 2*v2)) + exp((mq2*1%u2)/(2.*v2%(v2u2)))*datum::pi%v3%(-mq2 + v2u2)%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(sqrt(2)*exp((mq2%(v2u2))/(2.*v2%(v2u2)))*pow(datum::pi,1.5)%v2%v2u2_25) ;
    zb_d2.col(2) += (twopi_05*(mq)*s%v%sqrt(v2u2) + exp((mq2%u2)/(2.*v2%(v2u2)))*datum::pi%u%(-mq2 + v2u2)%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(sqrt(2)*exp((mq2%(v2u2))/(2.*v2%(v2u2)))*pow(datum::pi,1.5)%v2u2_25) ;
    zb_d2.col(3) += (s*u%(mq2*1%pow(u,8) + mq2*(1 + 4*1)%pow(u,6)%v2 - 2*u4%(-(mq2*(3 + 2*1)) + u2)%v4 + u2%(mq2*(2 + 7*1) - 6*u2)%pow(v,6) + 3*(mq2 - 2*u2)%pow(v,8) - 2*pow(v,10)))/(exp((mq2%(1 + (u2)/v2))/(2.*(v2u2)))*datum::pi%v3%v2u2_3%v2u2_2) + ((mq)%(u4 - (mq + u)%(-m + q + u)%v2 - 2*v4)%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(exp(mq2/(2.*(v2u2)))*twopi_05%v2u2_35) ;
    zb_d2.col(4) += (s*(pow(u,6)%(-(mq2*(-1 + 1)) - u2) + u4%(-(mq2*(-1 + 2*1)) + (1 - 3*1)*u2)%v2 + u2%(-(mq2*(-1 + 3*1)) - 3*(-1 + 1)*u2)%v4 - (mq2 + (-3 + 1)*u2)%pow(v,6) + pow(v,8)))/(exp((mq2%(v2u2))/(2.*v2%(v2u2)))*datum::pi%v2u2_3%v2u2_2) + ((mq)%u%v%(mq2 - 3*(v2u2))%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(exp(mq2/(2.*(v2u2)))*twopi_05%v2u2_35) ;
    zb_d2.col(5) +=  (s*u%v%(u4%(-mq2 + 2*u2) + u2%(mq2*(-1 - 2*1 + 1) + 6*u2)%v2 + (mq2*(-2 + 1) + 6*u2)%v4 + 2*pow(v,6)))/(exp((mq2%(v2u2))/(2.*v2%(v2u2)))*datum::pi%v2u2_3%v2u2_2) + ((mq)%(mq2%u2 - 2*u4 - u2%v2 + v4)%erf(((-m + q)*s%u)/(sqrt(2)*v%sqrt(v2u2))))/(exp(mq2/(2.*(v2u2)))*twopi_05%v2u2_35) ;
    
    if(deriv_order>2){
      //Third derivatives of za wrt mmm,mmv,mmu,mvv,mvu,muu,vvv,vvu,vuu,uuu
      za_d3.col(0) += 0 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,vvuu,vuuu,uuuu
      za_d4.col(0) += 0 ;
      
      //Third derivatives of zb wrt mmm,mmv,mmu,mvv,mvu, (5) muu,vvv,vvu,vuu,uuu
      zb_d3.col(0) += 0 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,(12) vvuu,vuuu,uuuu
      zb_d4.col(0) += 0 ;
    }
  }
  
  NumericMatrix out = sumrule(List::create(derivs_transform(list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), deriv_order), "pnorm", {0}, tri, deriv_order),
                                           list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), deriv_order)),
                                           tri, deriv_order) ;      
  
  if(logp){
    out = derivs_transform(out, "log", {0}, tri, deriv_order) ;//chainrule(List::create(transform(as<arma::mat>(out), "exp", {0}, deriv_order), out), tri, deriv_order) ;
  }
  
  return out ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix dnormexp_cpp (arma::vec x, arma::vec m, arma::vec v, arma::vec u, int s, int deriv_order, List tri, bool logp){
  
  // Initialize
  int n = x.size() ;
  //arma::mat out (n, 1) ;
  arma::vec mx = m-x ;
  arma::vec mx2 = pow(mx,2) ;
  
  arma::vec v2 = arma::pow(v,2) ;
  arma::vec v3 = arma::pow(v,3) ;
  arma::vec v4 = arma::pow(v,4) ;
  
  arma::vec u2 = arma::pow(u,2) ;
  arma::vec u3 = arma::pow(u,3) ;
  arma::vec u4 = arma::pow(u,4) ;
    
  arma::vec za = u%(m*s + 0.5*u%v2 - s*x) + log(u/2) ;
  za.reshape(n,1) ;
  arma::mat za_d1 (n, 3)  ;
  arma::mat za_d2 (n, 6)  ;
  arma::mat za_d3 (n, 10)  ;
  arma::mat za_d4 (n, 15)  ;
  
  arma::vec zb = -((m*s + u%v2 - s*x)/v) ;
  arma::mat zb_d1 (n, 3)  ;
  arma::mat zb_d2 (n, 6)  ;
  arma::mat zb_d3 (n, 10)  ;
  arma::mat zb_d4 (n, 15)  ;
  
  //arma::mat zeta_zc_arma = as<arma::mat>(zeta_zc) ;
  
  if(deriv_order>0){
    //First derivatives of za wrt m,v,u
    za_d1.col(0) += u*s ;
    za_d1.col(1) += u2%v;
    za_d1.col(2) += 1/u + m*s + u%v2 - s*x;
    
    //Second derivatives of za wrt mm,mv,mu,vv,vu,uu
    za_d2.col(0) += 0 ;
    za_d2.col(1) += 0 ;
    za_d2.col(2) += s ;
    za_d2.col(3) += u2 ;
    za_d2.col(4) += 2*u%v ;
    za_d2.col(5) += -pow(u,-2) + v2 ;
    
    //First derivatives of zb wrt m,v,u
    zb_d1.col(0) += -(s/v) ;
    zb_d1.col(1) += (m*s - u%v2 - s*x)/v2 ;
    zb_d1.col(2) += -v ;
    
    //Second derivatives of zb wrt mm,mv,mu,vv,vu,uu
    zb_d2.col(0) += 0 ;
    zb_d2.col(1) += s/v2 ;
    zb_d2.col(2) += 0 ;
    zb_d2.col(3) += (2*s*(-m + x))/v3 ;
    zb_d2.col(4) += -1 ;
    zb_d2.col(5) +=  0 ;
    
    if(deriv_order>2){
      //Third derivatives of za wrt mmm,mmv,mmu,mvv,mvu,muu,vvv,vvu,(8)vuu,uuu
      za_d3.col(0) += 0 ;
      za_d3.col(1) += 0 ;
      za_d3.col(2) += 0 ;
      za_d3.col(3) += 0 ;
      za_d3.col(4) += 0 ;
      za_d3.col(5) += 0 ;
      za_d3.col(6) += 0 ;
      za_d3.col(7) += 2.*u ;
      za_d3.col(8) += 2.*v ;
      za_d3.col(9) += 2/u3 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,vvuu,vuuu,uuuu
      za_d4.col(0) += 0 ;
      za_d4.col(1) += 0 ;
      za_d4.col(2) += 0 ;
      za_d4.col(3) += 0 ;
      za_d4.col(4) += 0 ;
      za_d4.col(5) += 0 ;
      za_d4.col(6) += 0 ;
      za_d4.col(7) += 0 ;
      za_d4.col(8) += 0 ;
      za_d4.col(9) += 0 ;
      za_d4.col(10) += 0 ;
      za_d4.col(11) += 0 ;
      za_d4.col(12) += 2 ;
      za_d4.col(13) += 0 ;
      za_d4.col(14) += -6/u4 ;
      
      //Third derivatives of zb wrt mmm,mmv,mmu,mvv,mvu, (5) muu,vvv,vvu,vuu,uuu
      zb_d3.col(0) += 0 ;
      zb_d3.col(1) += 0 ;
      zb_d3.col(2) += 0 ;
      zb_d3.col(3) += (-2*s)/v3 ;
      zb_d3.col(4) += 0 ;
      zb_d3.col(5) += 0 ;
      zb_d3.col(6) += (6*s*(mx))/v4 ;
      zb_d3.col(7) += 0 ;
      zb_d3.col(8) += 0 ;
      zb_d3.col(9) += 0 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,(12) vvuu,vuuu,uuuu
      zb_d4.col(0) += 0 ;
      zb_d4.col(1) += 0 ;
      zb_d4.col(2) += 0 ;
      zb_d4.col(3) += 0 ;
      zb_d4.col(4) += 0 ;
      zb_d4.col(5) += 0 ;
      zb_d4.col(6) += (6*s)/v4 ;
      zb_d4.col(7) += 0 ;
      zb_d4.col(8) += 0 ;
      zb_d4.col(9) += 0 ;
      zb_d4.col(10) += (24*s*(-m + x))/pow(v,5) ;
      zb_d4.col(11) += 0 ;
      zb_d4.col(12) += 0 ;
      zb_d4.col(13) += 0 ;
      zb_d4.col(14) += 0 ;
    }
    
  }
  
  zb.reshape(n,1) ;
  
  NumericMatrix out = sumrule(List::create(list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), deriv_order),
                                           derivs_transform(list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), deriv_order), "zeta", {0}, tri, deriv_order)),
                                           tri, deriv_order) ;      
  
  if(!logp){
    out = derivs_transform(out, "exp", {0}, tri, deriv_order) ;//chainrule(List::create(transform(as<arma::mat>(out), "exp", {0}, deriv_order), out), tri, deriv_order) ;
  }
  
  return out ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix pnormexp_cpp (arma::vec q, arma::vec m, arma::vec v, arma::vec u, int s, int deriv_order, List tri, bool logp){
  
  // Initialize
  int n = q.size() ;
  
  arma::vec mq = m - q ;
  arma::vec mq2 = pow(mq,2) ;
  
  arma::vec v2 = arma::pow(v,2) ;
  arma::vec v3 = arma::pow(v,3) ;
  arma::vec v4 = arma::pow(v,4) ;
  
  arma::vec za =  (-m + q)/v ;
  za.reshape(n,1) ;
  arma::mat za_d1 (n, 3)  ;
  arma::mat za_d2 (n, 6)  ;
  arma::mat za_d3 (n, 10)  ;
  arma::mat za_d4 (n, 15)  ;
  
  NumericMatrix zb_NM = dnormexp_cpp(q, m, v, u, s, deriv_order, tri, TRUE) ;
  NumericMatrix zc_NM = transform(u, "log", {0}, deriv_order) ;
  
  arma::vec zb = as<arma::mat>(zb_NM) - as<arma::mat>(zc_NM);
  zb.reshape(n,1) ;
  arma::mat zb_d1 (n, 3)  ;
  arma::mat zb_d2 (n, 6)  ;
  arma::mat zb_d3 (n, 10)  ;
  arma::mat zb_d4 (n, 15)  ;
  
  if(deriv_order>0){
    //First derivatives of za wrt m,v,u
    za_d1.col(0) += -(1/v) ;
    za_d1.col(1) += (mq)/v2 ;
    za_d1.col(2) += 0 ;
    
    //Second derivatives of za wrt mm,mv,mu,vv,vu,uu
    za_d2.col(0) += 0 ;
    za_d2.col(1) += pow(v,-2) ;
    za_d2.col(2) += 0 ;
    za_d2.col(3) += (2*(-m + q))/v3 ;
    za_d2.col(4) += 0 ;
    za_d2.col(5) += 0 ;
    
    //First derivatives of za wrt m,v,u
    zb_d1 += as<arma::mat>(zb_NM.attr("d1")) ;
    zb_d1.col(0) += 0 ;
    zb_d1.col(1) += 0 ;
    zb_d1.col(2) += -as<arma::mat>(zc_NM.attr("d1")) ;
    
    //Second derivatives of za wrt mm,mv,mu,vv,vu,uu
    zb_d2 += as<arma::mat>(zb_NM.attr("d2")) ;
    zb_d2.col(0) += 0 ;
    zb_d2.col(1) += 0 ;
    zb_d2.col(2) += 0 ;
    zb_d2.col(3) += 0 ;
    zb_d2.col(4) += 0 ;
    zb_d2.col(5) += -as<arma::mat>(zc_NM.attr("d2")) ;
    
    if(deriv_order>2){
      //Third derivatives of za wrt mmm,mmv,mmu,mvv,mvu,muu,vvv,vvu,(8)vuu,uuu
      za_d3.col(0) += 0 ;
      za_d3.col(1) += 0 ;
      za_d3.col(2) += 0 ;
      za_d3.col(3) += 0 ;
      za_d3.col(4) += 0 ;
      za_d3.col(5) += 0 ;
      za_d3.col(6) += -2/v3 ;
      za_d3.col(7) += (6*(mq))/v4 ;
      za_d3.col(8) += 0 ;
      za_d3.col(9) += 0 ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,vvuu,vuuu,uuuu
      za_d4.col(0) += 0 ;
      za_d4.col(1) += 0 ;
      za_d4.col(2) += 0 ;
      za_d4.col(3) += 0 ;
      za_d4.col(4) += 0 ;
      za_d4.col(5) += 0 ;
      za_d4.col(6) += 6/v4 ;
      za_d4.col(7) += 0 ;
      za_d4.col(8) += 0 ;
      za_d4.col(9) += 0 ;
      za_d4.col(10) += (24*(-m + q))/pow(v,5) ;
      za_d4.col(11) += 0 ;
      za_d4.col(12) += 0 ;
      za_d4.col(13) += 0 ;
      za_d4.col(14) += 0 ;
      
      //Third derivatives of zb wrt mmm,mmv,mmu,mvv,mvu, (5) muu,vvv,vvu,vuu,uuu
      zb_d3 += as<arma::mat>(zb_NM.attr("d3")) ;
      zb_d3.col(0) += 0 ;
      zb_d3.col(1) += 0 ;
      zb_d3.col(2) += 0 ;
      zb_d3.col(3) += 0 ;
      zb_d3.col(4) += 0 ;
      zb_d3.col(5) += 0 ;
      zb_d3.col(6) += 0 ;
      zb_d3.col(7) += 0 ;
      zb_d3.col(8) += 0 ;
      zb_d3.col(9) += -as<arma::mat>(zc_NM.attr("d3")) ;
      
      //Fourth derivatives of za wrt mmmm,mmmv,mmmu,mmvv,mmvu,mmuu,(6)mvvv,mvvu,mvuu,muuu,(10) vvvv,vvvu,(12) vvuu,vuuu,uuuu
      zb_d4 += as<arma::mat>(zb_NM.attr("d4")) ;
      zb_d4.col(0) += 0 ;
      zb_d4.col(1) += 0 ;
      zb_d4.col(2) += 0 ;
      zb_d4.col(3) += 0 ;
      zb_d4.col(4) += 0 ;
      zb_d4.col(5) += 0 ;
      zb_d4.col(6) += 0 ;
      zb_d4.col(7) += 0 ;
      zb_d4.col(8) += 0 ;
      zb_d4.col(9) += 0 ;
      zb_d4.col(10) += 0 ;
      zb_d4.col(11) += 0 ;
      zb_d4.col(12) += 0 ;
      zb_d4.col(13) += 0 ;
      zb_d4.col(14) += -as<arma::mat>(zc_NM.attr("d4")) ;
    }
  }
  
  za_d1 *= s ;
  za_d2 *= s ;
  za_d3 *= s ;
  za_d4 *= s ;
  
  //NumericMatrix out = sumrule(List::create(derivs_transform(list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), 4), "pnorm", {0}, tri, deriv_order),
  //                                         derivs_transform(list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), 4), "mexp", {0}, tri, deriv_order)),
  //                                         tri, deriv_order) ;      
  
  
  
  
    
  NumericMatrix za_mat = list2derivs(List::create(za, za_d1, za_d2, za_d3, za_d4), 4) ;
  NumericMatrix zb_mat = list2derivs(List::create(zb, zb_d1, zb_d2, zb_d3, zb_d4), 4) ;
  
  NumericMatrix out = sumrule(List::create(derivs_transform(za_mat, "pnorm", {0}, tri, deriv_order),
                                           derivs_transform(zb_mat, "mexp", {0}, tri, deriv_order)),
                                           tri, deriv_order) ;
  
  if(logp){
    out = derivs_transform(out, "log", {0}, tri, deriv_order) ;
  }         
  
  return out ;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix dcomper_cpp (arma::vec x, arma::vec m, arma::vec v, arma::vec u, int s, Rcpp::String distr ,int deriv_order, List tri, bool logp){
  
  Rcpp::StringVector distributions = {"normhnorm","normexp"} ;
  // Initialize
  NumericMatrix out ;
    
    //Normhnorm distr
    if(distr==distributions(0)){
      out = dnormhnorm_cpp (x, m, v, u, s, deriv_order, tri, logp) ;
    }
    
    //Normexp distr
    if(distr==distributions(1)){
      out = dnormexp_cpp (x, m, v, u, s, deriv_order, tri, logp) ;
    }
 
 return out ; 
}


// [[Rcpp::export]]
Rcpp::NumericMatrix pcomper_cpp (arma::vec q, arma::vec m, arma::vec v, arma::vec u, int s, Rcpp::String distr ,int deriv_order, List tri, bool logp){
  
  Rcpp::StringVector distributions = {"normhnorm","normexp"} ;
  // Initialize
  Rcpp::NumericMatrix out ;
  
  //Normhnorm distr
  if(distr==distributions(0)){
    out = pnormhnorm_cpp (q, m, v, u, s, deriv_order, tri, logp) ;
  }
  
  //Normexp distr
  if(distr==distributions(1)){
    out = pnormexp_cpp (q, m, v, u, s, deriv_order, tri, logp) ;
  }
  
  return out ; 
}

// [[Rcpp::export]]
NumericMatrix dcomper_mv_cpp (arma::mat x, arma::mat m, arma::mat v, arma::mat u, arma::vec delta, arma::vec s, Rcpp::StringVector distr ,int rot ,int deriv_order, List tri, bool logp){
  
  //Marginal log pdfs
  Rcpp::NumericMatrix f1 = dcomper_cpp (x.col(0), m.col(0), v.col(0), u.col(0), s(0), distr(0), deriv_order, tri[0], TRUE) ;
  Rcpp::NumericMatrix f2 = dcomper_cpp (x.col(1), m.col(1), v.col(1), u.col(1), s(1), distr(1), deriv_order, tri[1], TRUE) ;
  
  //Marginal cdfs
  Rcpp::NumericMatrix F1 = pcomper_cpp (x.col(0), m.col(0), v.col(0), u.col(0), s(0), distr(0), deriv_order, tri[0], FALSE) ;
  Rcpp::NumericMatrix F2 = pcomper_cpp (x.col(1), m.col(1), v.col(1), u.col(1), s(1), distr(1), deriv_order, tri[1], FALSE) ;
  
  //Cop log pdf
  Rcpp::NumericMatrix cop = dcop_cpp (as<arma::mat>(F1), as<arma::mat>(F2), delta, distr(2), rot, deriv_order, tri[0], TRUE) ;
  
  //Create dummy derivs objects
  int n = as<arma::mat>(f1).n_rows ; 
  arma::mat constant_zeros(n,1,fill::zeros) ; //zeros
  arma::mat constant_ones(n,1,fill::ones) ; //ones
  
  //Combine independent marginal log pdfs and dummy zeroes
  Rcpp::NumericMatrix pdf_extended = ind2joint(List::create(f1, f2, transform(constant_zeros, "constant", {0}, deriv_order)),
                                               List::create(tri[0],tri[1],tri[2]),
                                               List::create(tri[0],tri[3],tri[4]),
                                               deriv_order) ;
  
  //Sum up marginal log pdfs
  Rcpp::NumericMatrix pdf_extended_sum(n, 1, rowSums(pdf_extended).begin()) ;
  pdf_extended_sum.attr("d1") = pdf_extended.attr("d1") ;
  pdf_extended_sum.attr("d2") = pdf_extended.attr("d2") ;
  pdf_extended_sum.attr("d3") = pdf_extended.attr("d3") ;
  pdf_extended_sum.attr("d4") = pdf_extended.attr("d4") ;
  
  //Combine independent marginal cdfs and dummy ones
  Rcpp::NumericMatrix cdf_extended = ind2joint(List::create(F1, F2, transform(constant_ones, "identity", {0}, deriv_order)),
                                               List::create(tri[0],tri[1],tri[2]),
                                               List::create(tri[0],tri[3],tri[4]),
                                               deriv_order) ;
  
  //Extend cop log pdf to fit size of marginal cdfs
  arma::mat cop_d0 = as<arma::mat>(cop) ;
  arma::mat cop_d1 = as<arma::mat>(cop.attr("d1")) ;
  arma::mat cop_d2 = as<arma::mat>(cop.attr("d2")) ;
  
  arma::uvec i1 = {0, 0, 0, 1, 1, 1, 2} ;
  
  arma::uvec i2 = {0, 0, 0, 1, 1, 1, 2,
                      0, 0, 1, 1, 1, 2,
                         0, 1, 1, 1, 2,
                            3, 3, 3, 4,
                               3, 3, 4,
                                  3, 4,
                                     5} ;
    
  arma::mat cop_extended_d1 = cop_d1.cols(i1) ;
  arma::mat cop_extended_d2 = cop_d2.cols(i2) ;
  arma::mat cop_extended_d3(n, 84) ;
  arma::mat cop_extended_d4(n, 210) ;
  
  NumericMatrix cop_extended = list2derivs(List::create(cop_d0 , cop_extended_d1,
                                                                 cop_extended_d2,
                                                                 cop_extended_d3,
                                                                 cop_extended_d4), deriv_order) ;
  
  //Summarize loglike
  NumericMatrix out = sumrule(List::create(chainrule(List::create(cop_extended, cdf_extended), tri[4], deriv_order), 
                                           pdf_extended_sum) ,tri[4], deriv_order) ;
  
  return out ; 
}

/*** R
dcomper_mv_cpp (matrix(c(5,6),ncol=2), matrix(c(5,6),ncol=2), matrix(c(5,6),ncol=2), matrix(c(5,6),ncol=2), c(0.5), c(-1,-1), distr=c("normhnorm","normhnorm","normal") , 2 , tri=list(trind_generator(3),trind_generator(3),trind_generator(1),trind_generator(6),trind_generator(7)), TRUE)
*/
