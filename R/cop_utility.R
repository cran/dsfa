#' Utility functions
#'
#' Internal utility functions for copulae
#'
#' @param U matrix of pseudo observations. Must have two columns.
#' @param Tau matrix of Kendall's tau.
#' @param family integer, defines the copula family:\cr
#' `1` = Gaussian copula \cr
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @return \code{dcop_copula} gives the density, \code{pcop_copula} gives the distribution function for a specified copula and \code{rcop_copula} generates random numbers, with given Tau.
#' These functions are not written for vector inputs. \code{dcop_copula_deriv} provides a wrapper function for the 'numDeriv' package functions, such that the gradient and hessian of the density can be evaluated.
#'
#' @export
#' @keywords internal
#'density of copula
dcop_copula<-function(U, Tau=0, family=1, log.p = FALSE){

  if(family==1){
    #Gaussian copula
    cop_dim<-nrow(copula::p2P(c(Tau)))
    cop_object<-copula::normalCopula(param=sin(Tau*pi/2), dim = cop_dim , dispstr = "un")
    #Wrapper for Copula package function
    out<-copula::dCopula(u=U, copula=cop_object)
  }

  if(family==3){
    #Clayton copula
    cop_object<-copula::claytonCopula(param=2*Tau/(1-Tau), dim = 2)
    #Wrapper for Copula package function
    out<-copula::dCopula(u=U, copula=cop_object)
  }

  if(log.p){
    out<-log(out)
  }

  return(out)
}

#' @describeIn dcop_copula Wrapper function of dcop_copula for usage with 'numDeriv' package
#' @export
#' @keywords internal
dcop_copula_deriv<-function(x, U=NULL, family=1, D=2){
  #Wrapper function for numDeriv
  if(!is.null(U)){
    #Disjoint
    out<-dcop_copula(U=U, Tau=x, family=family, log.p = TRUE)
  } else {
    #Joint
    out<-dcop_copula(U=x[,c(1:D)], Tau=x[,-c(1:D)], family=family, log.p = TRUE)
  }

  return(out)
}

#' @describeIn dcop_copula distribution function of copula
#' @export
#' @keywords internal
pcop_copula<-function(U, Tau=0, family=1, log.p = FALSE){

  if(family==1){
    #Gaussian copula
    cop_dim<-nrow(copula::p2P(c(Tau)))
    cop_object<-copula::normalCopula(param=sin(Tau*pi/2), dim = cop_dim , dispstr = "un")
    #Wrapper for Copula package function
    out<-copula::pCopula(u=U, copula=cop_object)
  }

  if(family==3){
    #Clayton copula
    cop_object<-copula::claytonCopula(param=2*Tau/(1-Tau), dim = 2)
    #Wrapper for Copula package function
    out<-copula::pCopula(u=U, copula=cop_object)
  }

  if(log.p){
    out<-log(out)
  }
  return(out)
}

#' @describeIn dcop_copula random number generation for copula
#' @export
#' @keywords internal
rcop_copula<-function(Tau=0, family=1){

  if(family==1){
    #Gaussian copula
    cop_dim<-nrow(copula::p2P(c(Tau)))
    cop_object<-copula::normalCopula(param=sin(Tau*pi/2), dim = cop_dim)# , dispstr = "un")
    #Wrapper for Copula package function
    out<-copula::rCopula(1, copula=cop_object)
  }

  if(family==3){
    #Clayton copula
    cop_object<-copula::claytonCopula(param=2*Tau/(1-Tau), dim = 2)
    #Wrapper for Copula package function
    out<-copula::rCopula(n=1, copula=cop_object)
  }

  return(out)
}
