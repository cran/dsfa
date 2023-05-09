#' Copula function
#'
#' Probablitiy density function, distribution and random number generation for copulas.
#' 
#'
#' @details 
#' A two-dimensional copula \eqn{C(w_1, w_2, \delta)} is a multivariate cumulative distribution function for which the marginal probability distribution of \eqn{w_1} and \eqn{w_1} are uniform on the interval \eqn{[0,1]}.
#' The parameter \eqn{\delta} specifies the copula.
#'
#' The functions \code{pcop()} and \code{rcop()} are wrapper functions for the \code{\link[copula:pCopula]{pCopula()}} and \code{\link[copula:rCopula]{rCopula()}}.
#'
#' @return \code{dcop} gives the density, \code{pcop} gives the distribution function for a specified copula and \code{rcop} generates random numbers, with given \code{delta}.
#' \code{dcop()} returns a \code{derivs} object. For more details see [trind()] and [trind_generator()].
#'
#' @param W numeric matrix of pseudo observations. Must have two columns.
#' @param delta numeric vector of copula parameter \eqn{\delta}.
#' @param distr_cop string, defines the copula family:\cr
#' `independent` = Independence copula \cr
#' `normal` = Gaussian copula \cr
#' `clayton` = Clayton copula \cr
#' `gumbel` = Gumbel copula \cr
#' `frank` = Frank copula \cr
#' `joe` = Joe copula \cr
#' `amh` = Ali-Mikhail-Haq copula \cr
#' @param rot, integer determining the rotation for Archimedian copulas. Can be \code{90}, \code{180} or \code{270}.
#' @inheritParams dcomper
#' 
#' @examples
#' u=0.3; v=0.7; p=0.5
#' pdf <- dcop(W=cbind(u,v), delta=p, distr_cop="normal")
#' cdf <- pcop(W=cbind(u,v), delta=p, distr_cop="normal")
#' r <- rcop(n=100, delta=p, distr_cop="normal")
#'
#' @references
#' \itemize{
#' \item \insertRef{schepsmeier2014derivatives}{dsfa}
#' \item \insertRef{hofert2018elements}{dsfa}
#' }
#' 
#' @family copula
#' 
#' @export
dcop<-function(W, delta, distr_cop="normal", rot=0, deriv_order=0, tri=NULL, log.p=FALSE){
  #Density of copula
  if (any(W <= 0)|any(W >= 1)){
    stop(paste("W must be in [-1, 1]", "\n", ""))
  }
  
  distr_cop<-match.arg(distr_cop,c("independent","normal","clayton","gumbel","frank","joe","amh"))
  
  minmax<-delta_bounds(distr_cop)
  if(any(delta<minmax[1])|any(delta>minmax[2])){
    stop(paste("delta must be in [",minmax[1],",",minmax[2],"]", "\n", ""))
  }
  
  X<-tryCatch(cbind(W, delta), warning=function(w) {
    stop("Input vectors have incompatible lengths")
    })
  
  if(!rot%in%c(0, 90,180,270)){
    stop(paste("rotation must be in {0, 90, 180, 270}", "\n", ""))
  }
  
  if(is.null(tri)){
    tri=trind_generator(3)
  }
  
  out<-dcop_cpp (X[,1], X[,2], X[,3], distr_cop, rot=rot, deriv_order, tri, log.p)
  
  #return out
  return(out)
}

#' @describeIn dcop distribution function for copula.
#' @export
pcop<-function(W, delta=0, distr_cop="normal", rot=0, log.p = FALSE){
  #Distribution function of bivariate copula
  if (any(W <= 0)|any(W >= 1)){
    stop(paste("W must be in [-1, 1]", "\n", ""))
  }
  
  distr_cop<-match.arg(distr_cop,c("independent","normal","clayton","gumbel","frank","joe","amh"))
  
  minmax<-delta_bounds(distr_cop)
  if(any(delta<minmax[1])|any(delta>minmax[2])){
    stop(paste("delta must be in [",minmax[1],",",minmax[2],"]", "\n", ""))
  }
  
  X<-tryCatch(cbind(W, delta), warning=function(w) {
    stop("Input vectors have incompatible lengths")
    })

  if(!rot%in%c(0,90,180,270)){
    stop(paste("rotation must be in {0, 90, 180, 270}", "\n", ""))
  }
  
  out<-sapply(1:nrow(X), function(n) pcop_copula(W=X[n,-3], delta=X[n,3], distr_cop=distr_cop, log.p=log.p))

  #Return out
  return(out)
}

#' @describeIn dcop random number generation for copula.
#' @inheritParams rcomper
#' @export
rcop<-function(n, delta=0, distr_cop="normal", rot=0){
  #Random number generation function of copula
  distr_cop<-match.arg(distr_cop,c("independent","normal","clayton","gumbel","frank","joe","amh"))
  
  minmax<-delta_bounds(distr_cop)
  if(any(delta<minmax[1])|any(delta>minmax[2])){
    stop(paste("delta must be in [",minmax[1],",",minmax[2],"]", "\n", ""))
  }
  
  X<-tryCatch(cbind(rep(0,n), delta), warning=function(w) {
    stop("Input vectors have incompatible lengths")
    })

  if(!rot%in%c(0,90,180,270)){
    stop(paste("rotation must be in {0, 90, 180, 270}", "\n", ""))
  }
  
  N<-n
  out<-lapply(1:N, function(i) rcop_copula(delta=X[i,2, drop=T], distr_cop=distr_cop))
  out<-matrix(unlist(out), byrow=TRUE, nrow=N)

  return(out)
}

pcop_copula<-function(W, delta=0, distr_cop="normal", rot=0, log.p = FALSE){
  #dcop wrapper function for distribution function for copula with scalar inputs

  if(distr_cop=="independent"){
    cop_object<-copula::indepCopula(param=0, dim = 2)
  }

  if(distr_cop=="normal"){
    #cop_dim<-nrow(copula::p2P(c(delta)))
    cop_object<-copula::normalCopula(param=delta, dim = 2)
  }

  if(distr_cop=="clayton"){
    cop_object<-copula::claytonCopula(param=delta, dim = 2)
  }

  if(distr_cop=="gumbel"){
    cop_object<-copula::gumbelCopula(param=delta, dim = 2)
  }

  if(distr_cop=="frank"){
    cop_object<-copula::frankCopula(param=delta, dim = 2)
  }

  if(distr_cop=="joe"){
    cop_object<-copula::joeCopula(param=delta, dim = 2)
  }

  if(distr_cop=="amh"){
    cop_object<-copula::amhCopula(param=delta, dim = 2)
  }
  
  if(rot==90){
    cop_object<-copula::rotCopula(cop_object, flip = c(TRUE, FALSE))
  }
  
  if(rot==180){
    cop_object<-copula::rotCopula(cop_object, flip = c(TRUE, TRUE))
  }
  
  if(rot==270){
    cop_object<-copula::rotCopula(cop_object, flip = c(FALSE, TRUE))
  }
  
  
  #Wrapper for Copula package function
  out<-copula::pCopula(u=W, copula=cop_object)

  if(log.p){
    out<-log(out)
  }

  return(out)
}


rcop_copula<-function(delta=0, distr_cop="normal", rot=0){
  #wrapper function for random number generation for copula with scalar inputs

  if(distr_cop=="independent"){
    cop_object<-copula::indepCopula(param=0, dim = 2)
  }

  if(distr_cop=="normal"){
    #cop_dim<-nrow(copula::p2P(c(delta)))
    cop_object<-copula::normalCopula(param=delta, dim = 2)
  }

  if(distr_cop=="clayton"){
    cop_object<-copula::claytonCopula(param=delta, dim = 2)
  }

  if(distr_cop=="gumbel"){
    cop_object<-copula::gumbelCopula(param=delta, dim = 2)
  }

  if(distr_cop=="frank"){
    cop_object<-copula::frankCopula(param=delta, dim = 2)
  }

  if(distr_cop=="joe"){
    cop_object<-copula::joeCopula(param=delta, dim = 2)
  }

  if(distr_cop=="amh"){
    cop_object<-copula::amhCopula(param=delta, dim = 2)
  }
  
  if(rot==90){
    cop_object<-copula::rotCopula(cop_object, flip = c(TRUE, FALSE))
  }
  
  if(rot==180){
    cop_object<-copula::rotCopula(cop_object, flip = c(TRUE, TRUE))
  }
  
  if(rot==270){
    cop_object<-copula::rotCopula(cop_object, flip = c(FALSE, TRUE))
  }
  
  #Wrapper for Copula package function
  out<-copula::rCopula(n=1, copula=cop_object)

  return(out)
}
