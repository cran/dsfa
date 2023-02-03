#' Normal-Exponential distribution
#'
#' Probablitiy density function, distribution, quantile function and random number generation for the normal-exponential distribution
#'
#' @return \code{dnormexp()} gives the density, \code{pnormexp()} give the distribution function, \code{qnormexp()} gives the quantile function, and \code{rnormexp()} generates random numbers, with given parameters.
#' \code{dnormexp()} and \code{pnormexp()} return a \code{derivs} object. For more details see [trind()] and [trind_generator()].
#'
#' @details A random variable \eqn{X} follows a normal-exponential distribution if \eqn{X = V + s \cdot U }, where \eqn{V \sim N(\mu, \sigma_V^2)} and \eqn{U \sim Exp(\sigma_u)}.
#' The density is given by \deqn{f_X(x)=\frac{\sigma_u}{2} \exp \{\sigma_u (s \mu) + \frac{1}{2} \sigma_u^2 \sigma_V^2-\sigma_u (s x) \} 2 \Phi(\frac{1}{\sigma_V} (-s \mu)-\sigma_u \sigma_V+\frac{1}{\sigma_V}(s x)) \qquad,}
#' where \eqn{s=-1} for production and \eqn{s=1} for cost function. In the latter case the distribution is equivalent to the Exponentially modified Gaussian distribution. 
#'
#' @aliases normexp
#'
#' @inheritParams dcomper
#'
#' @examples
#' pdf <- dnormexp(x=5, mu=1, sigma_v=2, sigma_u=3, s=-1)
#' cdf <- pnormexp(q=5, mu=1, sigma_v=2, sigma_u=3, s=-1)
#' q <- qnormexp(p=seq(0.1, 0.9, by=0.1), mu=1, sigma_v=2, sigma_u=3, s=-1)
#' r <- rnormexp(n=10, mu=1, sigma_v=2, sigma_u=3, s=-1)
#'
#' @references
#' \itemize{
#' \item \insertRef{aigner1977formulation}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{schmidt2020analytic}{dsfa}
#' \item \insertRef{gradshteyn2014table}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' }
#' 
#' @family distribution
#' 
#' @export
dnormexp <- function(x, mu=0, sigma_v=1, sigma_u=1, s=-1, deriv_order=0, tri=NULL, log.p = FALSE){
  #Density function of the normexp distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(x, mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  if(is.null(tri)){
    tri=trind_generator(3)
  }
  
  out<-dnormexp_cpp(x=X[,1, drop=T], m=X[,2, drop=T], v=X[,3, drop=T], u=X[,4, drop=T], s=s, deriv_order=deriv_order, tri=tri, logp=log.p)
  
  #Return ouptut
  return(out)
}

#' @describeIn dnormexp distribution function for the normal-exponential distribution.
#' @inheritParams pcomper
#' @export
pnormexp <- function(q, mu=0, sigma_v=1, sigma_u=1, s=-1, deriv_order=0, tri=NULL, log.p = FALSE){
  #Probability function of the normexp distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(q, mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  if(is.null(tri)){
    tri=trind_generator(3)
  }
  
  out<-pnormexp_cpp(q=X[,1, drop=T], m=X[,2, drop=T], v=X[,3, drop=T], u=X[,4, drop=T], s=s, deriv_order=deriv_order, tri=tri, logp=log.p)
  
  #Return ouptut
  return(out)
}

#' @describeIn dnormexp quantile function for the normal-exponential distribution.
#' @inheritParams qcomper
#' @export
qnormexp <- function(p, mu=0, sigma_v=1, sigma_u=1, s=-1, log.p = FALSE){
  #Quantile function of the normexp distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(p, mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  out<-sapply(1:nrow(X), function(i) cdf2quantile(X[i,1, drop=F], pnormexp, mu=X[i,2, drop=T], sigma_v=X[i,3, drop=T], sigma_u=X[i,4, drop=T], s=s, interval=c(-3*(X[i,1, drop=T]*sqrt(X[i,2, drop=T]^2+1/X[i,3, drop=T]^2)+X[i,2, drop=T]+1/X[i,3, drop=T]),3*(X[i,1, drop=T]*sqrt(X[i,2, drop=T]^2+1/X[i,3, drop=T]^2)+X[i,2, drop=T]+1/X[i,3, drop=T])), tol=.Machine$double.eps))
  
  #Return output
  return(out)
}

#' @describeIn dnormexp random number generation for the normal-exponential distribution.
#' @inheritParams rcomper
#' @export
rnormexp <- function(n, mu=0, sigma_v=1, sigma_u=1, s=-1){
  #Function to generate n random numbers from the normexp distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(rep(0,n), mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  out<-stats::rnorm(n, mean=X[,2], sd=X[,3])+s*stats::rexp(n, X[,4])
  
  #Return output
  return(out)
}