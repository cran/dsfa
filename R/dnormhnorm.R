#' Normal-halfnormal distribution
#'
#' Probablitiy density function, distribution, quantile function and random number generation for the normal-halfnormal distribution
#'
#' @return \code{dnormhnorm()} gives the density, \code{pnormhnorm()} give the distribution function, \code{qnormhnorm()} gives the quantile function, and \code{rnormhnorm()} generates random numbers, with given parameters.
#' \code{dnormhnorm()} and \code{pnormhnorm()} return a \code{derivs} object. For more details see [trind()] and [trind_generator()].
#'
#' @details A random variable \eqn{X} follows a normal-halfnormal distribution if \eqn{X = V + s \cdot U }, where \eqn{V \sim N(\mu, \sigma_V^2)} and \eqn{U \sim HN(\sigma_U^2)}.
#' The density is given by \deqn{f_X(x)=\frac{1}{\sqrt{\sigma_V^2+\sigma_U^2}} \phi(\frac{x-\mu}{\sqrt{\sigma_V^2+\sigma_U^2}}) \Phi(s \frac{\sigma_U}{\sigma_V} \frac{x-\mu}{\sqrt{\sigma_V^2+\sigma_U^2}}) \qquad,}
#' where \eqn{s=-1} for production and \eqn{s=1} for cost function.
#'
#' @aliases normhnorm
#'
#' @inheritParams dcomper
#'
#' @examples
#' pdf <- dnormhnorm(x=5, mu=1, sigma_v=2, sigma_u=3, s=-1)
#' cdf <- pnormhnorm(q=5, mu=1, sigma_v=2, sigma_u=3, s=-1)
#' q <- qnormhnorm(p=seq(0.1, 0.9, by=0.1), mu=1, sigma_v=2, sigma_u=3, s=-1)
#' r <- rnormhnorm(n=10, mu=1, sigma_v=2, sigma_u=3, s=-1)
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
dnormhnorm <- function(x, mu=0, sigma_v=1, sigma_u=1, s=-1, deriv_order=0, tri=NULL, log.p = FALSE){
  #Density function of the normhnorm distribution
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
  
  out<-dnormhnorm_cpp(x=X[,1, drop=T], m=X[,2, drop=T], v=X[,3, drop=T], u=X[,4, drop=T], s=s, deriv_order=deriv_order, tri=tri, logp=log.p)

  #Return ouptut
  return(out)
}

#' @describeIn dnormhnorm distribution function for the normal-halfnormal distribution.
#' @inheritParams pcomper
#' @export
pnormhnorm <- function(q, mu=0, sigma_v=1, sigma_u=1, s=-1, deriv_order=0, tri=NULL, log.p = FALSE){
  #Probability function of the normhnorm distribution
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
  
  out<-pnormhnorm_cpp(q=X[,1, drop=T], m=X[,2, drop=T], v=X[,3, drop=T], u=X[,4, drop=T], s=s, deriv_order=deriv_order, tri=tri, logp=log.p)
  
  #Return ouptut
  return(out)
}

#' @describeIn dnormhnorm quantile function for the normal-halfnormal distribution.
#' @inheritParams qcomper
#' @export
qnormhnorm <- function(p, mu=0, sigma_v=1, sigma_u=1, s=-1){
  #Quantile function of the normhnorm distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(p, mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  out<-sapply(1:nrow(X), function(i) cdf2quantile(X[i,1, drop=F], pnormhnorm, mu=X[i,2, drop=T], sigma_v=X[i,3, drop=T], sigma_u=X[i,4, drop=T], s=s, interval=c(-3*(X[i,1, drop=T]*sqrt(X[i,2, drop=T]^2+X[i,3, drop=T]^2)+X[i,2, drop=T]),3*(X[i,1, drop=T]*sqrt(X[i,2, drop=T]^2+X[i,3, drop=T]^2)+X[i,2, drop=T])), tol=.Machine$double.eps))
  
  #Return output
  return(out)
}

#' @describeIn dnormhnorm random number generation for the normal-halfnormal distribution.
#' @inheritParams rcomper
#' @export
rnormhnorm <- function(n, mu=0, sigma_v=1, sigma_u=1, s=-1){
  #Function to generate n random numbers from the normhnorm distribution
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  X<-tryCatch(cbind(rep(0,n), mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  omega<-sqrt(X[,3]^2+X[,4]^2)
  alpha<-s*X[,4]/X[,3]
  delta<-alpha/sqrt(1+alpha^2)
  
  tn <- matrix(stats::rnorm(2 * n), 2, n, byrow = FALSE)
  chi <- c(abs(tn[1, ]))
  nrv <- c(tn[2, ])
  z <- delta * chi + sqrt(1 - delta^2) * nrv
  out <- as.vector(X[,2] + omega * z)

  #Return output
  return(out)
}