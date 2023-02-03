#' Inverse cumulative distribution function
#' 
#' @inheritParams dcomper
#' 
#' @param cdf function, cumulative distribution function which to invert.
#' @param interval numeric vector of length 2, determining the lower and upper bound of the uniroot interval
#' @param ... other arguments for the cdf, e.g. mu, sigma_v, sigma_u, s...
#' @return Numeric vector of \code{p} evaluated in the inverse cdf.
#' 
#' @details Code is a clone from the package 'gbutils'.
#' 
#' @examples 
#' q=5
#' cdf <- pnorm(q=q, mean=1, sd=2)
#' q_numeric <- cdf2quantile(p=cdf, cdf=pnorm, mean=1, sd=2)
#' all.equal(q,q_numeric)
#' 
#' @export
cdf2quantile<-function (p, cdf, interval = c(-3, 3),...){
  f <- function(x, ...) {
    cdf(x, ...) - p
  }
  
  out <- stats::uniroot(f, lower = min(interval), upper = max(interval), extendInt = "upX",...)$root
  
  return(out)
}
