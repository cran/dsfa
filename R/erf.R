#' erf function
#'
#' Error function, complementary error function, inverse error function,inverse complementary error function and first and second derivative of quantile function of the standard normal distribution.
#'
#' @param x vector of quantiles.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param deriv derivative of order \code{deriv}. Available are \code{0} and \code{2}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}  otherwise, \eqn{P[X > x]}.
#'
#' @details \code{erf} is the error function. \code{erfc} is the complementary error function.
#' \code{erfinv} is the inverse error function and and \code{erfcinv} is the inverse complementary error function. \code{pnorm} and \code{qnorm} are clones from
#' the \code{stats} package, but do have the argument for derivatives.
#'
#' @return Returns a vector of the corresponding function evaluated at x.
#'
#' @examples
#' erf(0.5)
#' erfc(0.5)
#' erfinv(0.5)
#' erfcinv(0.5)
#'
#' qnorm(0.5, deriv=2)
#' pnorm(0.5, deriv=2)
#'
#' @export
#' erf
erf <- function(x){
  2 * stats::pnorm(x * sqrt(2)) - 1
}

#' @describeIn erf complementary error function
#' @export
erfc <- function(x){
  2 * stats::pnorm(x * sqrt(2), lower = FALSE)
}

#' @describeIn erf inverse error function
#' @export
erfinv <- function (x){
  stats::qnorm((1 + x)/2)/sqrt(2)
}

#' @describeIn erf inverse complementary error function
#' @export
erfcinv <- function (x) {
  stats::qnorm(x/2, lower = FALSE)/sqrt(2)
}

#' @describeIn erf Distribution function of the standard normal distribution with derivatives.
#' @param q vector of quantiles.
#' @export
pnorm<-function(q, mean=0, sd=1, deriv=0, log.p=FALSE, lower.tail=TRUE){
  z<-(q-mean)/sd
  out<-stats::pnorm(q=z, lower.tail = lower.tail, log.p=log.p)
  if(deriv>0){
    attr(out,"gradient")<-stats::dnorm(z)
    attr(out,"hessian")<--z*stats::dnorm(z)
  }

  return(out)
}

#' @describeIn erf Quantile function of the standard normal distribution with derivatives.
#' @param p vector of probabilites.
#' @export
qnorm<-function(p, mean=0, sd=1, deriv=0, log.p=FALSE, lower.tail=TRUE){
  z<-(p-mean)/sd
  out<-stats::qnorm(p=z, lower.tail=lower.tail, log.p=log.p)
  if(deriv>0){
    attr(out,"gradient")<-sqrt(2*pi)/exp(-stats::qnorm(z)^2/2)
    attr(out,"hessian")<-stats::qnorm(z)*attr(out,"gradient")^2
  }

  return(out)
}


