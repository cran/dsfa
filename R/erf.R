#' erf function
#'
#' Error function, complementary error function, inverse error function,inverse complementary error function and first and second derivative of quantile function of the standard normal distribution.
#'
#' @param x vector of quantiles.
#'
#' @details \code{erf} is the error function. \code{erfc} is the complementary error function.
#' \code{erfinv} is the inverse error function and and \code{erfcinv} is the inverse complementary error function.
#' \code{d1qnormdx1} evaluates the first derivative of quantile function of the standard normal distribution w.r.t. \code{x}.
#' \code{d2qnormdx2} evaluates the second derivative of quantile function of the standard normal distribution w.r.t. \code{x}.
#'
#' @return Returns a vector of the corresponding function evaluated at x.
#'
#' @examples
#' erf(0.5)
#' erfc(0.5)
#' erfinv(0.5)
#' erfcinv(0.5)
#' d1qnormdx1(0.5)
#' d2qnormdx2(0.5)
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

#' @describeIn erf First derivative of quantile function of the standard normal distribution
#' @export
#' @keywords internal
d1qnormdx1<-function(x){
  sqrt(2*pi)/exp(-stats::qnorm(x)^2/2)
}

#' @describeIn erf Second derivative of quantile function of the standard normal distribution
#' @export
#' @keywords internal
d2qnormdx2<-function(x){
  stats::qnorm(x)*d1qnormdx1(x)^2
}

