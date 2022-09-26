#' OwenT
#'
#' Evaluates the Owen T-function
#'
#' @param x numeric vector of quantiles. Missing values (NAs) and Inf are allowed.
#' @param a numeric vector. Inf is allowed.
#' @param jmax an integer scalar value which regulates the accuracy of the result.
#' @param cut.point a scalar value which regulates the behaviour of the algorithm.
#' @param deriv derivative of order \code{deriv} of the T.Owen function. Available are \code{0} and \code{2}.
#'
#' @return \code{OwenT} evaluates the OwenT function with given parameters \code{x} and \code{a}.  If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian}.
#'
#' @details The OwenT function is defined as \deqn{T(x,a)=\frac{1}{2 \pi} \int_0^{a} \frac{\exp\{-x^2 (1+t^2)/2 \}}{1+t^2} d t}. If \code{a}>1 and 0<\code{x}<=\code{cut.point}, a series expansion is used, truncated after \code{jmax} terms. If \code{a}>1 and \code{x}>\code{cut.point},
#'  an asymptotic approximation is used. In the other cases, various reflection properties of the function are exploited. For \code{deriv}=0, the function is a clone of \code{\link[sn:T.Owen]{T.Owen}}.
#'
#'
#' @examples
#' OwenT(x=1, a=1, jmax = 50, cut.point = 8, deriv=0)
#'
#' @references
#' \itemize{
#' \item \insertRef{Owen1956TablesFC}{dsfa}
#' }
#' @export
#'
#OwenT
OwenT<-base::Vectorize(function (x, a, jmax = 50, cut.point = 8, deriv=0){
  T.int <- function(x, a, jmax, cut.point) {
    fui <- function(x, i) (x^(2 * i))/((2^i) * gamma(i +
                                                       1))
    seriesL <- seriesH <- NULL
    i <- 0:jmax
    low <- (x <= cut.point)
    hL <- x[low]
    hH <- x[!low]
    L <- length(hL)
    if (L > 0) {
      b <- outer(hL, i, fui)
      cumb <- apply(b, 1, cumsum)
      b1 <- exp(-0.5 * hL^2) * t(cumb)
      matr <- matrix(1, jmax + 1, L) - t(b1)
      jk <- rep(c(1, -1), jmax)[1:(jmax + 1)]/(2 * i +
                                                 1)
      matr <- t(matr * jk) %*% a^(2 * i + 1)
      seriesL <- (atan(a) - as.vector(matr))/(2 * pi)
    }
    if (length(hH) > 0)
      seriesH <- atan(a) * exp(-0.5 * (hH^2) * a/atan(a)) *
      (1 + 0.00868 * (hH * a)^4)/(2 * pi)
    series <- c(seriesL, seriesH)
    id <- c((1:length(x))[low], (1:length(x))[!low])
    series[id] <- series
    series
  }
  if (!is.vector(a) | length(a) > 1)
    stop("'a' must be a vector of length 1")
  if (!is.vector(x))
    stop("'x' must be a vector")
  aa <- abs(a)
  ah <- abs(x)
  if (is.na(aa))
    stop("parameter 'a' is NA")
  if (aa == Inf)
    return(sign(a) * 0.5 * pnorm(-ah))
  if (aa == 0)
    return(rep(0, length(x)))
  na <- is.na(x)
  inf <- (ah == Inf)
  ah <- replace(ah, (na | inf), 0)
  if (aa <= 1)
    owen <- T.int(ah, aa, jmax, cut.point)
  else owen <- (0.5 * pnorm(ah) + pnorm(aa * ah) * (0.5 - pnorm(ah)) -
                  T.int(aa * ah, (1/aa), jmax, cut.point))
  owen <- replace(owen, na, NA)
  owen <- replace(owen, inf, 0)
  out<-owen * sign(a)

  if(deriv>0){
    gradient<-matrix(0,ncol=2,nrow=length(x))
    hessian<-matrix(0,ncol=4,nrow=length(x))

    attr(out,"gradient")<-gradient
    attr(out,"hessian")<-hessian
  }

  #Return out
  return(out)
}, vectorize.args=c("x","a"))


