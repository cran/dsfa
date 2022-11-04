#' Zeta
#'
#' Evaluates the zeta function, i.e. \code{log(2*pnorm(x))}.
#'
#' @param x numeric vector of quantiles. Missing values (NAs) and Inf are allowed.
#' @param deriv derivative of order \code{deriv} of the zeta function. Available are \code{0}, \code{2} and \code{4}.
#'
#' @return Evaluates the zeta function.  If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian}, \code{l3} and \code{l4}.
#'
#' @details The zeta function is defined as \eqn{log\{ 2 \Phi(x)\}}. The function is a clone of \code{\link[sn:zeta]{zeta()}}.
#'
#'
#' @examples
#' zeta(1,1)
#'
#' @references
#' \itemize{
#' \item \insertRef{abramowitz1964handbook}{dsfa}
#' \item \insertRef{azzalini2014collaboration}{dsfa}
#' }
#' @export
#'
#zeta
zeta<-function (x, deriv=0) {
  na <- is.na(x)
  x <- replace(x, na, 0)
  x2 <- x^2
  out<-stats::pnorm(x, log.p = TRUE) + log(2)

  if(deriv>0){
    gradient<-ifelse(x > (-50), exp(stats::dnorm(x, log = TRUE) - stats::pnorm(x, log.p = TRUE)),
             -x/(1 - 1/(x2 + 2) + 1/((x2 + 2) * (x2 + 4)) - 5/((x2 + 2) * (x2 + 4) * (x2 + 6)) + 9/((x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8)) - 129/((x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8) * (x2 + 10))))
    hessian<-(-gradient * (x + gradient))

    attr(out,"gradient")<-gradient
    attr(out,"hessian")<-hessian

    if(deriv>2){
        l3<-(-hessian * (x + gradient) - gradient * (1 + hessian))
        l4<-(-l3 * (x + 2 * gradient) - 2 * hessian * (1 + hessian))
        attr(out,"l3")<-l3
        attr(out,"l4")<-l4
    }
  }

  neg.inf <- (x == -Inf)
  if (any(neg.inf)){
    if(deriv>0){
      attr(out,"gradient")<-replace(attr(out,"gradient"), neg.inf, Inf)
      attr(out,"hessian")<-replace(attr(out,"hessian"), neg.inf, -1)
      if(deriv>2){
        attr(out,"l3")<-replace(attr(out,"l3"), neg.inf, 0)
        attr(out,"l4")<-replace(attr(out,"l4"), neg.inf, 0)
      }
    }
  }

  if(deriv>0){
    attr(out,"gradient")<-replace(attr(out,"gradient"), x == Inf, 0)
    attr(out,"hessian")<-replace(attr(out,"hessian"), x == Inf, 0)
    if(deriv>2){
      attr(out,"l3")<-replace(attr(out,"l3"), x == Inf, 0)
      attr(out,"l4")<-replace(attr(out,"l4"), x == Inf, 0)
    }
  }

  replace(out, na, NA)

  #Return
  return(out)
}
