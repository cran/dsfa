#' Utility functions
#'
#' Collection of utility functions.
#'
#' @keywords internal
#'
#' @aliases erfcinv, d1qnormdx1, d2qnormdx2, vec_T.Owen, vec_qsn, pos_eigen, remove_attr
#' @param x vector of quantiles.
#'
#' @details Collection of utility functions which are for internal usage.
#'
#' @return erfcinv returns a vector of the inverse complementary error function evaluated at x.
#' @examples
#' erfcinv(0.5)
#' d1qnormdx1(0.5)
#' d2qnormdx2(0.5)
#' vec_T.Owen(c(1,1),c(1,2))
#' vec_qsn(c(0.5,0.5))
#' x<-1
#' attr(x,"gradient")<-2
#' x<-remove_attr(x)
#'
#' @export
erfcinv <- function (x) {
  stats::qnorm(x/2, lower = FALSE)/sqrt(2)
}

#' @describeIn erfcinv error function
#' @export
#' @keywords internal
#' @return erf returns a vector of the error function evaluated at x.
erf <- function(x){
  2 * stats::pnorm(x * sqrt(2)) - 1
}

#' @describeIn erfcinv complementary error function
#' @export
#' @keywords internal
#' @return erfc returns a vector of complementary error function evaluated at x.
erfc <- function(x){
  2 * stats::pnorm(x * sqrt(2), lower = FALSE)
}

#' @describeIn erfcinv inverse error function
#' @export
#' @keywords internal
#' @return erfinv returns a vector of the inverse error function evaluated at x.
erfinv <- function (x){
  stats::qnorm((1 + x)/2)/sqrt(2)
}


#' @describeIn erfcinv First derivative of quantile function of the standard normal distribution
#' @export
#' @keywords internal
#' @return d1qnormdx1 returns a vector of first derivative of quantile function of the standard normal distribution evaluated at x.
d1qnormdx1<-function(x){
  sqrt(2*pi)/exp(-stats::qnorm(x)^2/2)
}

#' @describeIn erfcinv Second derivative of quantile function of the standard normal distribution
#' @export
#' @keywords internal
#' @return d2qnormdx2 returns a vector of second derivative of quantile function of the standard normal distribution evaluated at x.
d2qnormdx2<-function(x){
  stats::qnorm(x)*d1qnormdx1(x)^2
}

#' @describeIn erfcinv Vectorized T.Owen function T(h,a)
#' @param h numeric vector, for more details see \code{\link[sn:T.Owen]{T.Owen}}.
#' @param a numeric vector, for more details see \code{\link[sn:T.Owen]{T.Owen}}.
#' @export
#' @keywords internal
#' @return vec_T.Owen returns a vector of the T.Owen function evaluated at h and a.
vec_T.Owen<-base::Vectorize(sn::T.Owen)

#' @describeIn erfcinv Vectorized quantile function of the skew normal distribution
#' @param p numeric vector.
#' @param xi numeric vector, for more details see \code{\link[sn:qsn]{qsn}}.
#' @param omega numeric vector, for more details see \code{\link[sn:qsn]{qsn}}.
#' @param alpha numeric vector, for more details see \code{\link[sn:qsn]{qsn}}.
#' @export
#' @keywords internal
#' @return vec_qsn returns a vector of quantile function of the skew-normal distribution.
vec_qsn<-base::Vectorize(sn::qsn)


#' @describeIn erfcinv Positive eigenvalue function
#' @export
#' @param H symmetric matrix of full rank.
#' @keywords internal
#' @return pos_eigen returns a list with three elements: 1) Matrix with positive eigenvalues, 2) Matrix with square root of positive eigenvalues and 3) Matrix with inverse positive eigenvalues.
pos_eigen<-function(H){
  #Takes a symmetric matrix H and returns matrix I with pos eigen values as well as squareroot and inverse squareroot
  I<-H

  tolI <- sqrt(.Machine$double.eps)
  I.eig <- eigen(I, symmetric = TRUE)
  if (min(I.eig$values) < tolI && sign(min(sign(I.eig$values))) ==  -1)    I.eig$values <- abs(I.eig$values)
  if (min(I.eig$values) < tolI) {
    pep <- which(I.eig$values < tolI)
    I.eig$values[pep] <- tolI
  }

  #Calculate squareroot of eigenvalues
  eig.sr <- sqrt(I.eig$val)


  #Calculate squareroot of I
  I.sr <- I.eig$vec %*% tcrossprod(diag(eig.sr, nrow = length(eig.sr),
                                        ncol = length(eig.sr)), I.eig$vec)

  #Calculate inverse squareroot of I
  I.invsr <- I.eig$vec %*% tcrossprod(diag(1/eig.sr, nrow = length(eig.sr),
                                           ncol = length(eig.sr)), I.eig$vec)

  I <- I.eig$vec %*% tcrossprod(diag(I.eig$val, nrow = length(I.eig$val),
                                     ncol = length(I.eig$val)), I.eig$vec)

  out<-list(I=I, I.sr=I.sr, I.invsr=I.invsr)
  return(out)
}

#' @describeIn erfcinv Remove attribute function
#' @param object vector or matrix.
#' @export
#' @keywords internal
#' @return remove_attr returns the input object without attributes.
remove_attr<-function(object){
  unlist(lapply(object, function(x) { attributes(x) <- NULL; x }))
}

