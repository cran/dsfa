#' Utility functions
#'
#' Collection of utility functions.
#'
#' @keywords internal
#'
#' @param H symmetric matrix of full rank.
#'
#' @details Collection of utility functions which are for internal usage. \code{pos_eigen} takes a symmetric matrix \code{H} and returns matrix \code{I} with positive eigenvalues as well as square root and inverse square root
#'
#' @return \code{pos_eigen} returns a list with three elements: 1) Matrix with positive eigenvalues, 2) Matrix with square root of positive eigenvalues and 3) Matrix with inverse positive eigenvalues.
#' \code{remove_attr} returns the input object without attributes.
#' @examples
#' H<-diag(3)*(-1)
#' pos_eigen(H)
#'
#' x<-1
#' attr(x,"gradient")<-2
#' x<-remove_attr(x)
#'
#' @export
#' @keywords internal
pos_eigen<-function(H){
  #Takes a symmetric matrix H and returns matrix I with positive eigenvalues as well as square root and inverse square root
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

#' @describeIn pos_eigen Remove attribute function
#' @param object vector or matrix.
#' @export
#' @keywords internal
remove_attr<-function(object){
  unlist(lapply(object, function(x) { attributes(x) <- NULL; x }))
}

