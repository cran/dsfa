#' Check inputs arguments
#'
#' Checks for valid inputs of the probablitiy density, distribution, quantile function and random number generation for the normal-halfnormal and normal-exponential distribution.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param q vector of quantiles.
#' @param n number of observations.
#' @param mu vector of \eqn{\mu}
#' @param sigma_v vector of \eqn{\sigma_V}. Must be positive.
#' @param par_u vector of parameter of the (in)efficiency term. Must be positive.
#' @param s \eqn{s=-1} for production and \eqn{s=1} for cost function.
#' @param family \eqn{normhnorm} for normal-halfnormal and \eqn{normexp} for normal-exponential distribution.
#' @param deriv derivative of order \code{deriv} of the log density. Available are \code{0},\code{2} and \code{4}.
#'
#' @details Mostly internal function. This function is written to automatically check the inputs and if required to recycle the inputs such that they have the same length. The length is determined by the longest input.
#'
#' @return Returns a list with five elements:\cr
#' `x` = x,p,q of the recycled length. If the number of observations are provided returns numeric `n` \cr
#' `mu` = mu of the recycled length  \cr
#' `sigma_v` = sigma_v of the recycled length \cr
#' `par_u` = par_u of the recycled length \cr
#' `recycle` = logical; if TRUE, the numerical arguments have been recycled to the length of the longest input. \cr
#' @examples
#'
#' check_arguments(x=c(-1,1), mu=1, sigma_v=2, par_u=3, s=-1, family="normhnorm", deriv=2)
#' @export
#' @keywords internal
check_arguments<-function(x=NULL, q=NULL, p=NULL, n=NULL, mu=NULL, sigma_v=NULL, par_u=NULL, s=NULL, family=NULL, deriv=NULL){
  #Check for correct inputs
  if (!is.null(x)){
    if (!is.numeric(x)){
      stop(paste("x must be numeric", "\n", ""))
    }
  }

  if (!is.null(q)){
    if (!is.numeric(q)){
      stop(paste("q must be numeric", "\n", ""))
    }
    x<-q
  }
  if (!is.null(p)){
    if (!is.numeric(p)|any(0>p)|any(p>1)){
      stop(paste("p must be in [0,1]", "\n", ""))
    }
    x<-p
  }
  if (!is.null(n)){
    if (!is.numeric(n)|n!=floor(n)){
      stop(paste("n must be an integer", "\n", ""))
    }
    x<-n
  }

  if (!is.numeric(mu)){
    stop(paste("mu must be numeric", "\n", ""))
  }

  if (any(sigma_v <= 0)){
    stop(paste("sigma_v must be positive", "\n", ""))
  }

  if (any(par_u <= 0)){
    par_u_name<-ifelse(family=="normhnorm","sigma_u","lambda")
    stop(paste(par_u_name,"must be positive", "\n", ""))
  }

  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }

  if (!family%in%c("normhnorm","normexp")){
    stop(paste("Incorrect family", "\n", ""))
  }

  if (!is.numeric(deriv)|deriv<0){
    stop(paste("deriv must be a numeric greater than 0", "\n", ""))
  }

  #Recycle parameters
  argument_list<-list(x, mu, sigma_v, par_u)
  argument_list_length<-sapply(argument_list,length)
  max_length<-max(argument_list_length)
  if(any(argument_list_length)<max_length){
    recycle<-TRUE
    if(is.null(n)){
      x<-rep(x,length.out=max_length)
    }
    mu<-rep(mu,length.out=max_length)
    sigma_v<-rep(sigma_v,length.out=max_length)
    par_u<-rep(par_u,length.out=max_length)
  } else {
    recycle<-FALSE
  }

  out<-list(x=x, mu=mu, sigma_v=sigma_v, par_u=par_u, recycle=recycle)

  #Return output
  return(out)
}
