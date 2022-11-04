#' rsp function
#'
#' Response function and inverse response function (link function).
#'
#' @param x vector of quantiles.
#' @param link specifying the type of link function.\cr
#' `identitity` = Identity link function \cr
#' `log` = Log link function \cr
#' `glogit` = Generalized logit link function with parameters \code{a} and \code{b} \cr
#' @param a numeric min value for \code{glogit} link function.
#' @param b numeric  max value for \code{glogit} link function.
#' @param inv logical; if TRUE, evaluate link function and its derivatives.
#' @param deriv derivative of order \code{deriv}. Available are \code{0}, \code{2} and \code{4}.
#'
#' @details The link functions are defined as follows:\cr
#' `identitity`: \eqn{g(x)=x} and thus the response function is \eqn{g^{-1}(x)=x}. \cr
#' `log`: \eqn{g(x)=log(x)} and thus the response function is \eqn{g^{-1}(x)=exp(x)}. \cr
#' `glogit`: \eqn{g(x)=log(\frac{x-a}{b-a}/(1-\frac{x-a}{b-a}))} and thus the response function is \eqn{g^{-1}(x)=\frac{exp(x)}{1+exp(x)} \cdot (b-a)+a}. \cr
#'
#' @return Mostly internal function. Returns a vector of the corresponding function evaluated at x with its derivatives.
#'
#' @examples
#' rsp(x=5, link="glogit", a=0, b=1, inv=FALSE, deriv=4)
#' rsp(x=0.5, link="glogit", a=0, b=1, inv=TRUE, deriv=4)
#'
#' @export
rsp<-function(x, link="identity", a=0, b=1, inv=FALSE, deriv=0){

  link<-match.arg(link,c("identity", "log", "glogit"))

  n<-length(x)

  #Identity link fun
  if(link=="idenitity"){
    out<-x

    if(deriv>0){
      attr(out,"gradient")<-matrix(1,ncol = 1, nrow=n)
      attr(out,"hessian")<-matrix(0,ncol = 1, nrow=n)
      if(deriv>2){
        attr(out,"l3")<-matrix(0,ncol = 1, nrow=n)
        attr(out,"l4")<-matrix(0,ncol = 1, nrow=n)
      }
    }
  }

  #Generalized logit fun
  if(link=="glogit"){
    if(!inv){
      out <- exp(x)/(1 + exp(x))
      out <- ifelse(is.na(out) & !is.na(x), 1, out)
      out <- out * (b - a) + a#exp(x)/(1+exp(x))*(b-a)+a

      if(deriv>0){
        attr(out,"gradient")<-matrix(exp(x)*(b-a)/(exp(x)+1)^2, ncol = 1, nrow=n)
        attr(out,"hessian")<-matrix(exp(x)*(exp(x)-1)*(a-b)/(exp(x)+1)^3, ncol = 1, nrow=n)
        if(deriv>2){
          attr(out,"l3")<-matrix(-((6 * (-a + b) * exp(4 *x))/(1 + exp(x))^4) + (12 * (-a + b) * exp(3 * x))/(1 + exp(x))^3 - (7 * (-a + b) * exp(2 * x))/(1 + exp(x))^2 + ((-a + b) * exp(x))/(1 + exp(x)), ncol = 1, nrow=n)
          attr(out,"l4")<-matrix((24 * (-a + b) * exp(5 * x))/(1 + exp(x))^5 - (60 * (-a + b) * exp(4 * x))/(1 + exp(x))^4 + (50 * (-a + b) * exp(3 * x))/(1 + exp(x))^3 - (15 * (-a + b) * exp(2 * x))/(1 + exp(x))^2 + ((-a + b) * exp(x))/(1 + exp(x)), ncol = 1, nrow=n)
        }
      }
    } else {
      out<-(x-a)/(b-a)
      out<-log(out/(1-out)) #log((x-a)/(b-a)/(1-(x-a)/(b-a)))

      if(deriv>0){
        attr(out,"gradient")<-matrix(1/(b - x) + 1/(-a + x), ncol = 1, nrow=n)
        attr(out,"hessian")<-matrix(-(1/(a - x)^2) + 1/(b - x)^2, ncol = 1, nrow=n)
        if(deriv>2){
          attr(out,"l3")<-matrix(2 * (1/(b - x)^3 + 1/(-a + x)^3), ncol = 1, nrow=n)
          attr(out,"l4")<-matrix(-(6/(a - x)^4) + 6/(b - x)^4, ncol = 1, nrow=n)
        }
      }
    }
  }

  #Log link fun
  if(link=="log"){
    if(!inv){
      out<-exp(x)

      if(deriv>0){
        attr(out,"gradient")<-matrix(exp(x),ncol = 1, nrow=n)
        attr(out,"hessian")<-matrix(exp(x),ncol = 1, nrow=n)
        if(deriv>2){
          attr(out,"l3")<-matrix(exp(x),ncol = 1, nrow=n)
          attr(out,"l4")<-matrix(exp(x),ncol = 1, nrow=n)
        }
      }

    } else {
      out<-log(x)

      if(deriv>0){
        attr(out,"gradient")<-matrix(1/x,ncol = 1, nrow=n)
        attr(out,"hessian")<-matrix(-1/x^2,ncol = 1, nrow=n)
        if(deriv>2){
          attr(out,"l3")<-matrix(2/x^3,ncol = 1, nrow=n)
          attr(out,"l4")<-matrix(-6/x^4,ncol = 1, nrow=n)
        }
      }
    }
  }

  #return
  return(out)
}

