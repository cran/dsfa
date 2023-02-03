#' Moments to Parameters 
#'
#' Calculates the parameters of composed-error distribution based on the provided moments. 
#'
#' @return Returns a matrix where the first column corresponds to \eqn{\mu}, the second to \eqn{\sigma_V} and the third to \eqn{\sigma_U}.
#' 
#' @details See [dcomper()] for details of the distribution. For the inverse transformation see [par2mom()].
#' 
#' @inheritParams dcomper
#' @param mean numeric vector of means.
#' @param sd numeric vector of standard deviations. Must be positive.
#' @param skew numeric vector of skewness. \code{s*skew} must be positive.
#' 
#' @examples 
#' mom2par(mean=0, sd=1, skew=-0.5, s=-1, distr="normhnorm")
#' mom2par(mean=0, sd=1, skew=-1, s=-1, distr="normexp")
#' 
#' @export
mom2par<-function(mean=0, sd=1, skew=0, s=-1, distr="normhnorm"){
  if (any(sd <= 0)|any(s*skew <= 0)){
    stop(paste("sd and s*skew must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  distr<-match.arg(distr,c("normhnorm","normexp"))
  
  X<-tryCatch(cbind(0,mean, sd, skew), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  if(distr=="normhnorm"){
    # max.skew<-0.5 * (4 - pi) * (2/(pi - 2))^1.5 - (.Machine$double.eps)^(1/4)
    # 
    # if (any(abs(skew) > max.skew)){
    #   skew[abs(skew)>max.skew]<-s * 0.9 * max.skew
    # }
    
    R<-(2*abs(skew)/(4-pi))^(1/3)
    alpha<-s*R/sqrt(2/pi-(1-2/pi)*R^2)
    delta<-alpha/sqrt(1+alpha^2)
    b<-sqrt(2/pi)
    omega<-sd/sqrt(1-b^2*delta^2)
    xi<-mean-b*omega*delta
    mu<-xi
    sigma_u<-omega*s*alpha/sqrt((s*alpha)^2+1)
    sigma_v<-sigma_u/(s*alpha)
  }
  
  if(distr=="normexp"){
    # max.skew<-2
    # 
    # if (any(abs(skew) > max.skew)){
    #   skew[abs(skew)>max.skew]<-s * 0.9 * max.skew
    # }
    
    mu<-mean-s*sd*(s*skew/2)^(1/3)
    sigma_v<-sqrt(sd^2*(1-(s*skew/2)^(2/3)))
    sigma_u<-1/(sd*(s*skew/2)^(1/3))
  }
  
  out<-cbind(mu, sigma_v, sigma_u)
  
  #Return output
  return(out)
}