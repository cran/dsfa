#' Parameter to Moments 
#'
#' Calculates the moments of composed-error distribution based on the provided parameters. 
#'
#' @return Returns a matrix where the first column corresponds to the mean, the second to the standard deviation and the third to the skewness.
#' 
#' @details See [dcomper()] for details of the distribution. For the inverse transformation see [mom2par()].
#'  
#' @inheritParams dcomper
#' 
#' @examples 
#' par2mom(mu=0, sigma_v=1, sigma_u=1, s=-1, distr="normhnorm")
#' par2mom(mu=0, sigma_v=1, sigma_u=1, s=-1, distr="normexp")
#' 
#' @export
par2mom<-function(mu=0, sigma_v=1, sigma_u=1, s=-1, distr="normhnorm"){
  if (any(sigma_v <= 0)|any(sigma_u <= 0)){
    stop(paste("sigma_v and sigma_u must be positive", "\n", ""))
  }
  
  if (any(!s%in%c(-1,1))){
    stop(paste("s must be {-1, 1}", "\n", ""))
  }
  
  distr<-match.arg(distr,c("normhnorm","normexp"))
  
  X<-tryCatch(cbind(0,mu, sigma_v, sigma_u), warning=function(w) {
    stop("Input vectors have incompatible lengths")})
  
  
  if(distr=="normhnorm"){
    omega<-sqrt(X[,3]^2+X[,4]^2)
    alpha<-s*X[,4]/X[,3]
    xi<-X[,2]
    b<-sqrt(2/pi)
    delta<-alpha/sqrt(1+alpha^2)
    mu_z<-b*delta
    sigma_z<-sqrt(1-mu_z^2)
    mean<-xi+omega*mu_z
    sd<-sqrt(omega^2*(1-mu_z^2))
    skew<-1/2*(4-pi)*mu_z^3/(1-mu_z^2)^(3/2)
  }
  
  if(distr=="normexp"){
    mean<-mu+s*1/sigma_u
    sd<-sqrt(sigma_v^2+1/sigma_u^2)
    skew<-2/(sigma_v^3*s^3*sigma_u^3)*(1+1/(sigma_v^2*sigma_u^2))^(-3/2)
  }
  
  out<-cbind(mean, sd, skew)
  
  #Return output
  return(out)
}