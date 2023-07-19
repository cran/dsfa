#' efficiency
#'
#' Calculates the expected technical (in)efficiency index.
#'
#'
#' @return Returns a matrix of the expected (in)efficiency estimates as well the lower and upper bound of the \eqn{(1-\alpha)\cdot 100\%} confidence interval.
#'
#' @param type default is "jondrow" for \eqn{E[U|\mathcal{E}]}, alternatively "battese" for \eqn{E[\exp(-U)|\mathcal{E}]}.
#' @param alpha for the \eqn{(1-\alpha) \cdot 100\%} confidence interval. Must be in (0,1).
#' @inheritParams elasticity
#' 
#' @examples
#' \donttest{
#' #Set seed, sample size and type of function
#' set.seed(1337)
#' N=500 #Sample size
#' s=-1 #Set to production function
#'
#' #Generate covariates
#' x1<-runif(N,-1,1); x2<-runif(N,-1,1); x3<-runif(N,-1,1)
#' x4<-runif(N,-1,1); x5<-runif(N,-1,1)
#'
#' #Set parameters of the distribution
#' mu=2+0.75*x1+0.4*x2+0.6*x2^2+6*log(x3+2)^(1/4) #production function parameter
#' sigma_v=exp(-1.5+0.75*x4) #noise parameter
#' sigma_u=exp(-1+sin(2*pi*x5)) #inefficiency parameter
#'
#' y<-rcomper(n=N, mu=mu, sigma_v=sigma_v, sigma_u=sigma_u, s=s, distr="normhnorm")
#' dat<-data.frame(y, x1, x2, x3, x4, x5)
#'
#' #Write formulae for parameters
#' mu_formula<-y~x1+x2+I(x2^2)+s(x3, bs="ps")
#' sigma_v_formula<-~1+x4
#' sigma_u_formula<-~1+s(x5, bs="ps")
#'
#' #Fit model
#' model<-dsfa(formula=list(mu_formula, sigma_v_formula, sigma_u_formula),
#'                  data=dat, family=comper(s=s, distr="normhnorm"), optimizer = c("efs"))
#'                                    
#' #Estimate efficiency
#' efficiency(model, type="jondrow")
#' efficiency(model, type="battese")
#' }
#' @references
#' \itemize{
#' \item \insertRef{schmidt2023multivariate}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' \item \insertRef{jondrow1982estimation}{dsfa}
#' \item \insertRef{battese1988prediction}{dsfa}
#' }
#' @export
efficiency<-function (object, alpha=0.05, type="jondrow"){
  #Takes fitted object and calculates the expected
  #technical efficiency of each production unit
  #based on Stochastic Frontier Analysis using Stata
  mu<-object$fitted.values[,1, drop=F]
  sigma_v<-object$fitted.values[,2, drop=F]
  sigma_u<-object$fitted.values[,3, drop=F]
  y<-as.matrix(object$y)

  if(length(object$family$distr)>1){
    mu<-cbind(mu,object$fitted.values[,4])
    sigma_v<-cbind(sigma_v,object$fitted.values[,5])
    sigma_u<-cbind(sigma_u,object$fitted.values[,6])
  }
  
  out<-c()#matrix(0,nrow(y),1)
  
  for(i in 1:ncol(y)){
    if(object$family$distr[i]=="normhnorm"){
      sigma_c<-sqrt((sigma_u[,i]^2*sigma_v[,i]^2)/(sigma_v[,i]^2+sigma_u[,i]^2))#(1/sigma_u^2+1/sigma_v^2)^(-1)
      mu_c<-object$family$s[i]*(y-mu)[,i]/sigma_v[,i]^2*sigma_c^2#object$family$s*mgcv::residuals.gam(object)/sigma_v^2*sigma_c^2
    }
    
    if(object$family$distr[i]=="normexp"){
      sigma_c<-sigma_v[,i]
      mu_c<-object$family$s[i]*((y-mu)[,i]+sigma_v[,i]^2*sigma_u[,i])
    }
    
    #If type Jondrow et al.
    if(type=="jondrow"){
      u<-mu_c+sigma_c*stats::dnorm(mu_c/sigma_c)/stats::pnorm(mu_c/sigma_c)
    }
    
    #Counter numerical over- and underflow of qnorm
    lower<-1-(1-alpha/2)*(1-stats::pnorm(-mu_c/sigma_c))
    lower[lower>=1-1e-16]<-1-1e-16
    upper<-1-alpha/2*(1-stats::pnorm(-mu_c/sigma_c))
    upper[upper>=0+1e-16]<-0+1e-16
    
    CI_lower<-mu_c+stats::qnorm(lower)*sigma_c
    CI_upper<-mu_c+stats::qnorm(upper)*sigma_c
    
    #If type Battese & Coelli
    if(type=="battese"){
      u<-exp(1/2*sigma_c^2-mu_c)*stats::pnorm(mu_c/sigma_c-sigma_c)/stats::pnorm(mu_c/sigma_c)
      CI_lower<-exp(-CI_upper)
      CI_upper<-exp(-CI_lower)
    }
    
    if(i==1){
      out<-cbind(u, CI_lower, CI_upper)
      colnames(out)<-c("Expected (In)efficiency", "CI_lower", "CI_upper")
    } else {
      out<-cbind(out, cbind(u, CI_lower, CI_upper))
      colnames(out)<-paste0(c("Expected (In)efficiency_", "CI_lower_", "CI_upper_"),rep(1:ncol(y), each=3))
    }
  }
  
  # if(length(object$family$distr)==1){
    
  # } else {
  #   colnames(out)<-paste0(c("Expected (In)efficiency", "CI_lower_", "CI_upper_"),rep(1:ncol(y), each=3))
  # }
  
  #Return output
  return(out)
}
