#' efficiency
#'
#' Calculates the expected technical (in)efficiency index.
#'
#' @param object fitted mgcv object with family \code{normhnorm} or \code{normexp}.
#' @param type default is "jondrow" for \eqn{E[u|\epsilon]}, alternatively "battese" for \eqn{E[\exp(-u)|\epsilon]}.
#' @param level for the \eqn{(1-level) \cdot 100\%} confidence interval. Must be in (0,1).
#'
#' @return Returns a matrix of the expected (in)efficiency estimates as well the lower and upper bound of the \eqn{(1-level)\cdot 100\%} confidence interval.
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
#' #Simulate responses and create dataset
#' y<-rnormhnorm(n=N, mu=mu, sigma_v=sigma_v, sigma_u=sigma_u, s=s)
#' dat<-data.frame(y, x1, x2, x3, x4, x5)
#'
#' #Write formulae for parameters
#' mu_formula<-y~x1+x2+I(x2^2)+s(x3, bs="ps")
#' sigma_v_formula<-~1+x4
#' sigma_u_formula<-~1+s(x5, bs="ps")
#'
#' #Fit model
#' model<-mgcv::gam(formula=list(mu_formula, sigma_v_formula, sigma_v_formula),
#'                  data=dat, family=normhnorm(s=s), optimizer = c("efs"))
#'
#' #Estimate efficiency
#' efficiency(model, type="jondrow")
#' efficiency(model, type="battese")
#' }
#' @references
#' \itemize{
#' \item \insertRef{schmidt2022mvdsfm}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' \item \insertRef{jondrow1982estimation}{dsfa}
#' \item \insertRef{battese1988prediction}{dsfa}
#' }
#' @export
efficiency<-function (object, level=0.05, type="jondrow"){
  #Takes fitted object and calculates the expected
  #technical efficiency of each production unit
  #based on Stochastic Frontier Analysis using Stata
  mu<-object$fitted.values[,1]
  sigma_v<-object$fitted.values[,2]
  y<-object$y

  out<-NULL

  if(object$family$family=="normhnorm"){
    sigma_u<-object$fitted.values[,3]
    sigma_c<-sqrt((sigma_u^2*sigma_v^2)/(sigma_v^2+sigma_u^2))#(1/sigma_u^2+1/sigma_v^2)^(-1)
    mu_c<-object$family$s*mgcv::residuals.gam(object)/sigma_v^2*sigma_c^2#object$family$s*mgcv::residuals.gam(object)/sigma_v^2*sigma_c^2
  }

  if(object$family$family=="normexp"){
    lambda<-object$fitted.values[,3]
    sigma_c<-sigma_v#(1/sigma_u^2+1/sigma_v^2)^(-1)
    mu_c<-object$family$s*(mgcv::residuals.gam(object)+sigma_v^2*lambda)
  }

  #If type Jondrow et al.
  if(type=="jondrow"){
    u<-mu_c+sigma_c*sn::zeta(1,mu_c/sigma_c)
  }

  #Counter numerical over- and underflow of qnorm
  lower<-1-(1-level/2)*(1-stats::pnorm(-mu_c/sigma_c))
  lower[lower>=1-1e-16]<-1-1e-16
  upper<-1-level/2*(1-stats::pnorm(-mu_c/sigma_c))
  upper[upper>=0+1e-16]<-0+1e-16

  CI_lower<-mu_c+stats::qnorm(lower)
  CI_upper<-mu_c+stats::qnorm(upper)

  # CI_lower<-mu_c+qnorm(1-(1-level/2)*(1-pnorm(-mu_c/sigma_c)))
  # CI_upper<-mu_c+qnorm(1-level/2*(1-pnorm(-mu_c/sigma_c)))
  #
  #If type Battese & Coelli
  if(type=="battese"){
    u<-exp(1/2*sigma_c^2-mu_c)*stats::pnorm(mu_c/sigma_c-sigma_c)/stats::pnorm(mu_c/sigma_c)
    CI_lower<-exp(-CI_upper)
    CI_upper<-exp(-CI_lower)
  }

  #Combine to data.frame
  out<-data.frame(u=u, CI_lower=CI_lower, CI_upper=CI_upper)

  if(is.null(out)){
    stop(paste("Wrong inputs for the efficiency", "\n", ""))
  }

  #Return output
  return(out)
}
