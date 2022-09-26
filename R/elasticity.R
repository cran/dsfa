#' elasticity
#'
#' Calculates and plots the elasticity of a smooth function.
#'
#' @param object fitted mgcv object with family \code{normhnorm()}, \code{normexp()} or \code{comperr_mv()}.
#' @param select specifying the smooth function for which the elasticity is calculated. If \code{term=NULL} the elasticities for all smooths of \eqn{\mu} are returned (excluding random and spatial effects).
#' @param plot logical; if TRUE, plots the elasticities. If FALSE, returns the average elasticity.
#' @param se logical; if TRUE, adds standard errors to the plot of elasticities.
#'
#' @return If plot is TRUE, plots the elasticities specified in select of the provided object. If plot is FALSE returns a named vector of the elasticity of the provided inputs.
#'
#' @details Calculates the marginal product for parametric terms. For smooth terms the average of the derivative is calculated.
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
#' model<-mgcv::gam(formula=list(mu_formula, sigma_v_formula, sigma_u_formula),
#'                  data=dat, family=normhnorm(s=s), optimizer = c("efs"))
#'
#' #Get elasticities
#' elasticity(model)
#' }
#' @references
#' \itemize{
#' \item \insertRef{schmidt2022mvdsfm}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{aigner1977formulation}{dsfa}
#' \item \insertRef{meeusen1977efficiency}{dsfa}
#' }
#' @export

elasticity<-function(object, select=NULL, plot=TRUE, se=TRUE){
  #Calculates the elasticity of the inputs.
  #If no select is provided elasticity for all non-parametric terms is returned.
  if(object$family$family=="joint"){
    par_index<-c(attr(object$nsdf,"pstart")[1]:attr(object$nsdf,"pstart")[2], attr(object$nsdf,"pstart")[4]:attr(object$nsdf,"pstart")[5])
  } else {
    par_index<-attr(object$nsdf,"pstart")[1]:attr(object$nsdf,"pstart")[2]
  }
  if(is.null(select)){
     out<-c()
     for(i in 1:gratia::n_smooths(object)){
       if((any(c("mrf.smooth", "random.effect")%in%attr(object$smooth[[i]],"class"))|object$smooth[[i]]$first.para%in%par_index)){
         if(i>1){
           invisible(readline(prompt="Hit <Return> to see next plot"))
         }
         ela<-gratia::derivatives(object, type = "central", term = i)
         term<-object$smooth[[i]]$term
         if(plot){
           plot(ela$data,ela$derivative, type="l", xlab=term, ylab=paste0("h'(",term,")"), main=paste0("Elasticity of ",term, " with average=",round(mean(ela$derivative),4)))
           if(se){
             graphics::lines(ela$data,ela$derivative-ela$se, lty=2)
             graphics::lines(ela$data,ela$derivative+ela$se, lty=2)
           }
           # invisible(readline(prompt="Hit <Return> to see next plot"))
         } else {
           out<-c(out, mean(ela$derivative))
           names(out)[i]<-term
         }
       }
     }
  } else {
    if((any(c("mrf.smooth", "random.effect")%in%attr(object$smooth[[i]],"class"))|object$smooth[[i]]$first.para%in%par_index)){
      i<-select
      ela<-gratia::derivatives(object, type = "central", term = i)
      term<-object$smooth[[i]]$term
      if(plot){
        plot(ela$data,ela$derivative, type="l", xlab=term, ylab=paste0("h'(",term,")"), main=paste0("Elasticity of ",term, " with average=",round(mean(ela$derivative),4)))
        if(se){
          graphics::lines(ela$data,ela$derivative-ela$se, lty=2)
          graphics::lines(ela$data,ela$derivative+ela$se, lty=2)
        }
      } else {
        out<-mean(ela$derivative)
        names(out)<-term
      }
    } else {
      stop(paste("Incorrect input for the argument select", "\n", ""))
    }
  }

  #Return average derivatives
  if(!plot){
   return(out)
  }
  # return(out)
}

