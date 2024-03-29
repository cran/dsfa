#' elasticity
#'
#' Calculates and plots the elasticity of a smooth function.
#'
#' @param object fitted mgcv object with family \code{comper()} or \code{comper_mv()}.
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
#' #Get elasticities
#' elasticity(model, plot=TRUE)
#' }
#' @references
#' \itemize{
#' \item \insertRef{schmidt2023multivariate}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{aigner1977formulation}{dsfa}
#' \item \insertRef{meeusen1977efficiency}{dsfa}
#' }
#' @export

elasticity<-function(object, select=NULL, plot=TRUE, se=TRUE){
  #Calculates the elasticity of the inputs.
  #If no select is provided elasticity for all non-parametric terms is returned.
  jj<-attr(stats::model.matrix(object),"lpi")
  
  prod_fun_index<-jj[[1]]
  
  if(length(object$family$distr)>1){
    prod_fun_index<-c(prod_fun_index,jj[[4]])
  }
  
  if(is.null(select)){
    select<-c()
    for(i in 1:length(object$smooth)){
      if((!any(c("mrf.smooth", "random.effect")%in%attr(object$smooth[[i]],"class"))&object$smooth[[i]]$first.para%in%prod_fun_index)){
        select<-c(select,i)
      }
    }
  } 
  
  if(!se){
    out<-matrix(rep(NA,1*length(select)), ncol=1)
    colnames(out)<-c("elasticity")
  } else {
    out<-matrix(rep(NA,3*length(select)), ncol=3)
    colnames(out)<-c("elasticity","CI_lower","CI_upper")
  }
  rownames(out)<-select
  
  for(i in select){
    if(i!=select[1]){
      if(plot){
        invisible(readline(prompt="Hit <Return> to see next plot"))
      }
    }
    
    ela<-gratia::derivatives(object, type = "central", term = i)
    term<-object$smooth[[i]]$term
    
    if(plot){
      plot(ela$data, ela$derivative, type="l", xlab=term, ylab=paste0("h'(",term,")"), main=paste0("Elasticity of ",term, " with average = ", round(mean(ela$derivative),4)))
      
      if(se){
        graphics::lines(ela$data,ela$derivative-ela$se, lty=2)
        graphics::lines(ela$data,ela$derivative+ela$se, lty=2)
      }
      
    } else {
      out[which(select%in%i),]<-c(mean(ela$derivative),mean(ela$lower),mean(ela$upper))
      rownames(out)[which(select%in%i)]<-term
    }
  }

  #Return average derivatives
  if(!plot){
   return(out)
  }
}

