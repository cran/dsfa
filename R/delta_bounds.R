#'  Bounds of Copula Parameter delta
#' 
#' Provides the minimum and maximum of the parameter space for \eqn{\delta}
#'  
#' @return Returns numeric vector of length two with first argument being the min and the second argument being the max of the parameter space.
#' 
#' @details Although the parameter space is larger in theory for some copulas, numeric under- and overflow limits the parameter space. The parameter space of \eqn{\delta} is specified for each copula below: 
#' \enumerate{
#' \item  `independent`, min=0 and max=1
#' \item  `normal`, min=-1 and max=1
#' \item  `clayton`, min=1e-16 and max=28
#' \item  `gumbel`, min=1 and max=17
#' \item  `frank`, min=-35 and max=35
#' \item  `joe`, min=1e-16 and max=30
#' }
#' 
#' @inheritParams dcop
#' 
#' 
#' @examples 
#' delta_bounds("normal")
#' 
#' @family copula
#' 
#' @export
delta_bounds<-function(distr_cop){
  #Delta parameter space boundaries
  if(distr_cop=="independent"){
    min<-0
    max<-1
  }
  
  if(distr_cop=="normal"){
    min<--1
    max<-1
  }
  
  if(distr_cop=="clayton"){
    min<-0+1e-16
    max<-28
  }
  
  if(distr_cop=="gumbel"){
    min<-1
    max<-17
  }
  
  if(distr_cop=="frank"){
    min<--35
    max<-35
    
    stop(paste("The pdf of the frank copula is not functional yet.", "\n", ""))
    
  }
  
  if(distr_cop=="joe"){
    min<-1+1e-6
    max<-30
  }
  
  out<-c(min,max)
  
  return(out)
}