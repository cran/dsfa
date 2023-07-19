#'  Bounds of Copula Parameter delta
#' 
#' Provides the minimum and maximum of the parameter space for \eqn{\delta}
#'  
#' @return Returns numeric vector of length two with first argument being the minimum and the second argument being the maximum of the parameter space.
#' 
#' @details Although the parameter space is larger in theory for some copulas, numeric under- and overflow limits the parameter space. The parameter space of \eqn{\delta} is specified for each copula below: 
#' \itemize{
#' \item  `independent`, min=0 and max=1
#' \item  `normal`, min=-1 and max=1
#' \item  `clayton`, min=1e-16 and max=28
#' \item  `gumbel`, min=1 and max=17
#' \item  `frank`, min=-35 and max=35
#' \item  `joe`, min=1e-16 and max=30
#' \item  `amh`, min=-1 and max=1
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
delta_bounds<-function(distr){
  #Delta parameter space boundaries
  if(distr=="independent"){
    min<-0
    max<-1
  }
  
  if(distr=="normal"){
    min<--1
    max<-1
  }
  
  if(distr=="clayton"){
    min<-0
    max<-28
  }
  
  if(distr=="gumbel"){
    min<-1
    max<-17
  }
  
  if(distr=="frank"){
    min<--35
    max<-35
  }
  
  if(distr=="joe"){
    min<-1
    max<-30
  }
  
  if(distr=="amh"){
    min<--1
    max<-1
  }
  
  out<-c(min,max)
  
  return(out)
}
