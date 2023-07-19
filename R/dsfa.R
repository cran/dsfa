#' Sets up and fitting of distributional stochastic frontier models 
#'
#' @return Returns a  \code{\link[mgcv:gam]{gam}} object.
#' 
#' @details This function is a wrapper for \code{\link[mgcv:gam]{gam}}. 
#'   
#' @param formula A list of formulas specifying the additive predictors. See \code{\link[mgcv:formula.gam]{formula.gam}} and \code{\link[mgcv:gam.models]{gam.models}} for more details.
#' @param family The family object specifies the (multivariate) composed-error distribution and link of the model. See \code{\link{comper}} and \code{\link{comper_mv}} for more details.
#' @param data A data frame or list containing the model response variable and covariates required by the formula. By default the variables are taken from environment(formula): typically the environment from which \code{dsfa} is called.
#' @param optimizer An array specifying the numerical optimization method to use to optimize the smoothing parameter estimation criterion (given by method). 
#' "outer" for the more stable direct approach. "outer" can use several alternative optimizers, specified in the second element of optimizer: "newton" (default), "bfgs", "optim", "nlm" and "nlm.fd" (the latter is based entirely on finite differenced derivatives and is very slow).
#' "efs" for the extended Fellner Schall method of Wood and Fasiolo (2017).
#' @param ... other parameters of \code{\link[mgcv:gam]{gam}}
#' 
#' @export
dsfa<-function(formula, family=comper(link = list("identity", "logshift", "logshift"), s = -1, distr = "normhnorm"), data=list(), optimizer="efs",...){
  #Check if the model determined by family can be fit by the optimizer
  if((family$family=="comper_mv") & (optimizer!="efs")){
    stop(paste("comper_mv can only be fit via optimizer efs", "\n", ""))
  }
  
  #Wrapper for mgcv::gam
  out<-mgcv::gam(formula=formula, family=family, data=data, optimizer = optimizer,...)
  
  #Return output
  return(out)
}