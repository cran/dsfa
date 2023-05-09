#' comper
#'
#' The comper implements the composed-error distribution in which the \eqn{\mu}, \eqn{\sigma_V} and \eqn{\sigma_U} can depend on additive predictors.
#' Useable with \code{mgcv::gam}, the additive predictors are specified via a list of formulae.
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the 'mgcv' and 'dsfa' package.
#'
#' @details Used with \code{\link[mgcv:gam]{gam()}} to fit distributional stochastic frontier model. The function is called with a list containing three formulae:
#' \enumerate{
#'   \item The first formula specifies the response on the left hand side and the structure of the additive predictor for \eqn{\mu} parameter on the right hand side. Link function is "identity".
#'   \item The second formula is one sided, specifying the additive predictor for the  \eqn{\sigma_V} on the right hand side. Link function is "log".
#'   \item The third formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_U} on the right hand side. Link function is "log".
#' }
#' The fitted values and linear predictors for this family will be three column matrices. The first column is the \eqn{\mu}, the second column is the \eqn{\sigma_V}, the third column is \eqn{\sigma_U}.
#' For more details of the distribution see \code{dcomper()}.
#'
#' @inheritParams dcomper
#' @param link three item list, specifying the link for the \eqn{\mu}, \eqn{\sigma_V} and \eqn{\sigma_U} parameters. See details.
#' 
#' @examples
#' ### First example with simulated data
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
#' y<-rcomper(n=N, mu=mu, sigma_v=sigma_v, sigma_u=sigma_u, s=s, distr="normhnorm")
#' dat<-data.frame(y, x1, x2, x3, x4, x5)
#' 
#' #Write formulae for parameters
#' mu_formula<-y~x1+x2+I(x2^2)+s(x3, bs="ps")
#' sigma_v_formula<-~1+x4
#' sigma_u_formula<-~1+s(x5, bs="ps")
#' 
#' #Fit model
#' model<-mgcv::gam(formula=list(mu_formula, sigma_v_formula, sigma_u_formula),
#'                  data=dat, family=comper(s=s, distr="normhnorm"), optimizer = c("efs"))
#' 
#' #Model summary
#' summary(model)
#' 
#' #Smooth effects
#' #Effect of x3 on the predictor of the production function
#' plot(model, select=1) #Estimated function
#' lines(x3[order(x3)], 6*log(x3[order(x3)]+2)^(1/4)-
#'         mean(6*log(x3[order(x3)]+2)^(1/4)), col=2) #True effect
#' 
#' #Effect of x5 on the predictor of the inefficiency
#' plot(model, select=2) #Estimated function
#' lines(x5[order(x5)], -1+sin(2*pi*x5)[order(x5)]-
#'         mean(-1+sin(2*pi*x5)),col=2) #True effect
#' 
#' ### Second example with real data
#' \donttest{
#' data("RiceFarms", package = "plm") #load data
#' RiceFarms[,c("goutput","size","seed", "totlabor", "urea")]<-
#'   log(RiceFarms[,c("goutput","size","seed", "totlabor", "urea")]) #log outputs and inputs
#' RiceFarms$id<-factor(RiceFarms$id) #id as factor
#' 
#' #Set to production function
#' s=-1 
#' 
#' #Write formulae for parameters
#' mu_formula<-goutput ~  s(size, bs="ps") + s(seed, bs="ps") + #non-linear effects
#'   s(totlabor, bs="ps") + s(urea, bs="ps") + #non-linear effects
#'   varieties + #factor
#'   s(id, bs="re") #random effect
#' sigma_v_formula<-~1 
#' sigma_u_formula<-~bimas
#' 
#' #Fit model with normhnorm dstribution
#' model_normhnorm<-mgcv::gam(formula=list(mu_formula, sigma_v_formula, sigma_u_formula),
#' data=RiceFarms, family=comper(s=-1, distr="normhnorm"), optimizer = "efs")
#' 
#' #Summary of model
#' summary(model_normhnorm)
#' 
#' #Plot smooths
#' plot(model_normhnorm)
#' }
#' 
#' ### Third example with real data of cost function
#' \donttest{
#' data("electricity", package = "sfaR") #load data
#' 
#' #Log inputs and outputs as in Greene 1990 eq. 46
#' electricity$lcof<-log(electricity$cost/electricity$fprice)
#' electricity$lo<-log(electricity$output)
#' electricity$llf<-log(electricity$lprice/electricity$fprice)
#' electricity$lcf<-log(electricity$cprice/electricity$fprice)
#' 
#' #Set to cost function
#' s=1
#' 
#' #Write formulae for parameters
#' mu_formula<-lcof ~ s(lo, bs="ps") + s(llf, bs="ps") + s(lcf, bs="ps") #non-linear effects
#' sigma_v_formula<-~1
#' sigma_u_formula<-~s(lo, bs="ps") + s(lshare, bs="ps") + s(cshare, bs="ps")
#' 
#' #Fit model with normhnorm dstribution
#' model_normhnorm<-mgcv::gam(formula=list(mu_formula, sigma_v_formula, sigma_u_formula),
#'                            data=electricity, family=comper(s=s, distr="normhnorm"),
#'                            optimizer = "efs")
#' 
#' #Summary of model
#' summary(model_normhnorm)
#' 
#' #Plot smooths
#' plot(model_normhnorm)
#' }
#' 
#' @references
#' \itemize{
#' \item \insertRef{schmidt2022mvdsfm}{dsfa}
#' \item \insertRef{wood2017generalized}{dsfa}
#' \item \insertRef{aigner1977formulation}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' \item \insertRef{schmidt2020analytic}{dsfa}
#' }
#' @export
#comper distribution object for mgcv
comper<- function(link = list("identity", "log", "log"), s = -1, distr="normhnorm"){
  #Object for mgcv::gam such that the composed-error distribution can be estimated.
  
  #Number of parameters
  npar <- 3
  
  #Link functions
  if (length(link) != npar) stop("comper requires 3 links specified as character strings")
  okLinks <- list("identity", "log", "log")
  stats <- list()
  param.names <- c("mu", "sigma_v", "sigma_u")
  for (i in 1:npar) { # Links for mu, sigma_v, sigma_u
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- stats::make.link(link[[i]]) else
      stop(link[[i]]," link not available for ", param.names[i]," parameter of comper")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- mgcv::fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  
  residuals <- function(object, type = c("deviance", "response")) {
    #Calculate residuals
    type <- match.arg(type)
    
    mom<-par2mom(mu=object$fitted[,1], sigma_v = object$fitted[,2], sigma_u=object$fitted[,3], s=object$family$s, distr=object$family$distr)
   
    y_hat<-mom[,1]
    
    rsd <- object$y-y_hat
    
    if (type=="response"){
      out<-rsd
    }  else {
      out<-rsd/mom[,2]
    }
    
    return(out)
  }
  
  
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    #Loglike function with derivatives
    #function defining the comper model loglike
    #deriv: 0 - eval
    #       1 - grad and Hess
    #       2 - diagonal of first deriv of Hess
    #       3 - first deriv of Hess
    #       4 - everything
    
    # If offset is not null or a vector of zeros, give an error
    if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")
    
    #Extract linear predictor index
    jj <- attr(X, "lpi") 
    
    #Number of parameters and observations
    npar <- 3
    n <- length(y)
    
    #Get additive predictors
    eta<-matrix(0,n,npar)
    for(i in 1:npar){
      eta[,i]<-drop( X[ , jj[[i]], drop=FALSE] %*% coef[jj[[i]]] )
    }
    
    #Additive predictors 2 parameters
    theta<-matrix(0,n,npar)
    for(i in 1:npar){
      theta[,i]<-family$linfo[[i]]$linkinv( eta[,i] )
    }
    
    deriv_order<-ifelse(deriv==1,2,4)
    
    #Evaluate ll
    l0<-dcomper_cpp(x=y, m=theta[,1], v=theta[,2], u=theta[,3], s=family$s, distr=family$distr, deriv_order=deriv_order, tri=family$tri_mat, logp = TRUE)
    
    #Assign sum of individual loglikehood contributions to l
    l<-sum(l0)
    
    if (deriv>0) {
      #First derivatives
      ig1<-matrix(0,n,npar)
      for(i in 1:npar){
        ig1[,i]<-family$linfo[[i]]$mu.eta(eta[,i])
      }
      
      #Second derivatives
      g2<-matrix(0,n,npar)
      for(i in 1:npar){
        g2[,i]<-family$linfo[[i]]$d2link(theta[,i])
      }
      
    }
    
    #Assign first and second derivative of ll to l1 and l2
    l1<-attr(l0,"d1")
    l2<-attr(l0,"d2")
    
    #Set default values 
    l3 <- l4 <- g3 <- g4 <- 0
    
    if (deriv>1) {
      #Third derivatives
      g3<-matrix(0,n,npar)
      for(i in 1:npar){
        g3[,i]<-family$linfo[[i]]$d3link(theta[,i])
      }
      
      #Assign third derivative of ll to l3
      l3<-attr(l0,"d3")
    }
    
    if (deriv>3) {
      #Fourth derivatives
      g4<-matrix(0,n,npar)
      for(i in 1:npar){
        g4[,i]<-family$linfo[[i]]$d4link(theta[,i])
      }  
      
      #Assign fourth derivative of ll to l4
      l4<-attr(l0,"d4")
    }
    
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      
      #Transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- mgcv::gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)
      
      #Calculate the gradient and Hessian...
      out <- mgcv::gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                             d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D)
      
    } else {
      out <- list()
    }
    out$l <- l
    
    return(out)
  } 
  
  
  initialize <- expression({
    #Function to calculate starting values of the coefficients
    #Idea is to get starting values utilizing the method of moments,
    
    n <- rep(1, nobs)
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))){
      TRUE
    }  else {
      FALSE
    }
    
    if (is.null(start)) {
      #Extract linear predictor index
      jj <- attr(x,"lpi")
      
      #Initialize start vector
      start <- rep(0,ncol(x))
      
      #Linear model for starting values
      m0<-lm(y~x)
      
      #Calculate skewness and bound the estimates
      skew<-family$s*abs(mean(m0$residuals^3))/(sum(m0$residuals^2)/(length(y)-1))^(3/2)
      
      if(family$distr=="normhnorm"){
        max.skew<-0.5 * (4 - pi) * (2/(pi - 2))^1.5 - (.Machine$double.eps)^(1/4)
      }
      
      if(family$distr=="normexp"){
        max.skew<-2- (.Machine$double.eps)
      }
      
      if (any(abs(skew) > max.skew)){
        skew[abs(skew)>max.skew]<-family$s * 0.9 * max.skew
      }
      
      #Transform moments to parameters
      par<-mom2par(mean = m0$fitted.values,
                   sd = mean(abs(m0$residuals)),
                   skew = skew,
                   s=family$s,
                   distr=family$distr)
      
      # 1) Ridge regression for mean
      #Extract covariates and square root of total penalty associated with additive predictor for mu
      x1 <- x[ , jj[[1]], drop=FALSE]
      e1 <- E[ , jj[[1]], drop=FALSE]
      
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startMu <- qr.coef(qr(x1), c(par[,1], rep(0,nrow(E))))
        startMu[ !is.finite(startMu) ] <- 0
      } else { 
        startMu <- pen.reg(x1, e1, par[,1])
      }
      start[jj[[1]]] <- startMu
      
      # 2) Ridge regression using log absolute residuals
      #Extract covariates and square root of total penalty associated with additive predictor for sigma_v
      x2 <- x[,jj[[2]],drop=FALSE]
      e2 <- E[,jj[[2]],drop=FALSE]
      
      if (use.unscaled) {
        x2 <- rbind(x2,e2)
        startSigma_v <- qr.coef(qr(x2),c(log(par[,2]),rep(0,nrow(E))))
        startSigma_v[!is.finite(startSigma_v)] <- 0
      } else { 
        startSigma_v <- pen.reg(x2,e2,log(par[,2])) 
      }
      start[jj[[2]]] <- startSigma_v
      
      # 3) Skewness
      #Extract covariates and square root of total penalty associated with additive predictor for sigma_u
      x3 <-  x[, jj[[3]],drop=FALSE]
      e3 <-  E[, jj[[3]],drop=FALSE]
      
      if (use.unscaled) {
        x3 <- rbind(x3,e3)
        startSigma_u <- qr.coef(qr(x3), c(log(par[,3]),rep(0,nrow(E))))
        startSigma_u[!is.finite(startSigma_u)] <- 0
      } else {
        startSigma_u <- pen.reg(x3,e3,log(par[,3])) 
      }
      start[jj[[3]]] <- startSigma_u
    }
  }) 
  
  rd <- function(mu, wt, scale) {
    #Random number generation for comper
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    
    return(rcomper(n=nrow(mu), mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=attr(mu,"s"), distr=attr(mu,"distr")))
  }
  
  qf <- function(p, mu, wt, scale) {
    #Quantile function of comper
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    
    return(qcomper(p=p, mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=attr(mu,"s"), distr=attr(mu,"distr")))
  }
  
  cdf <- function(q, mu, wt, scale) {
    #Cumulative distribution function of comper
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    
    return(pcomper(q=q, mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=attr(mu,"s"), distr=attr(mu,"distr")))
  }
  
  predict <- function(family, se=FALSE,eta=NULL,y=NULL,X=NULL,beta=NULL,off=NULL,Vb=NULL) {
    #Prediction function
    # optional function to give predicted values - idea is that
    # predict.gam(...,type="response") will use this, and that
    # either eta will be provided, or {X, beta, off, Vb}. family$data
    # contains any family specific extra information.
    # if se = FALSE returns one item list containing matrix otherwise
    # list of two matrices "fit" and "se.fit"...
    
    if (is.null(eta)) {
      if (is.null(off)){
        off <- list(0,0,0)
      } 
      off[[4]] <- 0
      
      for (i in 1:family$npar) {
        if (is.null(off[[i]])){
          off[[i]] <- 0
        } 
      } 
      
      #Extract linear predictor index
      jj <- attr(X,"lpi")
      
      #Calculate additive predictors for mu, sigma_v, sigma_u
      eta<-matrix(0, nrow(X), family$npar)
      for(i in 1:family$npar){
        eta[,i]<-drop(X[,jj[[i]],drop=FALSE]%*%beta[jj[[i]]] + off[[i]])
      }
      
      if (se) {
        stop("se still not available for this family")
      }
    } else {
      se <- FALSE
    }
    
    #Additive predictors 2 parameters
    theta<-matrix(0,nrow(X),family$npar)
    for(i in 1:family$npar){
      theta[,i]<-family$linfo[[i]]$linkinv( eta[,i] )
    }
    
    #Calculate moments from parameters
    mom<-par2mom(mu = theta[,1], sigma_v = theta[,2], sigma_u = theta[,3], s=family$s, distr=family$distr)
    
    #Assign mean
    fv <- list(mom[,1])
    if (!se) return(fv) else {
      stop("se not still available for this family")
      #Assign mean and standard deviation
      fv <- list(fit=mom[,1], se.fit=mom[,2])
      return(fv)
    }
  }
  
  postproc <- expression({
    attr(object$fitted.values,"s")<-object$family$s
    attr(object$fitted.values,"distr")<-object$family$distr
    object$fitted.values
  })
  
  structure(list(family="comper", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), # symmetric indices for accessing derivative arrays
                 tri_mat = trind_generator(npar), # symmetric indices for accessing derivative matrices
                 initialize=initialize, #initial parameters
                 s = s, #production or cost function
                 distr = distr, #specifiying distribution
                 residuals=residuals, #residual function
                 rd=rd, #random number generation
                 qf=qf, #quantile function
                 cdf=cdf, #cdf function
                 predict=predict, #prediction function for mgcv
                 postproc=postproc, #Assigning attributes such that other functions work
                 linfo = stats, # link information list
                 d2link=1, d3link=1, d4link=1, # signals to fix.family.link that all done
                 ls=1, # signals that ls not needed here
                 available.derivs = 2), # can use full Newton here
  class = c("general.family","extended.family","family"))
}

