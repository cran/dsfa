#' comper
#'
#' The comper implements the multivariate composed-error distribution in which the \eqn{\mu_1}, \eqn{\sigma_{V1}}, \eqn{\sigma_{U2}}, \eqn{\mu_2}, \eqn{\sigma_{V2}},  \eqn{\sigma_{U2}} and \eqn{\delta} can depend on additive predictors.
#' Useable with \code{mgcv::gam}, the additive predictors are specified via a list of formulae.
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the \pkg{mgcv} and \pkg{dsfa} package.
#'
#' @details Used with \code{\link[mgcv:gam]{gam()}} to fit distributional stochastic frontier model. The function is called with a list containing three formulae:
#' \enumerate{
#'   \item The first formula specifies the response of marginal one on the left hand side and the structure of the additive predictor for \eqn{\mu_1} parameter on the right hand side. Link function is "identity".
#'   \item The second formula is one sided, specifying the additive predictor for the  \eqn{\sigma_{V1}} on the right hand side. Link function is "logshift", e.g. \eqn{\log \{ \sigma_{V1} \}  + b }.
#'   \item The third formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_{U1}} on the right hand side. Link function is "logshift", e.g. \eqn{\log \{ \sigma_{U1} \}  + b }.
#'   \item The fourth formula specifies the response of marginal two on the left hand side and the structure of the additive predictor for \eqn{\mu_2} parameter on the right hand side. Link function is "identity".
#'   \item The fifth formula is one sided, specifying the additive predictor for the  \eqn{\sigma_{V2}} on the right hand side. Link function is "logshift", e.g. \eqn{\log \{ \sigma_{V2} \}  + b }.
#'   \item The sixth formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_{U2}} on the right hand side. Link function is "logshift", e.g. \eqn{\log \{ \sigma_{U2} \}  + b }.
#'   \item The seventh formula  is one sided, specifying the additive predictor for the  \eqn{\delta} on the right hand side. Link function is "glogit".
#' }
#' The fitted values and linear predictors for this family will be seven column matrices.
#' For more details of the distribution see \code{dcomper()}.
#'
#' @inheritParams dcomper_mv
#' @param link seven item list, specifying the links for \eqn{\mu_1}, \eqn{\sigma_{V1}}, \eqn{\sigma_{U2}}, \eqn{\mu_2}, \eqn{\sigma_{V2}},  \eqn{\sigma_{U2}} and \eqn{\delta}. See details.
#' @param b positive parameter of the logshift link function.
#' @examples
#' \donttest{
#' #Set seed, sample size and type of function
#' set.seed(1337)
#' N=1000 #Sample size
#' s<-c(-1,-1) #Set to production function for margin 1 and set to cost function for margin 2
#' 
#' distr_cop="normal"
#' distr_marg1="normhnorm"
#' distr_marg2="normhnorm"
#' 
#' #Generate covariates
#' x1<-runif(N,-1,1); x2<-runif(N,-1,1); x3<-runif(N,-1,1)
#' x4<-runif(N,-1,1); x5<-runif(N,-1,1); x6<-runif(N,-1,1)
#' x7<-runif(N,-1,1)
#' 
#' mu1=6+2*x1+(-2/3)*x1^2 #production function parameter 1
#' sigma_v1=exp(-1.5+sin(2*pi*x2)) #noise parameter 1
#' sigma_u1=exp(-1) #inefficiency parameter 1
#' mu2=5*x4^2+4*log(x4+2)^(1/4) #cost function parameter 2
#' sigma_v2=exp(-1.5) #noise parameter 2
#' sigma_u2=exp(-1+sin(2*pi*x6)) #inefficiency parameter 2
#' delta=transform(x=matrix(1+2.5*cos(4*x7)),
#'       type="glogitinv",
#'       par=delta_bounds(distr_cop), deriv_order = 0)
#' 
#' #Simulate responses and create dataset
#' Y<-rcomper_mv(n=N, mu=cbind(mu1,mu2),
#'               sigma_v=cbind(sigma_v1, sigma_v2),
#'               sigma_u = cbind(sigma_u1, sigma_u2), s=s,
#'               delta=delta,
#'               distr = c(distr_marg1,distr_marg2,distr_cop))
#' dat<-data.frame(y1=Y[,1],y2=Y[,2], x1, x2, x3, x4, x5, x6, x7)
#' 
#' #Write formulae for parameters
#' mu_1_formula<-y1~s(x1,bs="ps")
#' sigma_v1_formula<-~s(x2,bs="ps")
#' sigma_u1_formula<-~1
#' mu_2_formula<-y2~s(x4,bs="ps")
#' sigma_v2_formula<-~1
#' sigma_u2_formula<-~s(x6,bs="ps")
#' delta_formula<-~s(x7,bs="ps")
#' 
#' #Fit model
#' model<-dsfa(formula=list(mu_1_formula,sigma_v1_formula,sigma_u1_formula,
#'                               mu_2_formula,sigma_v2_formula,sigma_u2_formula,
#'                               delta_formula),  data=dat,
#'                  family=comper_mv(s=s, distr=c(distr_marg1,distr_marg2,distr_cop)),
#'                  optimizer="efs")
#' 
#' #Model summary
#' summary(model)
#' 
#' #Smooth effects
#' #Effect of x1 on the predictor of the production function of margin 1
#' plot(model, select=1) #Estimated function
#' lines(x1[order(x1)], 2*x1[order(x1)]+(-1/3)*x1[order(x1)]^2-
#'         mean(2*x1+(-1/3)*x1^2), col=2) #True effect
#' 
#' #Effect of x2 on the predictor of the noise of margin 1
#' plot(model, select=2) #Estimated function
#' lines(x2[order(x2)], -1.5+sin(2*pi*x2[order(x2)])-
#'         mean(-1.5+sin(2*pi*x2)),col=2) #True effect
#' 
#' #Effect of x4 on the predictor of the production function of margin 2
#' plot(model, select=3) #Estimated function
#' lines(x4[order(x4)], 3+5*x4[order(x4)]^2+4*log(x4[order(x4)]+2)^(1/4)-
#'         mean(3+5*x4^2+4*log(x4+2)^(1/4)), col=2) #True effect
#' 
#' #Effect of x6 on the predictor of the inefficiency of margin 2
#' plot(model, select=4) #Estimated function
#' lines(x6[order(x6)], -1+sin(2*pi*x6[order(x6)])-
#'         mean(-1+sin(2*pi*x6)),col=2) #True effect
#' 
#' #Effect of x7 on the predictor of the copula
#' plot(model, select=5) #Estimated function
#' lines(x7[order(x7)], 2.5*cos(4*x7[order(x7)])-
#'         mean(2.5*cos(4*x7)),col=2) #True effect
#'
#' efficiency(model)
#' elasticity(model)
#' 
#' #' ### Second example with real data
#' 
#' data(BurkinaFarms)
#' data(BurkinaFarms_polys)
#' 
#' #Write formulae for parameters
#' mu_1_formula<-qharv_millet~s(land_millet, bs="ps")+s(labour_millet, bs="ps")+
#'                            s(material, bs="ps")+s(fert_millet, bs="ps")+
#'                            s(adm1, bs="mrf",xt=BurkinaFarms_polys)
#' sigma_v1_formula<-~1
#' sigma_u1_formula<-~farmtype+s(pest_millet, bs="ps")
#' 
#' mu_2_formula<-qharv_sorghum~s(land_sorghum, bs="ps")+s(labour_sorghum, bs="ps")+
#'                             s(material, bs="ps")+s(fert_sorghum, bs="ps")+
#'                             s(adm1, bs="mrf",xt=BurkinaFarms_polys)
#' sigma_v2_formula<-~1
#' sigma_u2_formula<-~farmtype+s(pest_sorghum, bs="ps")
#' 
#' delta_formula<-~1
#' 
#' model<-dsfa(formula=list(mu_1_formula, sigma_v1_formula, sigma_u1_formula,
#'                                 mu_2_formula, sigma_v2_formula, sigma_u2_formula,
#'                                 delta_formula),
#'                         data=BurkinaFarms,
#'                         family=comper_mv(s=c(-1,-1),
#'                         distr=c("normhnorm","normhnorm","normal")),
#'                         optimizer="efs")
#' plot(model) 
#'}
#'
#' @references
#' \itemize{
#' \item \insertRef{schmidt2023multivariate}{dsfa}
#' \item \insertRef{wood2017generalized}{dsfa}
#' \item \insertRef{aigner1977formulation}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' \item \insertRef{schmidt2020analytic}{dsfa}
#' }
#' @export
#comper distribution object for mgcv
comper_mv<- function(link = list("identity", "logshift", "logshift","identity", "logshift", "logshift","glogit"), s = c(-1,-1), distr=c("normhnorm","normhnorm","normal"), rot=0, b=1e-2){
  #Object for mgcv::gam such that the composed-error distribution can be estimated.
  
  #Number of parameters
  npar <- 7
  
  #Get boundaries of parameter space
  minmax<-delta_bounds(distr[3])
  
  #Link functions
  if (length(link) != npar) stop("comperr_mv requires 7 links specified as character strings")
  okLinks <- list("identity", "logshift", "logshift","identity", "logshift", "logshift", "glogit")
  stats <- list()
  param.names <- c("mu1", "sigma_v1", "sigma_u1","mu2", "sigma_v2", "sigma_u2","delta")
  
  for (i in 1:npar) { # Links for mu, sigma_v, sigma_u
    if (link[[i]] %in% okLinks[[i]]) {
      
      if(link[[i]]%in%c("identity", "log")){
        
        stats[[i]] <- stats::make.link(link[[i]])
        fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                              mu.eta=stats[[i]]$mu.eta),
                         class="family")
        fam <- mgcv::fix.family.link(fam)
        stats[[i]]$d2link <- fam$d2link
        stats[[i]]$d3link <- fam$d3link
        stats[[i]]$d4link <- fam$d4link
      }
      
      if(link[[i]]%in%c("logshift")){
        
        stats[[i]] <- list()
        stats[[i]]$valideta <- function(eta) TRUE 
        stats[[i]]$link = link[[i]]
        stats[[i]]$linkfun <- eval(parse(text=paste("function(mu) log(mu)  + ",b)))
        stats[[i]]$linkinv <- eval(parse(text=paste("function(eta) exp(eta - ",b,")")))
        stats[[i]]$mu.eta <-  eval(parse(text=paste("function(eta) exp(eta - ",b,")")))
        stats[[i]]$d2link <-  eval(parse(text=paste("function(mu)  -1/mu^2",sep='')))
        stats[[i]]$d3link <-  eval(parse(text=paste("function(mu)  2/mu^3",sep='')))
        stats[[i]]$d4link <-  eval(parse(text=paste("function(mu)  -6/mu^4",sep='')))
      } 
      
      if(link[[i]]%in%c("glogit")){
        
        stats[[i]] <- list()
        stats[[i]]$valideta <- function(eta) TRUE
        stats[[i]]$link = link[[i]]
        
        stats[[i]]$linkfun <- eval(parse(text = paste("function(mu) log((-",minmax[1],"+mu)/(",minmax[2],"-mu))")))#eval(parse(text = paste("function(mu) log((mu+1)/(1-mu))")))
        stats[[i]]$linkinv <- eval(parse(text = paste("function(eta) exp(eta)/(1+exp(eta))*(",minmax[2],"-",minmax[1],")+",minmax[1])))
        stats[[i]]$mu.eta <- eval(parse(text = paste("function(eta) exp(eta)*(",minmax[2],"-",minmax[1],")/(exp(eta)+1)^2")))
        stats[[i]]$d2link <- eval(parse(text = paste("function(mu) 1/(",minmax[2],"-mu)^2-1/(",minmax[1],"-mu)^2")))#eval(parse(text = paste("function(mu) 4*mu/(mu^2-1)^2"))) #-((2 * exp(eta) * (-1 + exp(eta)))/(1 + exp(eta))^3)
        stats[[i]]$d3link <- eval(parse(text = paste("function(mu) 2*(1/(",minmax[2],"-mu)^3+1/(-",minmax[1],"+mu)^3)"))) #(2 * exp(eta) * (1 - 4 * exp(eta) + exp(2 * eta)))/(1 + exp(eta))^4
        stats[[i]]$d4link <- eval(parse(text = paste("function(mu) 6/(",minmax[2],"-mu)^4-6/(",minmax[1],"-mu)^4"))) #-((2 * exp(eta)*  (-1 + 11 * exp(eta) - 11 * exp(2 * eta) + exp(3 * eta)))/(1 + exp(eta))^5)
        
      } 
    } else {
      stop(link[[i]]," link not available for ", param.names[i]," parameter of comper")
    }
  }
  
  # for (i in 1:6) { # Links for "mu1", "sigma_v1", "sigma_u1","mu2", "sigma_v2", "sigma_u2","delta"
  #   if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- stats::make.link(link[[i]]) else
  #     stop(link[[i]]," link not available for ", param.names[i]," parameter of comper_mv")
  #   fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
  #                         mu.eta=stats[[i]]$mu.eta),
  #                    class="family")
  #   fam <- mgcv::fix.family.link(fam)
  #   stats[[i]]$d2link <- fam$d2link
  #   stats[[i]]$d3link <- fam$d3link
  #   stats[[i]]$d4link <- fam$d4link
  # }
  # if (link[[7]] %in% okLinks[[7]]) {
  #   stats[[7]] <- list()
  #   stats[[7]]$valideta <- function(eta) TRUE
  #   stats[[7]]$link = link[[7]]
  #   
  #   stats[[7]]$linkfun <- eval(parse(text = paste("function(mu) log((-",minmax[1],"+mu)/(",minmax[2],"-mu))")))#eval(parse(text = paste("function(mu) log((mu+1)/(1-mu))")))
  #   stats[[7]]$linkinv <- eval(parse(text = paste("function(eta) exp(eta)/(1+exp(eta))*(",minmax[2],"-",minmax[1],")+",minmax[1])))
  #   stats[[7]]$mu.eta <- eval(parse(text = paste("function(eta) exp(eta)*(",minmax[2],"-",minmax[1],")/(exp(eta)+1)^2")))
  #   stats[[7]]$d2link <- eval(parse(text = paste("function(mu) 1/(",minmax[2],"-mu)^2-1/(",minmax[1],"-mu)^2")))#eval(parse(text = paste("function(mu) 4*mu/(mu^2-1)^2"))) #-((2 * exp(eta) * (-1 + exp(eta)))/(1 + exp(eta))^3)
  #   stats[[7]]$d3link <- eval(parse(text = paste("function(mu) 2*(1/(",minmax[2],"-mu)^3+1/(-",minmax[1],"+mu)^3)"))) #(2 * exp(eta) * (1 - 4 * exp(eta) + exp(2 * eta)))/(1 + exp(eta))^4
  #   stats[[7]]$d4link <- eval(parse(text = paste("function(mu) 6/(",minmax[2],"-mu)^4-6/(",minmax[1],"-mu)^4"))) #-((2 * exp(eta)*  (-1 + 11 * exp(eta) - 11 * exp(2 * eta) + exp(3 * eta)))/(1 + exp(eta))^5)
  # }
  # else stop(link[[7]], " link not available comper_mv distribution")
  # 
  residuals <- function(object, type = c("deviance", "response", "normalized")) {
    #Calculate residuals
    type <- match.arg(type)
    
    if(type%in%c("deviance", "response")){
      mom1<-par2mom(mu=object$fitted[,1], sigma_v = object$fitted[,2], sigma_u=object$fitted[,3], s=object$family$s[1], distr=object$family$distr[1])
      mom2<-par2mom(mu=object$fitted[,4], sigma_v = object$fitted[,5], sigma_u=object$fitted[,6], s=object$family$s[2], distr=object$family$distr[2])
      
      y_hat<-cbind(mom1[,1],mom2[,1])
      
      rsd <- object$y-y_hat
      
      if (type=="response"){
        out<-rsd
      }  else {
        out<-rsd/cbind(mom1[,2],mom2[,2])
      }
    } else{
      out<-cbind(stats::qnorm(pcomper(q=object$y[,1], mu=object$fitted[,1], sigma_v = object$fitted[,2], sigma_u=object$fitted[,3], s=object$family$s[1], distr=object$family$distr[1])),
                 stats::qnorm(pcomper(q=object$y[,2], mu=object$fitted[,4], sigma_v = object$fitted[,5], sigma_u=object$fitted[,6], s=object$family$s[2], distr=object$family$distr[2])))
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
    npar <- 7
    n <- nrow(y)
    
    #Get additive predictors
    eta<-matrix(0,n,npar)
    for(i in 1:npar){
      eta[,i]<-drop( X[ ,jj[[i]], drop=FALSE] %*% coef[jj[[i]]] )
    }
    
    #Additive predictors 2 parameters
    theta<-matrix(0,n,npar)
    for(i in 1:npar){
      theta[,i]<-family$linfo[[i]]$linkinv(eta[,i])
    }
    
    #Evaluate ll
    l0<-dcomper_mv_cpp(x=y, m=theta[,c(1,4), drop=FALSE], v=theta[,c(2,5), drop=FALSE], u=theta[,c(3,6), drop=FALSE], delta=theta[,7, drop=TRUE], s=family$s, distr=family$distr, rot=family$rot, deriv_order=2, tri=family$tri_mat, logp = TRUE)
    
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
    #Idea is to get starting values by estimation of the margins independently
    #then probability transform observations and estimate via cop.
    
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
      # start <- c()
      
      #Margin 1
      y1 <- y[,1]
      x1 <- x[ , jj[[1]], drop=FALSE]
      x2 <- x[ , jj[[2]], drop=FALSE]
      x3 <- x[ , jj[[3]], drop=FALSE]
      
      m1<-try(mgcv::gam(formula=list(y1 ~ x1 - 1, ~x2 - 1, ~x3 - 1),
                        family=comper(s=family$s[1], distr=family$distr[1]), optimizer = c("efs")), silent = TRUE)
      
      if (any(class(m1) == "try-error")) {
        m1<-try(mgcv::gam(formula=list(y1 ~ x1 - 1, ~x2 - 1, ~x3 - 1),
                          family=comper(s=family$s[1], distr=family$distr[1]), optimizer = c("outer","newton")), silent = TRUE)
      }

      if (any(class(m1) == "try-error")) {
        stop(paste("Algorithm to fit starting values for marginal 1 did not converge", "\n", ""))
      }
      
      #Margin 2
      y2 <- y[,2]
      x4 <- x[ , jj[[4]], drop=FALSE]
      x5 <- x[ , jj[[5]], drop=FALSE]
      x6 <- x[ , jj[[6]], drop=FALSE]
      
      m2<-try(mgcv::gam(formula=list(y2 ~ x4 - 1, ~x5 - 1, ~x6 - 1),
                        family=comper(s=family$s[2], distr=family$distr[2]), optimizer = c("efs")), silent = TRUE)
      
      if (any(class(m2) == "try-error")) {
        m2<-try(mgcv::gam(formula=list(y2 ~ x4 - 1, ~x5 - 1, ~x6 - 1),
                        family=comper(s=family$s[2], distr=family$distr[2]), optimizer = c("outer","newton")), silent = TRUE)
      }

      if (any(class(m2) == "try-error")) {
        stop(paste("Algorithm to fit starting values for marginal 2 did not converge", "\n", ""))
      }
      
      #Probability integral transform
      F1<-pcomper(q=m1$y, mu=m1$fitted.values[,1], sigma_v=m1$fitted.values[,2], sigma_u=m1$fitted.values[,3], s=family$s[1], distr=family$distr[1], log.p=FALSE)
      F2<-pcomper(q=m2$y, mu=m2$fitted.values[,1], sigma_v=m2$fitted.values[,2], sigma_u=m2$fitted.values[,3], s=family$s[2], distr=family$distr[2], log.p=FALSE)
      
      #Correct for numerical inaccuracy
      F1[F1>=1]<-1-1e-16
      F1[F1<=0]<-0+1e-16
      
      F2[F2>=1]<-1-1e-16
      F2[F2<=0]<-0+1e-16
      
      #Fit copula model with pseudo observations
      x7<-x[,jj[[7]],drop=FALSE]
                
      m3<-try(mgcv::gam(y1 ~ x7-1, family=cop(W=cbind(F1,F2),
              distr=family$distr[3], rot=family$rot), optimizer = "efs"), silent = TRUE)
      
      if (any(class(m3) == "try-error")) {
        stop(paste("Algorithm to fit starting values for copula did not converge", "\n", ""))
      }
      
      start<-c(m1$coefficients, m2$coefficients, m3$coefficients)
    }
  }) 
  
  rd <- function(mu, wt, scale) {
    #Random number generation for comper
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    
    return(rcomper_mv(n=nrow(mu), mu = mu[,c(1,4)], sigma_v = mu[,c(2,5)], sigma_u = mu[,c(3,6)], delta = mu[,7], s=attr(mu,"s"), distr=attr(mu,"distr"), rot=attr(mu,"rot")))
  }
  
  cdf <- function(q, mu, wt, scale) {
    #Cumulative distribution function of comper
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }
    
    return(pcomper_mv(q=q, mu = mu[,c(1,4)], sigma_v = mu[,c(2,5)], sigma_u = mu[,c(3,6)], delta = mu[,7], s=attr(mu,"s"), distr=attr(mu,"distr"), rot=attr(mu,"rot")))
  }
  
  predict <- function(family, se=FALSE, eta=NULL,y=NULL,X=NULL,beta=NULL,off=NULL,Vb=NULL) {
    #Prediction function
    # optional function to give predicted values - idea is that
    # predict.gam(...,type="response") will use this, and that
    # either eta will be provided, or {X, beta, off, Vb}. family$data
    # contains any family specific extra information.
    # if se = FALSE returns one item list containing matrix otherwise
    # list of two matrices "fit" and "se.fit"...
    
    #Extract linear predictor index
    jj <- attr(X, "lpi") 
    
    #Number of parameters and observations
    npar <- 7
    n <- nrow(y)
    
    if (is.null(eta)) {
      if (is.null(off)){
        off <- list(0,0,0,0,0,0,0)
      } 
      off[[8]] <- 0
      
      for (i in 1:npar) {
        if (is.null(off[[i]])){
          off[[i]] <- 0
        } 
      } 
      
      #Extract linear predictor index
      jj <- attr(X,"lpi")
      
      #Calculate additive predictors for mu, sigma_v, sigma_u
      eta<-matrix(0, nrow(X), npar)
      for(i in 1:npar){
        eta[,i]<-drop(X[,jj[[i]],drop=FALSE]%*%beta[jj[[i]]] + off[[i]])
      }
      
      if (se) {
        stop("se still available for this family")
      }
      
    } else {
      se <- FALSE
    }
    
    #Additive predictors 2 parameters
    theta<-matrix(0,nrow(X),npar)
    for(i in 1:npar){
      theta[,i]<-family$linfo[[i]]$linkinv(eta[,i])
    }
    
    #Calculate moments from parameters
    mom1<-par2mom(mu = theta[,1], sigma_v = theta[,2], sigma_u = theta[,3], s=family$s[1], distr=family$distr[1])
    mom2<-par2mom(mu = theta[,4], sigma_v = theta[,5], sigma_u = theta[,6], s=family$s[2], distr=family$distr[2])
    
    #Assign mean
    fv <- list(cbind(mom1[,1],mom2[,1]))
    if (!se) return(fv) else {
      stop("se not still available for this family")
      #Assign mean and standard deviation
      fv <- list(fit=cbind(mom1[,1], mom2[,1]), se.fit=cbind(mom1[,2], mom2[,2]))
      return(fv)
    }
  }
  
  # preinitialize <- function(G) {
  #   G$family$formula<-G$formula
  #   G$family$data<-G$mf
  #   
  #   return(list(family=G$family))
  # }
  
  postproc <- expression({
    # object$family$family<-NULL
    # object$family$data<-NULL
    
    attr(object$fitted.values,"s")<-object$family$s
    attr(object$fitted.values,"distr")<-object$family$distr
    attr(object$fitted.values,"rot")<-object$family$rot
    
    object$fitted.values
  })
  
  structure(list(family="comper_mv", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), # symmetric indices for accessing derivative arrays
                 tri_mat = list(trind_generator(3),trind_generator(3),trind_generator(1),trind_generator(6),trind_generator(7)), # symmetric indices for accessing derivative matrices
                 initialize = initialize, #initial parameters
                 s = s, #production or cost function
                 distr = distr, #specifiying distribution
                 rot = rot, #rotation
                 b=b,
                 residuals = residuals, #residual function
                 rd=rd, #random number generation
                 cdf=cdf, #cdf function
                 predict=predict, #prediction function for mgcv
                 #preinitialize=preinitialize,
                 postproc=postproc, #Assigning attributes such that other functions work
                 linfo = stats, # link information list
                 d2link=1, d3link=1, d4link=1, # signals to fix.family.link that all done
                 ls=1, # signals that ls not needed here
                 available.derivs = 2), # can use full Newton here
            class = c("general.family","extended.family","family"))
}

