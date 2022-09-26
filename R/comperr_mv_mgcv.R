#' Composed error multivariate distribution object for mgcv
#'
#' The comperr_mv family implements the composed error multivariate distribution in which the \eqn{\mu_1}, \eqn{\sigma_{V1}}, \eqn{\sigma_{U1}}, (or \eqn{\lambda_1}), \eqn{\mu_2}, \eqn{\sigma_{V2}}, \eqn{\sigma_{U2}}, (or \eqn{\lambda_2}) and \eqn{\delta} can depend on additive predictors. Useable only with \code{\link[mgcv:gam]{gam()}}, the additive predictors are specified via a list of formulae.
#'
#'
#' @param link seven item list specifying the link for the \eqn{\mu_1}, \eqn{\sigma_{V1}}, \eqn{\sigma_{U1}}, \eqn{\lambda_1}, \eqn{\mu_2}, \eqn{\sigma_{V2}}, \eqn{\sigma_{U2}}, \eqn{\lambda_2} and \eqn{\delta} parameters. See details.
#' @param s1 \eqn{s_1=-1} for production and \eqn{s_1=1} for cost function for margin 1.
#' @param s2 \eqn{s_2=-1} for production and \eqn{s_2=1} for cost function for margin 2.
#' @param family_mv vector of length three, specifying the bivariate distribution.  First element is the name of the first marginal distribution. Second element is the name of the second marginal distribution. Third element specifies the copula. See \code{dcop()} for more details.
#'
#' @return An object inheriting from class \code{general.family} of the 'mgcv' package, which can be used in the dsfa package.
#'
#' @details Used with gam to fit distributional stochastic frontier model. The function \code{gam} is from the mgcv package is called with a list containing seven formulae:
#' \enumerate{
#'   \item The first formula specifies the response of margin 1 on the left hand side and the structure of the additive predictor for \eqn{\mu_1} parameter on the right hand side. Link function is "identity".
#'   \item The second formula is one sided, specifying the additive predictor for the  \eqn{\sigma_{V1}} on the right hand side. Link function is "log".
#'   \item The third formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_{U1}} or \eqn{\lambda_1} on the right hand side. Link function is "log".
#'   \item The fourth formula specifies the response of margin 2 on the left hand side and the structure of the additive predictor for \eqn{\mu_2} parameter on the right hand side. Link function is "identity".
#'   \item The fifth formula is one sided, specifying the additive predictor for the  \eqn{\sigma_{V2}} on the right hand side. Link function is "log".
#'   \item The sixth formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_{U2}} or \eqn{\lambda_2} on the right hand side. Link function is "log".
#'   \item The seventh formula  is one sided, specifying the additive predictor for the  \eqn{\delta} on the right hand side. Link function is "glogit".
#' }
#' The fitted values and linear predictors for this family will be seven column matrices. The columns correspond with the order of the formulae for the parameters.
#'
#' @examples
#' \donttest{
#' #Set seed, sample size and type of function
#' set.seed(1337)
#' N=1000 #Sample size
#' s1<--1 #Set to production function for margin 1
#' s2<-1 #Set to cost function for margin 2
#'
#' #Generate covariates
#' x1<-runif(N,-1,1); x2<-runif(N,-1,1); x3<-runif(N,-1,1)
#' x4<-runif(N,-1,1); x5<-runif(N,-1,1); x6<-runif(N,-1,1)
#' x7<-runif(N,-1,1)
#'
#' mu1=4+x1 #production function parameter 1
#' sigma_v1=exp(-1.5+0.75*x2) #noise parameter 1
#' sigma_u1=exp(-1+1.25*x3) #inefficiency parameter 1
#' mu2=3+2*x4 #cost function parameter 2
#' sigma_v2=exp(-1.5+0.75*x5) #noise parameter 2
#' sigma_u2=exp(-1+.75*x6) #inefficiency parameter 2
#' delta<-(exp(1+2.5*x7)-1)/(exp(1+2.5*x7)+1) #delta
#'
#' #Simulate responses and create dataset
#' Y<-rcomperr_mv(n=N,
#'                mu1=mu1, sigma_v1=sigma_v1, par_u1 = sigma_u1, s1=s1,
#'                mu2=mu2, sigma_v2=sigma_v2, par_u2 = sigma_u2, s2=s2,
#'                delta=delta, family=c("normhnorm","normhnorm","normal"))
#' dat<-data.frame(y1=Y[,1],y2=Y[,2], x1, x2, x3, x4, x5, x6, x7)
#'
#' #Write formulae for parameters
#' mu_1_formula<-y1~x1
#' sigma_v1_formula<-~x2
#' sigma_u1_formula<-~x3
#' mu_2_formula<-y2~x4
#' sigma_v2_formula<-~x5
#' sigma_u2_formula<-~x6
#' delta_formula<-~x7
#'
#' #Fit model
#' model<-mgcv::gam(formula=list(mu_1_formula,sigma_v1_formula,sigma_u1_formula,
#'                               mu_2_formula,sigma_v2_formula,sigma_u2_formula,
#'                               delta_formula),
#'                  data=dat,
#'                  family=comperr_mv(s1=s1, s2=s2, family_mv=c("normhnorm","normhnorm","normal")),
#'                  optimizer="efs")
#'
#' #Model summary
#' summary(model)
#' }
#' @references
#' \itemize{
#' \item \insertRef{schmidt2022mvdsfm}{dsfa}
#' \item \insertRef{wood2017generalized}{dsfa}
#' }
#' @export
comperr_mv <- function(link = list("identity", "log", "log", "identity", "log", "log", "glogit"), s1 = -1, s2 = -1, family_mv=c("normhnorm","normhnorm","normal")){
  #Object for mgcv::gam such that the composed error multivariate distribution can be estimated.

  #Number of parameters
  npar <- 7

  #Delta parameter space boundaries
  if(family_mv[3]=="independent"){
    a<-0
    b<-1
  }

  if(family_mv[3]=="normal"){
    a<--1
    b<-1
  }

  if(family_mv[3]=="clayton"){
    a<-0+1e-16
    b<-28
  }

  if(family_mv[3]=="gumbel"){
    a<-1
    b<-17
  }

  if(family_mv[3]=="frank"){
    a<--35
    b<-35
  }

  if(family_mv[3]=="joe"){
    a<-1+1e-16
    b<-30
  }

  #Link functions
  if (length(link) != npar) stop("comperr_mv requires 7 links specified as character strings")
  okLinks <- list("identity", "log", "log","identity", "log", "log", "glogit")
  stats <- list()
  param.names <- c("mu1", "sigma_v1", "par_u1","mu2", "sigma_v2", "par_u2","delta")
  for (i in 1:6) { # Links for "mu1", "sigma_v1", "sigma_u1","mu2", "sigma_v2", "sigma_u2","delta"
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- stats::make.link(link[[i]]) else
      stop(link[[i]]," link not available for ", param.names[i]," parameter of comperr_mv")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- mgcv::fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  if (link[[7]] %in% okLinks[[7]]) {
    stats[[7]] <- list()
    stats[[7]]$valideta <- function(eta) TRUE
    stats[[7]]$link = link[[7]]

    stats[[7]]$linkfun <- eval(parse(text = paste("function(mu) log((-",a,"+mu)/(",b,"-mu))")))#eval(parse(text = paste("function(mu) log((mu+1)/(1-mu))")))
    stats[[7]]$linkinv <- eval(parse(text = paste("function(eta) exp(eta)/(1+exp(eta))*(",b,"-",a,")+",a)))
    stats[[7]]$mu.eta <- eval(parse(text = paste("function(eta) exp(eta)*(",b,"-",a,")/(exp(eta)+1)^2")))
    stats[[7]]$d2link <- eval(parse(text = paste("function(mu) 1/(",b,"-mu)^2-1/(",a,"-mu)^2")))#eval(parse(text = paste("function(mu) 4*mu/(mu^2-1)^2"))) #-((2 * exp(eta) * (-1 + exp(eta)))/(1 + exp(eta))^3)
    stats[[7]]$d3link <- eval(parse(text = paste("function(mu) 2*(1/(",b,"-mu)^3+1/(-",a,"+mu)^3)"))) #(2 * exp(eta) * (1 - 4 * exp(eta) + exp(2 * eta)))/(1 + exp(eta))^4
    stats[[7]]$d4link <- eval(parse(text = paste("function(mu) 6/(",b,"-mu)^4-6/(",a,"-mu)^4"))) #-((2 * exp(eta)*  (-1 + 11 * exp(eta) - 11 * exp(2 * eta) + exp(3 * eta)))/(1 + exp(eta))^5)
  }
  else stop(link[[7]], " link not available comperr_mv distribution")

  #Calculate residuals
  residuals <- function(object, type = c("deviance", "response")) {

    type <- match.arg(type)
    para1<-reparametrize(mu = object$fitted[,1], sigma_v = object$fitted[,2], par_u = object$fitted[,3],
                         s = object$family$s1, family = object$family$family_mv[1])
    para2<-reparametrize(mu = object$fitted[,4], sigma_v = object$fitted[,5], par_u = object$fitted[,6],
                         s = object$family$s2, family = object$family$family_mv[2])


    rsd <- object$y-cbind(para1$mean,para2$mean)

    if (type=="response") return(rsd) else
      return(rsd/cbind(para1$sd,para2$sd))
  }

  #Loglike function with derivatives
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv=0, d1b=0, d2b=0, Hp=NULL, rank=0, fh=NULL, D=NULL) {
    ## function defining the normhnorm model log lik.
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.


    # If offset is not null or a vector of zeros, give an error
    if( !is.null(offset[[1]]) && sum(abs(offset)) )  stop("offset not still available for this family")

    jj <- attr(X, "lpi") ## extract linear predictor index

    npar <- 7
    n <- nrow(y)

    eta <-  drop( X[ , jj[[1]], drop=FALSE] %*% coef[jj[[1]]] )
    eta1 <- drop( X[ , jj[[2]], drop=FALSE] %*% coef[jj[[2]]] )
    eta2 <- drop( X[ , jj[[3]], drop=FALSE] %*% coef[jj[[3]]] )
    eta3 <- drop( X[ , jj[[4]], drop=FALSE] %*% coef[jj[[4]]] )
    eta4 <- drop( X[ , jj[[5]], drop=FALSE] %*% coef[jj[[5]]] )
    eta5 <- drop( X[ , jj[[6]], drop=FALSE] %*% coef[jj[[6]]] )
    eta6 <- drop( X[ , jj[[7]], drop=FALSE] %*% coef[jj[[7]]] )

    mu1      <- family$linfo[[1]]$linkinv( eta )
    sigma_v1 <- family$linfo[[2]]$linkinv( eta1 )
    par_u1   <- family$linfo[[3]]$linkinv( eta2 )
    s1       <- family$s1
    mu2      <- family$linfo[[4]]$linkinv( eta3 )
    sigma_v2 <- family$linfo[[5]]$linkinv( eta4 )
    par_u2   <- family$linfo[[6]]$linkinv( eta5 )
    s2       <- family$s2
    delta      <- matrix(family$linfo[[7]]$linkinv( eta6 ), nrow=n)

    ##Define parameters
    l0<-dcomperr_mv(x1=y[,1], mu1=mu1, sigma_v1=sigma_v1, par_u1=par_u1, s1=s1,
                    x2=y[,2], mu2=mu2, sigma_v2=sigma_v2, par_u2=par_u2, s2=s2,
                    delta=delta, family_mv=family$family_mv, deriv=2, tri=family$tri, log.p = TRUE)

    l<-sum(l0)

    if (deriv>0) {
      ## the first derivatives

      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2),
                   family$linfo[[4]]$mu.eta(eta3),
                   family$linfo[[5]]$mu.eta(eta4),
                   family$linfo[[6]]$mu.eta(eta5),
                   family$linfo[[7]]$mu.eta(eta6))

      g2 <- cbind(family$linfo[[1]]$d2link(mu1),
                  family$linfo[[2]]$d2link(sigma_v1),
                  family$linfo[[3]]$d2link(par_u1),
                  family$linfo[[4]]$d2link(mu2),
                  family$linfo[[5]]$d2link(sigma_v2),
                  family$linfo[[6]]$d2link(par_u2),
                  family$linfo[[7]]$d2link(delta))
    }

    l1<-attr(l0,"gradient")
    l2<-attr(l0,"hessian")

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      stop("deriv>2 not implemented")
      g3 <- cbind(family$linfo[[1]]$d3link(mu1),
                  family$linfo[[2]]$d3link(sigma_v1),
                  family$linfo[[3]]$d3link(par_u1),
                  family$linfo[[4]]$d3link(mu2),
                  family$linfo[[5]]$d3link(sigma_v2),
                  family$linfo[[6]]$d3link(par_u2),
                  family$linfo[[7]]$d3link(delta))

      l3<-attr(l0,"l3")

    }

    if (deriv>3) {
      g4 <- cbind(family$linfo[[1]]$d4link(mu1),
                  family$linfo[[2]]$d4link(sigma_v1),
                  family$linfo[[3]]$d4link(par_u1),
                  family$linfo[[4]]$d4link(mu2),
                  family$linfo[[5]]$d4link(sigma_v2),
                  family$linfo[[6]]$d4link(par_u2),
                  family$linfo[[7]]$d4link(delta))

      l4<-attr(l0,"l4")
    }

    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4

      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- mgcv::gamlss.etamu(l1,l2,l3,l4,ig1,g2,g3,g4,i2,i3,i4,deriv-1)

      ## get the gradient and Hessian...
      ret <- mgcv::gamlss.gH(X,jj,de$l1,de$l2,i2,l3=de$l3,i3=i3,l4=de$l4,i4=i4,
                             d1b=d1b,d2b=d2b,deriv=deriv-1,fh=fh,D=D)

    } else ret <- list()
    ret$l <- l; ret
  } ## end ll

  #Function to calculate starting values of the coefficients
  initialize <- expression({
    ## idea is to get starting values utilizing the method of moments,

    n <- rep(1, nobs)
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x,"lpi")
      start <- rep(0,ncol(x))

      #Margin 1
      if(family$family_mv[1]=="normhnorm"){
        m1<-mgcv::gam(formula=list(family$formula[[1]], family$formula[[2]], family$formula[[3]]),
                      data=family$data, family=normhnorm(s=family$s1), optimizer = c("efs"))
      }
      if(family$family_mv[1]=="normexp"){
        m1<-mgcv::gam(formula=list(family$formula[[1]], family$formula[[2]], family$formula[[3]]),
                      data=family$data, family=normexp(s=family$s1), optimizer = c("efs"))
      }

      #Margin 2
      if(family$family_mv[2]=="normhnorm"){
        m2<-mgcv::gam(formula=list(family$formula[[4]], family$formula[[5]], family$formula[[6]]),
                      data=family$data, family=normhnorm(s=family$s2), optimizer = c("efs"))
      }
      if(family$family_mv[2]=="normexp"){
        m2<-mgcv::gam(formula=list(family$formula[[4]], family$formula[[5]], family$formula[[6]]),
                      data=family$data, family=normexp(s=family$s2), optimizer = c("efs"))
      }

      start[-c(jj[[7]])]<-c(m1$coefficients,m2$coefficients)

      u1<-pcomperr(q=y[,1], mu=m1$fitted.values[,1], sigma_v=m1$fitted.values[,2], par_u=m1$fitted.values[,3], s=family$s1, family=family$family_mv[1], log.p=FALSE)
      u2<-pcomperr(q=y[,2], mu=m2$fitted.values[,1], sigma_v=m2$fitted.values[,2], par_u=m2$fitted.values[,3], s=family$s2, family=family$family_mv[2], log.p=FALSE)

      family$data$y3<-rep(1,nrow(Y))
      adj_par_formula<-update.formula(family$formula[[7]], y3 ~ .)
      m3<-mgcv::gam(list(adj_par_formula),
                    data=family$data, family=cop(W=as.matrix(cbind(u1,u2)), family_cop=family$family_mv[3]), optimizer = c("efs"))

      start[jj[[7]]] <- m3$coefficients
    }
  }) ## initialize

  #Random number generation for comperr_mv
  rd <- function(mu, wt, scale, family) {
    ## random number generation
    mu <- as.matrix(mu)
    if(ncol(mu)==2){ mu <- t(mu) }

    return(rcomperr_mv(ncol(mu), mu1 = mu[,1], sigma_v1 = mu[,2], par_u1 = mu[,3], s1=family$s1,
                                 mu2 = mu[,4], sigma_v2 = mu[,5], par_u2 = mu[,6], s2=family$s2,
                                 delta= mu[,7], family_mv=family$family_mv))
  } ## random number generation

  #Cumulative distribution function of comperr_mv
  cdf <- function(q, mu, wt, scale, logp, family) {
    ##cumulative distribution function
    mu <- as.matrix(mu)
    if(ncol(mu)==2){ mu <- t(mu) }

    return(pcomperr_mv(q1=q[,1], mu1 = mu[,1], sigma_v1 = mu[,2], par_u1 = mu[,3], s1=family$s1,
                       q2=q[,2], mu2 = mu[,4], sigma_v2 = mu[,5], par_u2 = mu[,6], s2=family$s2,
                       delta=mu[,7], family_mv=family$family_mv),
                       log.p = logp)
  } ##cumulative distribution function

  #Prediction function
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.
    ## if se = FALSE returns one item list containing matrix otherwise
    ## list of two matrices "fit" and "se.fit"...

    if (is.null(eta)) {
      if (is.null(off)) off <- list(0,0,0)
      off[[4]] <- 0
      for (i in 1:7) if (is.null(off[[i]])) off[[i]] <- 0
      lpi <- attr(X,"lpi")
      X1 <- X[,lpi[[1]],drop=FALSE]
      X2 <- X[,lpi[[2]],drop=FALSE]
      X3 <- X[,lpi[[3]],drop=FALSE]
      X4 <- X[,lpi[[4]],drop=FALSE]
      X5 <- X[,lpi[[5]],drop=FALSE]
      X6 <- X[,lpi[[6]],drop=FALSE]
      X7 <- X[,lpi[[7]],drop=FALSE]

      mu1      <- drop(X1%*%beta[lpi[[1]]] + off[[1]]) ## linear predictor for mu
      sigma_v1 <- drop(X2%*%beta[lpi[[2]]] + off[[2]])  ## linear predictor for sigma_v parameter
      par_u1 <- drop(X3%*%beta[lpi[[3]]] + off[[3]]) ## linear predictor for sigma_u parameter
      mu2      <- drop(X4%*%beta[lpi[[4]]] + off[[4]]) ## linear predictor for mu
      sigma_v2 <- drop(X5%*%beta[lpi[[5]]] + off[[5]])  ## linear predictor for sigma_v parameter
      par_u2 <- drop(X6%*%beta[lpi[[6]]] + off[[6]]) ## linear predictor for sigma_u parameter

      if (se) {
        stop("se still available for this family")
      }
    } else {
      se <- FALSE
      mu1 <- eta[,1]
      sigma_v1 <- eta[,2]
      par_u1 <- eta[,3]
      mu2 <- eta[,4]
      sigma_v2 <- eta[,5]
      par_u2 <- eta[,6]
    }


    para1<-reparametrize(mu = mu1, sigma_v = sigma_v1, par_u = par_u1, s=family$s1, family=family$family_mv[1])
    para2<-reparametrize(mu = mu2, sigma_v = sigma_v2, par_u = par_u2, s=family$s2, family=family$family_mv[2])


    fv <- list(cbind(para1$mean,para2$mean))
    if (!se) return(fv) else {
      stop("se not still available for this family")
      fv <- list(fit=cbind(para1$mean,para2$mean), se.fit=cbind(para1$sd,para2$sd))
      return(fv)
    }
  } ##predict

  preinitialize <- function(G) {
    G$family$formula<-G$formula
    G$family$data<-G$mf

    return(list(family=G$family))
  }

  structure(list(family="comperr_mv", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 preinitialize=preinitialize,
                 s1 = s1,
                 s2 = s2,
                 family_mv = family_mv,
                 residuals=residuals,
                 rd=rd,
                 cdf=cdf,
                 predict=predict,
                 linfo = stats, ## link information list
                 d2link=1, d3link=1, d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## comperr_mv

