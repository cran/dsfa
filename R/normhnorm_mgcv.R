#' normhnorm family
#'
#' The normhnorm family implements the normal-halfnormal distribution in which the \eqn{\mu}, \eqn{\sigma_V} and \eqn{\sigma_U} can depend on additive predictors. Useable only with \code{mgcv::gam}, the additive predictors are specified via a list of formulae.
#'
#'
#' @param link three item list specifying the link for the \eqn{\mu}, \eqn{\sigma_V} and \eqn{\sigma_U} parameters. See details.
#' @param s \eqn{s=-1} for production and \eqn{s=1} for cost function.
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the dsfa package.
#'
#' @details Used with gam to fit distributional stochastic frontier model. The function \code{gam} is from the mgcv package is called with a list containing three formulae:
#' \enumerate{
#'   \item The first formula specifies the response on the left hand side and the structure of the additive predictor for \eqn{\mu} parameter on the right hand side. Link function is "identity".
#'   \item The second formula is one sided, specifying the additive predictor for the  \eqn{\sigma_V} on the right hand side. Link function is "log".
#'   \item The third formula  is one sided, specifying the additive predictor for the  \eqn{\sigma_U} on the right hand side. Link function is "log".
#' }
#' The fitted values and linear predictors for this family will be three column matrices. The first column is the \eqn{\mu}, the second column is the \eqn{\sigma_V}, the third column is \eqn{\sigma_U}.
#'
#' @examples
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
#' sigma_V_formula<-~1+x4
#' sigma_U_formula<-~1+s(x5, bs="ps")
#'
#' #Fit model
#' model<-mgcv::gam(formula=list(mu_formula, sigma_V_formula, sigma_U_formula),
#'                  data=dat, family=normhnorm(s=s), optimizer = c("efs"))
#'
#' #Model summary
#' summary(model)
#'
#' #Smooth effects
#' #Effect of x3 on the predictor of the production function
#' plot(model, select=1) #Estimated function
#' lines(x3[order(x3)], 6*log(x3[order(x3)]+2)^(1/4)-
#'       mean(6*log(x3[order(x3)]+2)^(1/4)), col=2) #True effect
#'
#' #Effect of x5 on the predictor of the inefficiency
#' plot(model, select=2) #Estimated function
#' lines(x5[order(x5)], -1+sin(2*pi*x5)[order(x5)]-
#'       mean(-1+sin(2*pi*x5)),col=2) #True effect
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
#'
#normhnorm distribution object for mgcv
normhnorm <- function(link = list("identity", "log", "log"), s = -1){
  #Object for mgcv::gam such that the normal-halfnormal distribution can be estimated.

  #Number of parameters
  npar <- 3

  #Link functions
  if (length(link) != npar) stop("normhnorm requires 3 links specified as character strings")
  okLinks <- list("identity", "log", "log")
  stats <- list()
  param.names <- c("mu", "sigma_v", "sigma_u")
  for (i in 1:npar) { # Links for mu, sigma_v, sigma_u
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- stats::make.link(link[[i]]) else
      stop(link[[i]]," link not available for ", param.names[i]," parameter of normhnorm")
    fam <- structure(list(link=link[[i]],canonical="none",linkfun=stats[[i]]$linkfun,
                          mu.eta=stats[[i]]$mu.eta),
                     class="family")
    fam <- mgcv::fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }

  #Calculate residuals
  residuals <- function(object, type = c("deviance", "response")) {

    type <- match.arg(type)
    mu<-object$fitted[,1]
    sigma_v<-object$fitted[,2]
    sigma_u<-object$fitted[,3]

    para<-reparametrize(mu = mu, sigma_v = sigma_v, sigma_u = sigma_u, s=object$family$s)

    y_hat<-para$mean

    rsd <- object$y-y_hat

    if (type=="response") return(rsd) else
      return(rsd/para$sd)
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

    npar <- 3
    n <- length(y)

    eta <-  drop( X[ , jj[[1]], drop=FALSE] %*% coef[jj[[1]]] )
    eta1 <- drop( X[ , jj[[2]], drop=FALSE] %*% coef[jj[[2]]] )
    eta2 <- drop( X[ , jj[[3]], drop=FALSE] %*% coef[jj[[3]]] )

    mu      <- family$linfo[[1]]$linkinv( eta )
    sigma_v <- family$linfo[[2]]$linkinv( eta1 )
    sigma_u <- family$linfo[[3]]$linkinv( eta2 )
    s       <- family$s
    ##Define parameters
    l0<-dnormhnorm(x=y, mu=mu, sigma_v=sigma_v, sigma_u=sigma_u, s=s, deriv=4, xg=family$tri, log.p = T)

    l<-sum(l0)

    if (deriv>0) {
      ## the first derivatives

      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2))

      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(sigma_v),
                  family$linfo[[3]]$d2link(sigma_u))

    }

    l1<-attr(l0,"gradient")
    l2<-attr(l0,"hessian")

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {

      g3 <- cbind(family$linfo[[1]]$d3link(mu),
                  family$linfo[[2]]$d3link(sigma_v),
                  family$linfo[[3]]$d3link(sigma_u))

      l3<-attr(l0,"l3")

    }

    if (deriv>3) {

      g4 <- cbind(family$linfo[[1]]$d4link(mu),
                  family$linfo[[2]]$d4link(sigma_v),
                  family$linfo[[3]]$d4link(sigma_u))

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
      x1 <- x[ , jj[[1]], drop=FALSE]
      e1 <- E[ , jj[[1]], drop=FALSE] ## square root of total penalty

      #QR decomposition of covariates
      qr.x1 <- qr(x1)
      qr.x_all<-qr(x[,-which(duplicated(colnames(x))), drop=FALSE])

      s<-family$s
      para<-reparametrize(mean = qr.fitted(qr.x_all, y),
                            sd = abs(qr.fitted(qr.x_all, abs(qr.resid(qr.x_all, y)))),#sqrt(sum(qr.resid(qr.x_all, y)^2)/nrow(x1)),
                          skew = qr.fitted(qr.x_all, qr.resid(qr.x_all, y)^3), dist="normhnorm", s=s)

      yt1<-para$mu

      # 1) Ridge regression for the location parameter
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startMu <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
        startMu[ !is.finite(startMu) ] <- 0
      } else { startMu <- pen.reg(x1, e1, yt1) }
      start[jj[[1]]] <- startMu

      # 2) Ridge regression using log absolute residuals
      x1 <- x[,jj[[2]],drop=FALSE]
      e1 <- E[,jj[[2]],drop=FALSE]

      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startSigma_v <- qr.coef(qr(x1),c(log(para$sigma_v),rep(0,nrow(E))))
        startSigma_v[!is.finite(startSigma_v)] <- 0
      } else { startSigma_v <- pen.reg(x1,e1,log(para$sigma_v)) }
      start[jj[[2]]] <- startSigma_v


      # 3) Skewness
      x1 <-  x[, jj[[3]],drop=FALSE]
      e1 <-  E[, jj[[3]],drop=FALSE]

      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startSigma_u <- qr.coef(qr(x1), c(log(para$sigma_u),rep(0,nrow(E))))
        startSigma_u[!is.finite(startSigma_u)] <- 0
      } else { startSigma_u <- pen.reg(x1,e1,log(para$sigma_u)) }
      start[jj[[3]]] <- startSigma_u
    }
  }) ## initialize

  #Random number generation for normhnorm
  rd <- function(mu, wt, scale, s) {
    ## random number generation
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }

    return(rnormhnorm(n=nrow(mu), mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=s))
  } ## random number generation

  #Quantile function of normhnorm
  qf <- function(p, mu, wt, scale, s) {
    ##quantile function
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }

    return(qnormhnorm(p=p, mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=s))
  } ##quantile function

  #Cumulative distribution function of normhnorm
  cdf <- function(q, mu, wt, scale, logp, s) {
    ##cumulative distribution function
    mu <- as.matrix(mu)
    if(ncol(mu)==1){ mu <- t(mu) }

    return(pnormhnorm(q=q, mu = mu[,1], sigma_v = mu[,2], sigma_u = mu[,3], s=s, log.p = logp))
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
      for (i in 1:3) if (is.null(off[[i]])) off[[i]] <- 0
      lpi <- attr(X,"lpi")
      X1 <- X[,lpi[[1]],drop=FALSE]
      X2 <- X[,lpi[[2]],drop=FALSE]
      X3 <- X[,lpi[[3]],drop=FALSE]


      mu <- drop(X1%*%beta[lpi[[1]]] + off[[1]]) ## linear predictor for mu
      sigma_v <- drop(X2%*%beta[lpi[[2]]] + off[[2]])  ## linear predictor for sigma_v parameter
      sigma_u <- drop(X3%*%beta[lpi[[3]]] + off[[3]]) ## linear predictor for sigma_u parameter

      if (se) {
        stop("se still available for this family")
      }
    } else {
      se <- FALSE
      mu <- eta[,1]
      sigma_v <- eta[,2]
      sigma_u <- eta[,3]
    }

    para<-reparametrize(mu = mu, sigma_v = sigma_v, sigma_u = sigma_u, s=family$s)

    fv <- list(para$mean)
    if (!se) return(fv) else {
      stop("se not still available for this family")
      fv <- list(fit=para$mean, se.fit=para$sd)
      return(fv)
    }
  } ##predict

  structure(list(family="normhnorm",ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 s = s,
                 residuals=residuals,
                 rd=rd,
                 qf=qf,
                 cdf=cdf,
                 predict=predict,
                 linfo = stats, ## link information list
                 d2link=1, d3link=1, d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## normhnorm

