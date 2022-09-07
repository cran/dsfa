#' cop family
#'
#' The cop family implements multiple copula distributions in which the \eqn{\tau} can depend on additive predictors. Useable only with \code{mgcv::gam}, the additive predictors is specified via a formulae.
#'
#'
#' @param link specifying the link for the \eqn{\tau}.
#' @param u1 probability integral transformed observations for margin 1.
#' @param u2 probability integral transformed observations for margin 2.
#' @param cop_family integer, defines the copula family:\cr
#' `1` = Gaussian copula \cr
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the dsfa package.
#'
#' @details Internal function. Used with gam to fit copula model, which in turn is used for starting values. The function \code{gam} is from the mgcv package is called with a formulae.
#' @export
#' @keywords internal
cop <- function(link = list("tanh"), u1, u2, cop_family=1){
  #Object for mgcv::gam such that the copula can be estimated.

  #Number of parameters
  npar <- 1

  #Link functions
  if (length(link) != npar) stop("cop requires 1 links specified as character strings")
  okLinks <- list("tanh")
  stats <- list()
  param.names <- c("Tau")

  if (link[[1]] %in% okLinks[[1]]) {
    stats[[1]] <- list()
    stats[[1]]$valideta <- function(eta) TRUE
    stats[[1]]$link = link[[1]]
    # y<-log((x+1)/(1-x))
    # x<-(exp(y)-1)/(exp(y)+1)

    stats[[1]]$linkfun <- eval(parse(text = paste("function(mu) log((mu+1)/(1-mu))")))
    stats[[1]]$linkinv <- eval(parse(text = paste("function(eta) (exp(eta)-1)/(exp(eta)+1)")))
    stats[[1]]$mu.eta <- eval(parse(text = paste("function(eta) 1/(1 + cosh(eta))")))
    stats[[1]]$d2link <- eval(parse(text = paste("function(mu) 4*mu/(mu^2-1)^2"))) #-((2 * exp(eta) * (-1 + exp(eta)))/(1 + exp(eta))^3)
    stats[[1]]$d3link <- eval(parse(text = paste("function(mu) -4*(3*mu^2+1)/(mu^2-1)^3"))) #(2 * exp(eta) * (1 - 4 * exp(eta) + exp(2 * eta)))/(1 + exp(eta))^4
    stats[[1]]$d4link <- eval(parse(text = paste("function(mu) 48*(mu^3+mu)/(mu^2-1)^4"))) #-((2 * exp(eta)*  (-1 + 11 * exp(eta) - 11 * exp(2 * eta) + exp(3 * eta)))/(1 + exp(eta))^5)
  }
  else stop(link[[1]], " link not available joint distribution")

  #Calculate residuals
  residuals <- function(object, type = c("deviance", "response")) {

    type <- match.arg(type)

    rsd <- cbind(object$family$u1,object$family$u2)

    if (type=="response") return(rsd) else
      return(rsd)
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

    npar <- 1
    n <- nrow(X)

    eta <-  drop( X[ , jj[[1]], drop=FALSE] %*% coef[jj[[1]]] )


    Tau      <- family$linfo[[1]]$linkinv( eta )

    ##Define parameters
    l0<-dcop(U=cbind(family$u1,family$u2), Tau=Tau, family=family$cop_family, deriv=2, disjoint=T, log.p = TRUE)

    l<-sum(l0)

    if (deriv>0) {
      ## the first derivatives

      ig1 <- matrix(family$linfo[[1]]$mu.eta(eta), ncol=1)

      g2 <- matrix(family$linfo[[1]]$d2link(Tau), ncol=1)
    }

    l1<-attr(l0,"gradient")[,3,drop=FALSE]
    l2<-attr(l0,"hessian")[,6,drop=FALSE]

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      g3 <- matrix(family$linfo[[1]]$d3link(Tau), ncol=1)

      l3<-attr(l0,"l3")

    }

    if (deriv>3) {
      g4 <- matrix(family$linfo[[1]]$d4link(Tau), ncol=1)

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


      # 7) Copula
      x1 <-  x[, jj[[1]],drop=FALSE]
      e1 <-  E[, jj[[1]],drop=FALSE]
      #
      # #Evaluate copula at probability integral transformed observations
      # tau_start<-cor(x=pcomperr(q=y2, mu=para1$mu, sigma_v=para1$sigma_v, par_u=para1$par_u, s=s1, dist=family$dist[1], log.p=TRUE),
      #                y=pcomperr(q=y1, mu=para2$mu, sigma_v=para2$sigma_v, par_u=para2$par_u, s=s2, dist=family$dist[2], log.p=TRUE),
      #                use = "pairwise.complete.obs")
      tau_start<-cor(family$u1,family$u2)
      # family$linfo[[7]]$linkinv(par_start)
      if (use.unscaled) {
        x1 <- rbind(x1,e1)
        startTau <- qr.coef(qr(x1), c(rep(family$linfo[[1]]$linkinv(tau_start),length(family$u1)),rep(0,nrow(E))))
        startTau[!is.finite(startTau)] <- 0
      } else { startTau <- pen.reg(x1,e1,rep(family$linfo[[1]]$linkinv(tau_start),length(family$u1))) }
      start[jj[[1]]] <- startTau
    }
    # print(start)
  }) ## initialize


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

      if (se) {
        stop("se still available for this family")
      }
    } else {
      se <- FALSE
    }

    fv <- list(cbind(family$u1,family$u2))
    if (!se) return(fv) else {
      stop("se not still available for this family")
      fv <- list(fit=cbind(family$u1,family$u2), se.fit=cbind(family$u1,family$u2))
      return(fv)
    }
  } ##predict


  structure(list(family="cop", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 u1=u1,
                 u2=u2,
                 cop_family=cop_family,
                 residuals=residuals,
                 predict=predict,
                 linfo = stats, ## link information list
                 d2link=1, d3link=1, d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 2 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## cop
