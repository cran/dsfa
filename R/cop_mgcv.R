#' cop family
#'
#' The cop family implements multiple copula distributions in which the \eqn{\delta} can depend on additive predictors. Useable only with \code{\link[mgcv:gam]{gam()}}, the additive predictors is specified via a formulae.
#'
#' @param link specifying the link for the \eqn{\delta}.
#' @param W matrix of pseudo observations. Must have at least two columns.
#' @param family_cop string, defines the copula family:\cr
#' `independent` = Independence copula \cr
#' `normal` = Gaussian copula \cr
#' `clayton` = Clayton copula \cr
#' `gumbel` = Gumbel copula \cr
#' `frank` = Frank copula \cr
#' `joe` = Joe copula \cr
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the dsfa package.
#'
#' @details Internal function. Used with gam to fit copula model, which in turn is used for starting values. The function \code{gam} is from the mgcv package is called with a formula.
#' The formula specifies the a dummy on the left hand side and the structure of the additive predictor for \eqn{\delta} parameter on the right hand side. Link function is "generalized logit", where for each family_cop argument there are specific \code{a} and \code{b} arguments, which are the boundaries of the parameter space. Although the parameter space is
#' larger in theory, numeric under- and overflow limit the parameter space.
#' \enumerate{
#'  \item  `independent`, a=0 and b=1
#' \item  `normal`, a=-1 and b=1
#' \item  `clayton`, a=1e-16 and b=28
#' \item  `gumbel`, a=1 and b=17
#' \item  `frank`, a=-35 and b=35
#' \item  `joe`, a=1e-16 and b=30
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1337)
#' N=1000 #Sample size
#' x1<-runif(N,-1,1)
#' eta<-1+2.5*x1
#' delta<-rsp(x=eta, link="glogit", a=-1, b=1)
#' dat<-as.data.frame(rcop(n=N, delta=delta, family="normal"))
#' dat$y<-1
#' model<-mgcv::gam(y~x1, data=dat,
#'                  family=cop(W=as.matrix(dat[,1:2, drop=FALSE]), family_cop="normal"),
#'                  optimizer="efs")
#'}
#'
#' @export
#' @keywords internal
cop <- function(link = list("glogit"), W, family_cop="normal"){
  #Object for mgcv::gam such that the copula can be estimated.

  #Number of parameters
  npar <- 1

  #Delta parameter space boundaries
  if(family_cop=="independent"){
    a<-0
    b<-1
  }

  if(family_cop=="normal"){
    a<--1
    b<-1
  }

  if(family_cop=="clayton"){
    a<-0+1e-16
    b<-28
  }

  if(family_cop=="gumbel"){
    a<-1
    b<-17
  }

  if(family_cop=="frank"){
    a<--35
    b<-35

    stop(paste("The pdf of the frank copula is not functional yet.", "\n", ""))

  }

  if(family_cop=="joe"){
    a<-1+1e-6
    b<-30
  }

  #Link functions
  if (length(link) != npar) stop("cop requires 1 links specified as character strings")
  okLinks <- list("glogit")
  stats <- list()
  param.names <- c("delta")

  if (link[[1]] %in% okLinks[[1]]) {
    stats[[1]] <- list()
    stats[[1]]$valideta <- function(eta) TRUE
    stats[[1]]$link = link[[1]]

    stats[[1]]$linkfun <- eval(parse(text = paste("function(mu) log((-",a,"+mu)/(",b,"-mu))")))#eval(parse(text = paste("function(mu) log((mu+1)/(1-mu))")))
    stats[[1]]$linkinv <- eval(parse(text = paste("function(eta) exp(eta)/(1+exp(eta))*(",b,"-",a,")+",a)))
    stats[[1]]$mu.eta <- eval(parse(text = paste("function(eta) exp(eta)*(",b,"-",a,")/(exp(eta)+1)^2")))
    stats[[1]]$d2link <- eval(parse(text = paste("function(mu) 1/(",b,"-mu)^2-1/(",a,"-mu)^2")))#eval(parse(text = paste("function(mu) 4*mu/(mu^2-1)^2"))) #-((2 * exp(eta) * (-1 + exp(eta)))/(1 + exp(eta))^3)
    stats[[1]]$d3link <- eval(parse(text = paste("function(mu) 2*(1/(",b,"-mu)^3+1/(-",a,"+mu)^3)"))) #(2 * exp(eta) * (1 - 4 * exp(eta) + exp(2 * eta)))/(1 + exp(eta))^4
    stats[[1]]$d4link <- eval(parse(text = paste("function(mu) 6/(",b,"-mu)^4-6/(",a,"-mu)^4"))) #-((2 * exp(eta)*  (-1 + 11 * exp(eta) - 11 * exp(2 * eta) + exp(3 * eta)))/(1 + exp(eta))^5)
  }
  else stop(link[[1]], " link not available cop distribution")

  #Calculate residuals
  residuals <- function(object, type = c("deviance", "response")) {

    type <- match.arg(type)

    rsd <- object$family$W

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

    jj <- list(1:ncol(X))#attr(X, "lpi") ## extract linear predictor index
    attr(jj,"overlap")<-FALSE

    npar <- 1
    n <- nrow(X)
    eta <-  drop( X %*% coef)
    delta      <- matrix(family$linfo[[1]]$linkinv( eta ), nrow=n)

    ##Define parameters
    l0<-dcop(W=as.matrix(family$W, nrow=n), delta=delta, family_cop=family$family_cop, deriv=2, log.p = TRUE)

    l<-sum(l0)

    if (deriv>0) {
      ## the first derivatives

      ig1 <- matrix(family$linfo[[1]]$mu.eta(eta), ncol=1)

      g2 <- matrix(family$linfo[[1]]$d2link(delta), ncol=1)
    }

    l1<-attr(l0,"gradient")[,3,drop=FALSE]
    l2<-attr(l0,"hessian")[,6,drop=FALSE]

    l3 <- l4 <- g3 <- g4 <- 0 ## defaults

    if (deriv>1) {
      g3 <- matrix(family$linfo[[1]]$d3link(delta), ncol=1)

      l3<-attr(l0,"l3")

    }

    if (deriv>3) {
      g4 <- matrix(family$linfo[[1]]$d4link(delta), ncol=1)

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
    ## idea is to get starting values utilizing the method of moments
    n <- rep(1, nobs)
    ## should E be used unscaled or not?..
    use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {

      if(family$family_cop=="independent"){
        cop_object<-copula::indepCopula(dim = 2)
      }

      if(family$family_cop=="normal"){
        # cop_dim<-sum(lower.tri(diag(ncol(family$W)), diag=T))#(copula::p2P(c(delta)))
        cop_object<-copula::normalCopula(dim = ncol(family$W))
      }

      if(family$family_cop=="clayton"){
        cop_object<-copula::claytonCopula(dim = 2)
      }

      if(family$family_cop=="gumbel"){
        cop_object<-copula::gumbelCopula(dim = 2)
      }

      if(family$family_cop=="frank"){
        cop_object<-copula::frankCopula(dim = 2)
      }

      if(family$family_cop=="joe"){
        cop_object<-copula::joeCopula(dim = 2)
      }


      #Wrapper for Copula package function
      delta_start<-suppressWarnings(copula::fitCopula(copula=cop_object, data=as.matrix(family$W))@estimate)
      start <- c(delta_start,rep(0,ncol(x)-1))
    }
    #
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
      # lpi <- attr(X,"lpi")
      lpi <- list(1:ncol(X))#attr(X, "lpi") ## extract linear predictor index
      attr(lpi,"overlap")<-FALSE
      X1 <- X[,lpi[[1]],drop=FALSE]

      if (se) {
        stop("se still available for this family")
      }
    } else {
      se <- FALSE
    }

    fv <- list(family$W)
    if (!se) return(fv) else {
      stop("se not still available for this family")
      fv <- list(fit=family$W, se.fit=family$W)
      return(fv)
    }
  } ##predict


  preinitialize <- function(G) {
    lpi<-list(1:ncol(G$X))#attr(X, "lpi") ## extract linear predictor index
    attr(lpi,"overlap")<-FALSE
    attr(G$X,"lpi")<-lpi

    return(list(X=G$X))
  }

  structure(list(family="cop", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), ## symmetric indices for accessing derivative arrays
                 initialize=initialize,
                 W=W,
                 a=a,
                 b=b,
                 family_cop=family_cop,
                 residuals=residuals,
                 predict=predict,
                 preinitialize=preinitialize,
                 linfo = stats, ## link information list
                 d2link=1, d3link=1, d4link=1, ## signals to fix.family.link that all done
                 ls=1, ## signals that ls not needed here
                 available.derivs = 0 ## can use full Newton here
  ),class = c("general.family","extended.family","family"))
} ## cop
