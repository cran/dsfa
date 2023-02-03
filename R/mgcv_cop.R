#' cop
#'
#' The cop implements multiple copula distributions in which the parameter \eqn{\delta} can depend on additive predictors.
#' Useable with \code{mgcv::gam}, the additive predictors are specified via a formula.
#'
#' @return An object inheriting from class \code{general.family} of the mgcv package, which can be used in the 'mgcv' and 'dsfa' package.
#'
#' @details Mostly internal function. Used with gam to fit copula model, which in turn is used for starting values. The function \code{gam} is from the mgcv package and is called with a formula.
#' The formula specifies a dummy on the left hand side and the structure of the additive predictor for the \eqn{\delta} parameter on the right hand side.
#' Link function is "generalized logit", where for each \code{distr_cop} argument there are specific \code{min} and \code{max} arguments, which are the boundaries of the parameter space.
#' Although the parameter space is larger in theory for some copulas, numeric under- and overflow limits the parameter space. The intervals for the parameter \code{delta} are provided by [delta_bounds()].
#' WARNING: Only the estimates of the coefficients are useful. The rest of the 'mgcv' object has no meaningful values,
#' as \code{\link[mgcv:gam]{gam()}} was more or less abused here.
#' 
#' @inheritParams dcop
#' @param link formula, specifying the link for \eqn{\delta} parameter. See details.
#' 
#' @examples 
#' #Set seed, sample size and type of copula
#' set.seed(1337)
#' N=500 #Sample size
#' 
#' #Generate covariates
#' x1<-runif(N,-1,1); x2<-runif(N,-1,1)
#' 
#' #Set parameters of the copula
#' eta<-matrix(1+2.5*x1+1.75*sin(pi*x2),nrow=N)
#' delta<-transform(x=eta, type="glogitinv", par=delta_bounds("normal"), deriv_order = 0)
#' 
#' #Simulate pseudo observations W and create dataset
#' dat<-as.data.frame(rcop(n=N, delta=delta, distr_cop="normal"))
#' dat$y<-1 #Add dummy response variable
#' 
#' #Write formulae for parameters
#' delta_formula<-y~x1+s(x2,bs="ps")
#' 
#' #Fit model
#' model<-mgcv::gam(delta_formula, data=dat, 
#'                  family=cop(W=dat[,1:2],
#'                             distr_cop="normal"), optimizer="efs")
#' 
#' #Smooth effects
#' #Effect of x2 on the predictor of delta
#' plot(model, select=1) #Estimated function
#' lines(x2[order(x2)], 1.75*sin(pi*x2[order(x2)])-
#'         mean(1.75*sin(pi*x2)), col=2) #True effect
#' 
#' @family copula
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
cop <- function(link = list("glogit"), W, distr_cop="normal"){
  #Object for mgcv::gam such that the copula can be estimated.
  
  #Number of parameters
  npar <- 1
  
  #Get boundaries of parameter space
  a<-delta_bounds(distr_cop)[1]
  b<-delta_bounds(distr_cop)[2]
  
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
  } else {
    stop(link[[1]], " link not available cop distribution")
  }
  
  #Calculate residuals
  residuals <- function(object, type = c("deviance", "response")) {
    
    type <- match.arg(type)
    
    out <- object$family$W
    
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
    jj <- list(1:ncol(X))#attr(X, "lpi") ## extract linear predictor index
    attr(jj,"overlap")<-FALSE
    
    #Number of parameters and observations
    npar <- 1
    n <- nrow(X)
    
    #Get additive predictors
    eta <-  drop( X %*% coef)
    
    #Additive predictors 2 parameters
    delta      <- matrix(family$linfo[[1]]$linkinv( eta ), nrow=n)
    
    #Evaluate ll
    l0<-dcop_cpp(u=family$W[,1], v=family$W[,2], p=delta, distr_cop=family$distr_cop, deriv_order=2, tri=family$tri_mat, logp = TRUE)
    
    #Assign sum of individual loglikehood contributions to l
    l<-sum(l0)
    
    if (deriv>0) {
      #First derivatives
      ig1 <- matrix(family$linfo[[1]]$mu.eta(eta),ncol=1)
      
      g2 <- matrix(family$linfo[[1]]$d2link(delta),ncol=1)
    }
    
    #Assign first and second derivative of ll to l1 and l2
    l1<-attr(l0,"d1")[,3, drop=FALSE]
    l2<-attr(l0,"d2")[,6, drop=FALSE]
    
    #Set default values 
    l3 <- l4 <- g3 <- g4 <- 0
    
    if (deriv>1) {
      g3 <- matrix(family$linfo[[1]]$d3link(delta),ncol=1)
      
      #Assign third derivative of ll to l3
      l3<-attr(l0,"d3")[,10, drop=FALSE]
    }
    
    if (deriv>3) {
      g4 <- matrix(family$linfo[[1]]$d4link(delta),ncol=1)
      
      #Assign fourth derivative of ll to l4
      l4<-attr(l0,"d4")[,15, drop=FALSE]
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
      
      if(family$distr_cop=="independent"){
        cop_object<-copula::indepCopula(dim = 2)
      }
      
      if(family$distr_cop=="normal"){
        # cop_dim<-sum(lower.tri(diag(ncol(family$W)), diag=T))#(copula::p2P(c(delta)))
        cop_object<-copula::normalCopula(dim = 2)
      }
      
      if(family$distr_cop=="clayton"){
        cop_object<-copula::claytonCopula(dim = 2)
      }
      
      if(family$distr_cop=="gumbel"){
        cop_object<-copula::gumbelCopula(dim = 2)
      }
      
      if(family$distr_cop=="frank"){
        cop_object<-copula::frankCopula(dim = 2)
      }
      
      if(family$distr_cop=="joe"){
        cop_object<-copula::joeCopula(dim = 2)
      }
      
      
      #Wrapper for Copula package function
      delta_start<-suppressWarnings(copula::fitCopula(copula=cop_object, data=as.matrix(family$W))@estimate)
      start <- c(delta_start,rep(0,ncol(x)-1))
    }
  }) 
  
  preinitialize <- function(G) {
    lpi<-list(1:ncol(G$X))#attr(X, "lpi") ## extract linear predictor index
    attr(lpi,"overlap")<-FALSE
    attr(G$X,"lpi")<-lpi
    
    return(list(X=G$X))
  }
  
  structure(list(family="cop", ll=ll, link=paste(link), nlp=npar,
                 tri = mgcv::trind.generator(npar), # symmetric indices for accessing derivative arrays
                 tri_mat = trind_generator(npar), # symmetric indices for accessing derivative matrices
                 initialize=initialize, #initial parameters
                 distr_cop = distr_cop, #specifiying copula distribution
                 W=W, #pseudo observations
                 preinitialize=preinitialize, # specify model matrix
                 residuals=residuals,
                 linfo = stats, # link information list
                 d2link=1, d3link=1, d4link=1, # signals to fix.family.link that all done
                 ls=1, # signals that ls not needed here
                 available.derivs = 2), # can use full Newton here
            class = c("general.family","extended.family","family"))
}

