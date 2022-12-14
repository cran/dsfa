#' Normal-exponential distribution
#'
#' Probablitiy density function, distribution, quantile function and random number generation for the normal-exponential distribution.
#'
#' @param x vector of quantiles.
#' @param mu vector of \eqn{\mu}
#' @param sigma_v vector of \eqn{\sigma_V}. Must be positive.
#' @param lambda vector of \eqn{\lambda}. Must be positive.
#' @param s \eqn{s=-1} for production and \eqn{s=1} for cost function.
#' @param deriv derivative of order \code{deriv} of the log density. Available are \code{0},\code{2} and \code{4}.
#' @param tri optional, index arrays for upper triangular matrices, generated by \code{\link[mgcv:trind.generator]{trind.generator()}} and supplied to \code{chainrule()}.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param check logical; if TRUE, check inputs.
#'
#' @return \code{dnormexp} gives the density, \code{pnormexp} give the distribution function, \code{qnormexp} gives the quantile function, and \code{rnormexp} generates random numbers, with given parameters.  If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian}, \code{l3} and \code{l4} of the output of the density.
#'
#' @details A random variable \eqn{\mathcal{E}} follows a normal-exponential distribution if \eqn{\mathcal{E} = V + s \cdot U }, where \eqn{V \sim N(\mu, \sigma_V^2)} and \eqn{U \sim Exp(\lambda)}.
#' The density is given by \deqn{f_\mathcal{E}(\epsilon)=\frac{\lambda}{2} \exp \{\lambda (s \mu) + \frac{1}{2} \lambda^2 \sigma_V^2-\lambda (s \epsilon) \} 2 \Phi(\frac{1}{\sigma_V} (-s \mu)-\lambda \sigma_V+\frac{1}{\sigma_V}(s \epsilon)) \qquad,}
#' where \eqn{s=-1} for production and \eqn{s=1} for cost function.
#'
#'
#' @examples
#' pdf <- dnormexp(x=seq(-3, 3, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
#' cdf <- pnormexp(q=seq(-3, 3, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
#' q <- qnormexp(p=seq(0.1, 0.9, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
#' r <- rnormexp(n=100, mu=1, sigma_v=2, lambda=1/3, s=1)
#'
#' @references
#' \itemize{
#' \item \insertRef{meeusen1977efficiency}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{schmidt2020analytic}{dsfa}
#' \item \insertRef{gradshteyn2014table}{dsfa}
#' }
#' @export
#'
#dnormexp
dnormexp <- function(x, mu=0, sigma_v=1, lambda=1, s=-1, deriv=0, tri=NULL, log.p = FALSE, check=TRUE){
  #Density function of the normhnorm distribution

  #Check inputs arguments
  if(check){
    arguments<-check_arguments(x=x, mu=mu, sigma_v=sigma_v, par_u=lambda, s=s, family="normexp", deriv=deriv)
    if(arguments$recycle){
      x<-arguments$x
      mu<-arguments$mu
      sigma_v<-arguments$sigma_v
      lambda<-arguments$par_u
    }
  }

  v<-sigma_v
  l<-lambda

  z<-1/v*(-s*mu)-l*v+1/v*s*x#1/v*mu-l*v+1/v*x #
  zeta_z<-zeta(z, deriv=deriv)

  #Get loglikelihood
  out  <-  log(l)-log(2)+l*s*mu+0.5*l^2*v^2-l*s*x+zeta_z

  # out <- dsn(x=x, xi=para$xi, omega=para$omega, alpha=para$alpha, log=log.p)

  #Exponentialize output
  if (!log.p){
    out <- exp(out)
  }

  if(deriv>0){
    #First derivatives of z wrt mu, v, l
    d1zdmu1<--s/v
    d1zdv1<--l + s*mu/v^2 - s*x/v^2
    d1zdl1<--v

    #Second derivatives of z wrt mu, v, l
    d2zdmu2<-0
    d2zdmu1dv1<-(s/v^2)
    d2zdmu1dl1<-0
    d2zdv2<--((2 * s*mu)/v^3) + (2 * s*x)/v^3
    d2zdv1dl1<--1
    d2zdl2<-0

    attr(z,"gradient")<-cbind(d1zdmu1, d1zdv1, d1zdl1)
    attr(z,"hessian")<-cbind(d2zdmu2, d2zdmu1dv1, d2zdmu1dl1, d2zdv2, d2zdv1dl1, d2zdl2)

    if(deriv>2){
      #Third derivatives of z wrt mu, v, l
      d3zdmu3<-0
      d3zdmu2dv1<-0
      d3zdmu2dl1<-0
      d3zdmu1dv2<--2*s/v^3
      d3zdmu1dv1dl1<-0
      d3zdmu1dl2<-0
      d3zdv3<-(6 * mu*s)/v^4 - (6 * s* x)/v^4
      d3zdv2dl1<-0
      d3zdv1dl2<-0
      d3zdl3<-0

      #Fourth derivatives of z wrt mu, v, l
      d4zdmu4<-0
      d4zdmu3dv1<-0
      d4zdmu3dl1<-0
      d4zdmu2dv2<-0
      d4zdmu2dv1dl1<-0
      d4zdmu2dl2<-0
      d4zdmu1dv3<-6*s/v^4
      d4zdmu1dv2dl1<-0
      d4zdmu1dv1dl2<-0
      d4zdmu1dl3   <-0
      d4zdv4<- -((24 * s*mu)/v^5) + (24 * s* x)/v^5
      d4zdv3dl1<-0
      d4zdv2dl2<-0
      d4zdv1dl3<-0
      d4zdl4<-0

      attr(z,"l3")<-cbind(d3zdmu3,d3zdmu2dv1,d3zdmu2dl1,d3zdmu1dv2,d3zdmu1dv1dl1,d3zdmu1dl2,d3zdv3,d3zdv2dl1,d3zdv1dl2,d3zdl3)
      attr(z,"l4")<-cbind(d4zdmu4,d4zdmu3dv1,d4zdmu3dl1,
                d4zdmu2dv2,d4zdmu2dv1dl1,d4zdmu2dl2,
                d4zdmu1dv3,d4zdmu1dv2dl1,d4zdmu1dv1dl2,d4zdmu1dl3,
                d4zdv4,d4zdv3dl1,d4zdv2dl2,d4zdv1dl3,
                d4zdl4)
    }

    l0<-chainrule(f=zeta_z, g=z, deriv=deriv, tri=tri)

    #Calculate l1,l2,l3,l4
    ## the first derivatives
    ## order mu, v, l
    l1 <- attr(l0,"gradient")

    l1[,1] <-l*s+l1[,1]
    l1[,2] <-l^2 * v+l1[,2]
    l1[,3] <-1/l + mu*s + l * v^2 - s* x + l1[,3]

    ## the second derivatives
    ## order mumu, muv, muu, vv, vu, uu
    l2 <- attr(l0,"hessian")

    l2[,1] <-0+l2[,1]
    l2[,2] <-0+l2[,2]
    l2[,3] <-s+l2[,3]
    l2[,4] <-l^2+l2[,4]
    l2[,5] <-2 * l * v +l2[,5]
    l2[,6] <--(1/l^2) + v^2 +l2[,6]

    #Set gradient and hessian as attributes
    attr(out,"gradient")<-l1
    attr(out,"hessian")<-l2

    if(deriv>2){
      # order mumumu, mumul, mumul,
      #       muvv, muvl, mull
      #       vvv, vvl, vll
      #       lll
      l3 <- attr(l0,"l3")

      l3[,1] <- 0+l3[,1]
      l3[,2] <- 0+l3[,2]
      l3[,3] <- 0+l3[,3]
      l3[,4] <- 0+l3[,4]
      l3[,5] <- 0+l3[,5]
      l3[,6] <- 0+l3[,6]
      l3[,7] <- 0+l3[,7]
      l3[,8] <- 2*l+l3[,8]
      l3[,9] <- 2*v+l3[,9]
      l3[,10]<- 2/l^3+l3[,10]

      #order mumumumu, mumumuv, mumumul,
      #      mumuvv, mumuvl, mumull
      #      muvvv, muvvl, muvll
      #      mulll, vvvv, vvvl,
      #      vvll, vlll, llll

      l4 <- attr(l0,"l4")

      l4[,1]<-l4[,1]
      l4[,2]<-l4[,2]
      l4[,3]<-l4[,3]
      l4[,4]<-l4[,4]
      l4[,5]<-l4[,5]
      l4[,6]<-l4[,6]
      l4[,7]<-l4[,7]
      l4[,8]<-l4[,8]
      l4[,9]<-l4[,9]
      l4[,10]<-l4[,10]
      l4[,11]<-l4[,11]
      l4[,12]<-l4[,12]
      l4[,13]<-2+l4[,13]
      l4[,14]<-l4[,14]
      l4[,15]<--(6/l^4)+l4[,15]

      #Set third and fourth derivatives as attributes
      attr(out,"l3")<-l3
      attr(out,"l4")<-l4
    }
  }

  #Return ouptut
  return(out)
}

#' @describeIn dnormexp distribution function for the normal-exponential distribution.
#' @param q vector of probabilities.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le x]}  otherwise, \eqn{P[X > x]}.
#' @export
pnormexp <- function(q, mu=0, sigma_v=1, lambda=1, s=-1, deriv=0, tri=NULL, lower.tail = TRUE, log.p = FALSE, check=TRUE){
  #Distribution function of the normexp distribution

  #Check inputs arguments
  if(check){
    arguments<-check_arguments(q=q, mu=mu, sigma_v=sigma_v, par_u=lambda, s=s, family="normexp", deriv=deriv)
    if(arguments$recycle){
      q<-arguments$x
      mu<-arguments$mu
      sigma_v<-arguments$sigma_v
      lambda<-arguments$par_u
    }
  }

  x<-q
  v<-sigma_v
  l<-lambda


  z<-s*(x-mu)/v
  out <- stats::pnorm(z)-1/l*dnormexp(x=x, mu=mu, sigma_v = v, lambda=l, s=s)

  if(deriv>0){
    n<-length(x)
    #Calculate l1,l2,l3,l4
    ## the first derivatives
    ## order mu, v, l
    l1 <- matrix(0, n, 3)

    l1[,1]<--(1/2)*exp(l * (mu * s + 0.5 * l * v^2 - s * x)) * l * s * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v))
    l1[,2]<-l * (0.398942 * exp(-((s^2 * (mu - x)^2)/(2 * v^2))) - 0.5 * exp(l * (mu * s + 0.5 * l * v^2 - s * x))*  l * v * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))
    l1[,3]<-1/2 * exp(l * (mu * s + 0.5 * l * v^2 - s * x)) * (exp(-((mu * s + l * v^2 - s * x)^2/(2 * v^2))) * sqrt(2/pi)* v + (-mu * s - l * v^2 + s * x) * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))

    ## the second derivatives
    ## order mumu, muv, mul, vv, vl, ll
    l2 <- matrix(0, n, 6)

    l2[,1]<-l * s^2 * ((0.398942 * exp(-((s^2 * (mu - x)^2)/(2 * v^2))))/v - 0.5 * exp(l *(mu * s + 0.5 * l * v^2 - s * x)) * l * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))
    l2[,2]<-(exp(-((mu^2 * s^2 + 0.5 * l^2 * v^4 + l *s *v^2 *x + s^2 *x^2 + mu *s *(l *v^2 - 2* s *x))/v^2)) *s *(exp(0.5 * l^2 * v^2 + (s^2 * (mu - x)^2)/(2 * v^2) +
                        l * s * (mu + x)) * (0.398942 * mu^2 * s^2 + 0.398942 * l^2 * v^4 + 0.398942 * s^2 * x^2 + mu * s * (-0.398942 * l * v^2 - 0.797885 * s * x) + v^2 * (-0.398942 + 0.398942 * l * s * x)) +
                        exp(2 *l *s* x + (mu * s + l * v^2 - s * x)^2/(2 * v^2)) * (0.398942 * v^2 + s^2 * (-0.398942 * mu^2 + 0.797885 * mu * x - 0.398942 * x^2)) -
                        0.5 * exp((mu^2 * s^2 + l^2 * v^4 + s^2 * x^2 + mu * s * (2 * l * v^2 - 2* s *x))/v^2)* l^3 * v^5 * erfc((mu * s + l * v^2 - s * x)/(sqrt(2)* v))))/v^4
    l2[,3]<-exp(l * (mu * s + 0.5 * l * v^2 -s * x)) * s * (0.398942 * exp(-((mu * s + l * v^2 - s * x)^2/(2 * v^2))) * l * v + (-0.5 + l * (-0.5 * mu * s - 0.5 * l * v^2 + 0.5 * s * x))*
            erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))
    l2[,4]<-(exp(-((s^2 * (mu - x)^2)/(2 * v^2))) * (1.11022*10^-16 * mu^3 * s^3 + 0.398942 * l^3 * v^6 + 0.398942 * l^2 * s * v^4 * x + 0.398942 * l * s^2 * v^2 * x^2 -
                  1.11022*10^-16 * s^3 * x^3 +  mu^2 * s^2 * (0.398942 * l * v^2 - 2.22045*10^-16 * s * x) +mu * s * (-0.398942 * l^2 * v^4 - 0.797885 * l * s * v^2 * x +
                  2.22045*10^-16 * s^2 * x^2) + exp(0.5 * l^2 * v^2 + l * s * (mu - x) + (s^2 * (mu - x)^2)/(2 * v^2)) * l^2 * v^5 * (-0.5 - 0.5 * l^2 * v^2) * erfc((mu * s + l * v^2 - s * x)/(
                  sqrt(2) *v))))/v^5
    l2[,5]<-exp((s^2 * (-0.5 * mu^2 + mu * x - 0.5 * x^2))/v^2) * (0.398942 + 0.398942 * l^2 * v^2 + exp((mu * s + l * v^2 - s * x)^2/(2 * v^2))*
                                                           l * v * (-1 + l * (-0.5 * mu * s - 0.5 * l * v^2 + 0.5 * s * x)) * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))
    l2[,6]<-1/2 * exp(l * (mu * s + 0.5 * l * v^2 - s * x)) * (exp(-((mu * s + l * v^2 - s * x)^2/(2 * v^2))) * (0.797885 * mu * s * v + 0.797885 * l * v^3 - 0.797885 * s * v *x) -
                            v^2 * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)) -(mu * s + l * v^2 - s * x)^2 * erfc((mu * s + l * v^2 - s * x)/(sqrt(2) * v)))
    #Set gradient and hessian as attributes
    attr(out,"gradient")<-l1
    attr(out,"hessian")<-l2
  }

  #Correct for numerical inaccuracy
  out[out>1-1e-14]<-1-1e-14
  out[out<0+1e-14]<-0+1e-14

  #Compute upper tail
  if(!lower.tail){
    out<-1-out
  }

  #Logarithmize output
  if(log.p){
    out<-log(out)
  }

  #Return output
  return(out)
}

#' @describeIn dnormexp quantile function for the normal-exponential distribution.
#' @param p vector of quantiles.
#' @export
qnormexp <- function(p, mu=0, sigma_v=1, lambda=1, s=-1, lower.tail = TRUE, log.p = FALSE, check=TRUE){
  #Quantile function of the normhnorm distribution

  #Check inputs arguments
  if(check){
    arguments<-check_arguments(p=p, mu=mu, sigma_v=sigma_v, par_u=lambda, s=s, family="normexp", deriv=0)
    if(arguments$recycle){
      p<-arguments$x
      mu<-arguments$mu
      sigma_v<-arguments$sigma_v
      lambda<-arguments$par_u
    }
  }

  if (s==-1){
    stop(paste("The quantile function of the normal-exponential distribution with s=-1 is not functional yet.", "\n", ""))
  }

  n<-length(p)

  #Transform parameters
  h1 <- function(q) {
    pnormexp(q, mu = mu[i], sigma_v = sigma_v[i], lambda = lambda[i], s=s) - p[i]
  }
  h <- function(q) {
    pnormexp(q, mu = mu[i], sigma_v = sigma_v[i], lambda = lambda[i], s=s)
  }

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  lp <- max(length(p), length(mu), length(sigma_v), length(lambda))
  p <- rep(p, length = lp)
  sigma_v <- rep(sigma_v, length = lp)
  mu <- rep(mu, length = lp)
  lambda <- rep(lambda, length = lp)
  q <- rep(0, lp)
  for (i in seq(along = p)) {
    if (h(mu[i]) < p[i]) {
      interval <- c(mu[i], mu[i] + sigma_v[i])
      j <- 2
      while (h(interval[2]) < p[i]) {
        interval[2] <- mu[i] + j * sigma_v[i]
        j <- j + 1
      }
    }
    else {
      interval <- c(mu[i] - sigma_v[i], mu[i])
      j <- 2
      while (h(interval[1]) > p[i]) {
        interval[1] <- mu[i] - j * sigma_v[i]
        j <- j + 1
      }
    }
    q[i] <- stats::uniroot(h1, interval)$root
  }
  out<-q

  #Return output
  return(out)
}

#' @describeIn dnormexp random number generation for the normal-exponential distribution.
#' @param n number of observations.
#' @export
rnormexp <- function(n, mu=0, sigma_v=1, lambda=1, s=-1, check=TRUE){
  #Function to generate n random numbers from the normexp distribution
  #Check inputs arguments
  if(check){
    arguments<-check_arguments(n=n, mu=mu, sigma_v=sigma_v, par_u=lambda, s=s, family="normexp", deriv=0)
    if(arguments$recycle){
      n<-arguments$x
      mu<-arguments$mu
      sigma_v<-arguments$sigma_v
      lambda<-arguments$par_u
    }
  }

  out<-stats::rnorm(n, mean=mu, sd=sigma_v)+s*stats::rexp(n, lambda)

  #Return output
  return(out)
}
