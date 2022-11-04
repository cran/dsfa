#' reparametrize
#'
#' Transforms the given inputs to the parameters and the first three moments of the corresponding distribution. For the normal-halfnormal distribution the parametrization of the classical stochastic frontier as well as the skew-normal and centred skew-normal specification ar provided. For the normal-exponential an the specification via \eqn{\lambda} is available.
#'
#' @param mu vector of \eqn{\mu}
#' @param sigma_v vector of \eqn{\sigma_V}. Must be positive.
#' @param sigma_u vector of \eqn{\sigma_U}. Must be positive.
#' @param s \eqn{s=-1} for production and \eqn{s=1} for cost function.
#' @param lambda vector of \eqn{\lambda}. Must be positive.
#' @param mean vector of mean of \eqn{\mathcal{E}}
#' @param sd vector of standard deviation of \eqn{\mathcal{E}}. Must be positive.
#' @param skew  vector of skewness of \eqn{\mathcal{E}}.
#' @param par_u vector of \eqn{\sigma_U} or \eqn{\lambda}. Must be positive.
#' @param family \code{normhnorm} for normal-halfnormal and \code{normexp} for normal-exponential distribution.
#' @return Returns a data.frame with the parameter values for all specification.
#'
#' @details The following input combinations are allowed for the normal-halfnormal distribution
#' \itemize{
#'   \item \code{mu}, \code{sigma_v}, \code{sigma_u}, \code{s}
#'   \item \code{mean}, \code{sd}, \code{skew}, \code{family="normhnorm"}  with optional \code{s}  \eqn{\qquad,}
#' } while for the normal-exponential distribution the feasible inputs are
#' \itemize{
#'   \item \code{mu}, \code{sigma_v}, \code{lambda}, \code{s}
#'   \item \code{mean}, \code{sd}, \code{skew}, \code{family="normexp"} with optional \code{s} \eqn{\qquad.}
#' } Other input combinations are not feasible.
#'
#' @examples
#' #Normal-halfnormal distribution
#' para<-reparametrize(mu=1, sigma_v=2, sigma_u=3,s=-1)
#' reparametrize(mean=para$mean, sd=para$sd, skew=para$skew, family="normhnorm")
#'

#' #Normal-exponential distribution
#' para<-reparametrize(mu=1, sigma_v=2, lambda=1/3,s=-1)
#' reparametrize(mean=para$mean, sd=para$sd, skew=para$skew, family="normexp")
#'
#' @references
#' \itemize{
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{azzalini2013skew}{dsfa}
#' }
#' @export

reparametrize<-function(mu=NULL, sigma_v=NULL, sigma_u=NULL, s=NULL,
                        lambda=NULL, par_u=NULL,
                        mean=NULL, sd=NULL, skew=NULL, family=NULL){
  out<-NULL

  if(is.null(family)){
    family<-0
  }

  if(family=="normhnorm"&!is.null(par_u)){
     sigma_u<-par_u
  }

  if(family=="normhnorm"|!is.null(sigma_u)){
    #Function which transforms parameters
    #Inputs a vector of either
    #1) mean, sd, skew (csn parametrization)
    #2) mu, sigma, lambda (sf stabilized parametrization) with tau=1/sigma
    #3) mu, sigma_v, sigma_u (sf parametrization)
    #4) xi, omega, alpha (sn parametrization)
    #Return dataframe with all parameters

    #Idea is to get the parameters of sn distribution
    #and then calculate the rest from this parametrization

    #Numerical constant
    b<-sqrt(2/pi)

    #Parameters of csn 2 sn
    if(is.numeric(mean) & is.numeric(sd) & is.numeric(skew)){
      if(any(sd <= 0)){
        stop(paste("sd must be positive", "\n", ""))
      }
      if(is.null(s)){
        s<-sign(sum(skew))
      }
      max.skew<-0.5 * (4 - pi) * (2/(pi - 2))^1.5 - (.Machine$double.eps)^(1/4)
      if (any(abs(skew) > max.skew)){
        skew[abs(skew)>max.skew]<-s * 0.9 * max.skew
        # skew <- sign(skew) * 0.9 * max.skew
      }
      R<-(2*abs(skew)/(4-pi))^(1/3)
      alpha<-s*R/sqrt(2/pi-(1-2/pi)*R^2)
      delta<-alpha/sqrt(1+alpha^2)
      omega<-sd/sqrt(1-b^2*delta^2)
      xi<-mean-b*omega*delta
      tau<-1/omega
    }

    #sf 2 sn
    if(is.numeric(mu) & is.numeric(sigma_v) & is.numeric(sigma_u) & is.numeric(s)){
      #Check for correct inputs
      if (any(sigma_v <= 0)){
        stop(paste("sigma_v must be positive", "\n", ""))
      }

      if (any(sigma_u <= 0)){
        stop(paste("sigma_u must be positive", "\n", ""))
      }

      if (any(!s%in%c(-1,1))){
        stop(paste("s must be {-1, 1}", "\n", ""))
      }

      omega<-sqrt(sigma_v^2+sigma_u^2)
      alpha<-s*sigma_u/sigma_v
      xi<-mu
    }

    #Check if inputs are correct
    if(!(is.numeric(xi) & is.numeric(omega) & is.numeric(alpha))){
      stop("Incorrect inputs for the normhnorm distribution")
    }

    #Parameters of sn 2 sigma, lambda, mu, sigma_v, sigma_u, mean, sd, skew
    if(is.null(s)){
      s<-sign(sum(alpha))
    }
    # sigma<-omega
    tau<-1/omega
    # lambda<-s*alpha
    mu<-xi
    sigma_u<-omega*s*alpha/sqrt((s*alpha)^2+1)
    sigma_v<-sigma_u/(s*alpha)
    delta<-alpha/sqrt(1+alpha^2)
    mu_z<-b*delta
    sigma_z<-sqrt(1-mu_z^2)
    mean<-xi+omega*mu_z
    sd<-sqrt(omega^2*(1-mu_z^2))
    skew<-1/2*(4-pi)*mu_z^3/(1-mu_z^2)^(3/2)

    #Return output dataframe
    out<-data.frame(mu=mu, sigma_v=sigma_v, sigma_u=sigma_u, s=s,
                    xi=xi, omega=omega, alpha=alpha, delta=delta, par_u=sigma_u,
                    mean=mean, sd=sd, skew=skew)
  }

  if(family=="normexp"&!is.null(par_u)){
    lambda<-par_u
  }

  if(family=="normexp"|!is.null(lambda)){

    if(is.numeric(mean) & is.numeric(sd) & is.numeric(skew)){
      if(is.null(s)){
        s<-sign(sum(sign(skew)))
      }

      max.skew<-ifelse(s==1, 2-1e-10, 0-1e-10)
      min.skew<-ifelse(s==1, 0+1e-10, -2+1e-10)

      if (any(skew > max.skew)){
        skew[skew>max.skew]<-max.skew
      }
      if (any(skew < min.skew)){
        skew[skew<min.skew]<-min.skew
      }

      mu<-mean-s*sd*(s*skew/2)^(1/3)
      sigma_v<-sqrt(sd^2*(1-(s*skew/2)^(2/3)))
      lambda<-1/(sd*(s*skew/2)^(1/3))
      nu<-1/lambda
      # s<-sign(sum(mu))
    } else {
      if(is.numeric(mu) & is.numeric(sigma_v) & is.numeric(lambda) & is.numeric(s)){
          #mu<-s*mu
          mean<-mu+s*1/lambda
          sd<-sqrt(sigma_v^2+1/lambda^2)
          skew<-2/(sigma_v^3*s^3*lambda^3)*(1+1/(sigma_v^2*lambda^2))^(-3/2)
      }

      #Check if inputs are correct
      if(!(is.numeric(mu) & is.numeric(sigma_v) & is.numeric(lambda))){
        stop("Incorrect inputs for the normexp distribution")
      }
    }

    out<-data.frame(mu=mu, sigma_v=sigma_v, lambda=lambda, s=s,
                    par_u=lambda,
                    mean=mean, sd=sd, skew=skew)
  }

  if(is.null(out)){
    stop("Incorrect inputs for reparametrize")
  }

  #Return output
  return(out)
}

