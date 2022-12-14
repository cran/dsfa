% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normexp.R
\name{dnormexp}
\alias{dnormexp}
\alias{pnormexp}
\alias{qnormexp}
\alias{rnormexp}
\title{Normal-exponential distribution}
\usage{
dnormexp(
  x,
  mu = 0,
  sigma_v = 1,
  lambda = 1,
  s = -1,
  deriv = 0,
  tri = NULL,
  log.p = FALSE,
  check = TRUE
)

pnormexp(
  q,
  mu = 0,
  sigma_v = 1,
  lambda = 1,
  s = -1,
  deriv = 0,
  tri = NULL,
  lower.tail = TRUE,
  log.p = FALSE,
  check = TRUE
)

qnormexp(
  p,
  mu = 0,
  sigma_v = 1,
  lambda = 1,
  s = -1,
  lower.tail = TRUE,
  log.p = FALSE,
  check = TRUE
)

rnormexp(n, mu = 0, sigma_v = 1, lambda = 1, s = -1, check = TRUE)
}
\arguments{
\item{x}{vector of quantiles.}

\item{mu}{vector of \eqn{\mu}}

\item{sigma_v}{vector of \eqn{\sigma_V}. Must be positive.}

\item{lambda}{vector of \eqn{\lambda}. Must be positive.}

\item{s}{\eqn{s=-1} for production and \eqn{s=1} for cost function.}

\item{deriv}{derivative of order \code{deriv} of the log density. Available are \code{0},\code{2} and \code{4}.}

\item{tri}{optional, index arrays for upper triangular matrices, generated by \code{\link[mgcv:trind.generator]{trind.generator()}} and supplied to \code{chainrule()}.}

\item{log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{check}{logical; if TRUE, check inputs.}

\item{q}{vector of probabilities.}

\item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]}  otherwise, \eqn{P[X > x]}.}

\item{p}{vector of quantiles.}

\item{n}{number of observations.}
}
\value{
\code{dnormexp} gives the density, \code{pnormexp} give the distribution function, \code{qnormexp} gives the quantile function, and \code{rnormexp} generates random numbers, with given parameters.  If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian}, \code{l3} and \code{l4} of the output of the density.
}
\description{
Probablitiy density function, distribution, quantile function and random number generation for the normal-exponential distribution.
}
\details{
A random variable \eqn{\mathcal{E}} follows a normal-exponential distribution if \eqn{\mathcal{E} = V + s \cdot U }, where \eqn{V \sim N(\mu, \sigma_V^2)} and \eqn{U \sim Exp(\lambda)}.
The density is given by \deqn{f_\mathcal{E}(\epsilon)=\frac{\lambda}{2} \exp \{\lambda (s \mu) + \frac{1}{2} \lambda^2 \sigma_V^2-\lambda (s \epsilon) \} 2 \Phi(\frac{1}{\sigma_V} (-s \mu)-\lambda \sigma_V+\frac{1}{\sigma_V}(s \epsilon)) \qquad,}
where \eqn{s=-1} for production and \eqn{s=1} for cost function.
}
\section{Functions}{
\itemize{
\item \code{pnormexp()}: distribution function for the normal-exponential distribution.

\item \code{qnormexp()}: quantile function for the normal-exponential distribution.

\item \code{rnormexp()}: random number generation for the normal-exponential distribution.

}}
\examples{
pdf <- dnormexp(x=seq(-3, 3, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
cdf <- pnormexp(q=seq(-3, 3, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
q <- qnormexp(p=seq(0.1, 0.9, by=0.1), mu=1, sigma_v=2, lambda=1/3, s=1)
r <- rnormexp(n=100, mu=1, sigma_v=2, lambda=1/3, s=1)

}
\references{
\itemize{
\item \insertRef{meeusen1977efficiency}{dsfa}
\item \insertRef{kumbhakar2015practitioner}{dsfa}
\item \insertRef{schmidt2020analytic}{dsfa}
\item \insertRef{gradshteyn2014table}{dsfa}
}
}
