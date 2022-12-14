% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_arguments.R
\name{check_arguments}
\alias{check_arguments}
\title{Check inputs arguments}
\usage{
check_arguments(
  x = NULL,
  q = NULL,
  p = NULL,
  n = NULL,
  mu = NULL,
  sigma_v = NULL,
  par_u = NULL,
  s = NULL,
  family = NULL,
  deriv = NULL
)
}
\arguments{
\item{x}{vector of quantiles.}

\item{q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}

\item{mu}{vector of \eqn{\mu}}

\item{sigma_v}{vector of \eqn{\sigma_V}. Must be positive.}

\item{par_u}{vector of parameter of the (in)efficiency term. Must be positive.}

\item{s}{\eqn{s=-1} for production and \eqn{s=1} for cost function.}

\item{family}{\eqn{normhnorm} for normal-halfnormal and \eqn{normexp} for normal-exponential distribution.}

\item{deriv}{derivative of order \code{deriv} of the log density. Available are \code{0},\code{2} and \code{4}.}
}
\value{
Returns a list with five elements:\cr
\code{x} = x,p,q of the recycled length. If the number of observations are provided returns numeric \code{n} \cr
\code{mu} = mu of the recycled length  \cr
\code{sigma_v} = sigma_v of the recycled length \cr
\code{par_u} = par_u of the recycled length \cr
\code{recycle} = logical; if TRUE, the numerical arguments have been recycled to the length of the longest input. \cr
}
\description{
Checks for valid inputs of the probablitiy density, distribution, quantile function and random number generation for the normal-halfnormal and normal-exponential distribution.
}
\details{
Mostly internal function. This function is written to automatically check the inputs and if required to recycle the inputs such that they have the same length. The length is determined by the longest input.
}
\examples{

check_arguments(x=c(-1,1), mu=1, sigma_v=2, par_u=3, s=-1, family="normhnorm", deriv=2)
}
\keyword{internal}
