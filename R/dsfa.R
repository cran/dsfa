#' dsfa: Distributional Stochastic Frontier Analysis
#'
#' @docType package
#' @name dsfa
#' @author \itemize{
#' \item Rouven Schmidt  \email{rouven.schmidt@tu-clausthal.de}
#' }
#'
#' @description
#' The \code{dsfa} package implements the specification, estimation and prediction of distributional stochastic frontier models via \code{mgcv}.
#' The basic distributional stochastic frontier model is given by: \deqn{Y_n = \eta^\mu(\boldsymbol{x}_n^\mu) + V_n + s \cdot U_n } where \eqn{n \in \{1,2,...,N\}}, \eqn{V_n} and \eqn{U_n} are the noise and (in)efficiency respectively.
#' \itemize{
#' \item For \eqn{s=-1}, \eqn{\eta^\mu(\cdot)} is the production function and \eqn{\boldsymbol{x}_n^\mu} are the log inputs. Alternatively, if \eqn{s=1}, \eqn{\eta^\mu(\cdot)} is the cost function and \eqn{\boldsymbol{x}_n^\mu} are the log cost.
#' \item The noise is represented as \eqn{V_n \sim N(0,\sigma_{Vn}^2)}, where \eqn{\sigma_{Vn}=\exp(\eta^{\sigma_{V}}(\boldsymbol{x}_n^{\sigma_{V}}))}. Here, \eqn{\boldsymbol{x}_n^{\sigma_{V}}} are the observed covariates which influence the parameter of the noise.
#' \item The inefficiency is represented either as \eqn{U_n \sim HN(\sigma_{Un}^2)} or as \eqn{U_n \sim Exp(\lambda_{n})}, where \eqn{\sigma_{Un}=\exp(\eta^{\sigma_{Un}}(\boldsymbol{x}_n^{\sigma_{U}}))} and \eqn{\lambda_{n}=\exp(\eta^{\lambda_{n}}(\boldsymbol{x}_n^{\lambda}))}. Here, \eqn{\boldsymbol{x}_n^{\sigma_{U}}} or \eqn{\boldsymbol{x}_n^{\lambda}} are the observed covariates which influence the respective parameter of the (in)efficiency.
#'  }
#' If \eqn{U_n \sim HN(\sigma_{Un}^2)}, then: \deqn{Y \sim normhnorm(\mu_n=\eta^\mu(\boldsymbol{x}_n^\mu), \sigma_{Vn}=\exp(\eta^{\sigma_{V}}(\boldsymbol{x}_n^{\sigma_{V}})), \sigma_{Un}=\exp(\eta^{\sigma_{U}}(\boldsymbol{x}_n^{\sigma_{U}})), s=s)}
#' Alternatively, if \eqn{U_n \sim Exp(\lambda_{n})} then: \deqn{Y \sim normexp(\mu_n=\eta^\mu(\boldsymbol{x}_n^\mu), \sigma_{Vn}=\exp(\eta^{\sigma_{V}}(\boldsymbol{x}_n^{\sigma_{V}})), \lambda_{n}=\exp(\eta^{\lambda}(\boldsymbol{x}_n^{\lambda})), s=s)}.
#'
#' The package fits \eqn{\eta^\mu(\boldsymbol{x}_n^\mu), \eta^{\sigma_{V}}(\boldsymbol{x}_n^{\sigma_{V}})} and \eqn{\eta^{\lambda}(\boldsymbol{x}_n^{\lambda})} or \eqn{\eta^{\sigma_{U}}(\boldsymbol{x}_n^{\sigma_{U}})}.
#'
#' @details
#' The \code{mgcv} packages provides a framework for fitting distributional regression models.
#' The formulae in \code{\link[mgcv:gam]{gam}} allow for smooth terms utilizing the function \code{\link[mgcv:s]{s}}.
#' These may be
#' \itemize{
#' \item linear effects
#' \item non-linear effects which can be modeled via penalized regression splines, e.g. \code{\link[mgcv:p.spline]{p.spline}}, \code{\link[mgcv:tprs]{tprs}}
#' \item random effects, \code{\link[mgcv:random.effects]{random.effects}},
#' \item spatial effects which can be modeled via \code{\link[mgcv:mrf]{mrf}}.
#' }
#' An overview is provided at \code{\link[mgcv:smooth.terms]{smooth.terms}}. The functions \code{\link[mgcv:gam]{gam}}, \code{\link[mgcv:predict.gam]{predict.gam}} and \code{\link[mgcv:plot.gam]{plot.gam}}, are alike to the basic S functions.
#' A number of other functions such as \code{\link[mgcv:summary.gam]{summary.gam}}, \code{\link[mgcv:residuals.gam]{residuals.gam}} and \code{\link[mgcv:anova.gam]{anova.gam}} are also provided, for extracting information from a fitted \code{\link[mgcv:gamObject]{gamOject}}.
#'
#' The main functions are:
#' \itemize{
#' \item \code{\link{normhnorm}}  Object which can be used to fit a normal-halfnormal stochastic frontier model with the \code{mgcv} package.
#' \item \code{\link{normexp}}  Object which can be used to fit a normal-exponential stochastic frontier model with the \code{mgcv} package.
#' \item \code{\link{comperr_mv}}  Object which can be used to fit a multivariate stochastic frontier model with the \code{mgcv} package.
#' \item \code{\link{elasticity}}  Calculates and plots the elasticity of a smooth function.
#' \item \code{\link{efficiency}}  Calculates the expected technical (in)efficiency index \eqn{E[u|\epsilon]} or \eqn{E[\exp(-u)|\epsilon]}.
#' }
#'
#' @references
#' \itemize{
#' \item \insertRef{schmidt2022mvdsfm}{dsfa}
#' \item \insertRef{wood2017generalized}{dsfa}
#' \item \insertRef{kumbhakar2015practitioner}{dsfa}
#' \item \insertRef{schmidt2020analytic}{dsfa}
#' }
NULL