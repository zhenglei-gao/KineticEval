#' Kinetic Evaluations using NLS, IRLS and MCMC.
#'
#' Roxygen is a Doxygen-like documentation system for R; allowing
#' in-source specification of Rd files, collation and namespace
#' directives.
#'
#' If you have existing Rd files, check out the \code{Rd2roxygen} package
#' for a convenient way of converting Rd files to roxygen comments.
#'
#' \tabular{ll}{
#' Package: \tab KineticEval\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1-1\cr
#' Date: \tab 2012-03-15\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @author
#' Zhenglei Gao \email{zhenglei.gao@@bayer.com}, Horatio Meyer and Walter Schmitt.
#'
#' Maintainer: Zhenglei Gao \email{zhenglei.gao@@bayer.com}
#' @name KineticEval-package
#' @docType package
#' @title Kinetic Evaluations in R
#' @keywords package
#' @import mkin
#' @importFrom deSolve ode
#' @importFrom FME modFit modMCMC
#' @importFrom Rsolnp solnp
#' @importFrom minqa bobyqa
#' @importFrom coda traceplot as.mcmc
#' @importFrom optimx optimx
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @note This is based on the package \code{mkin} by Johannes Reinken.
#'
#' You are welcome to report other serious issues and bugs via
#' \url{http://dropbox.com}.
#' @examples
#' ## see the package vignette: vignette('KineticEval')
NULL