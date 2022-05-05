#' detectseparation: Methods for Detecting and Checking for Separation
#' and Infinite Maximum Likelihood Estimates
#'
#' \pkg{detectseparation} provides pre-fit and post-fit methods for
#' the detection of separation and of infinite maximum likelihood
#' estimates in binomial response generalized linear models.
#'
#' The key methods are \code{\link{detect_separation}} and
#' \code{\link{check_infinite_estimates}}.
#'
#' @seealso \code{\link{detect_separation}}, \code{\link{check_infinite_estimates}}
#'
#' @docType package
#' @name detectseparation
#' @import ROI.plugin.lpsolve
#' @importFrom stats coef coefficients gaussian update vcov sd binomial nobs
#' @importFrom graphics matplot
#'
NULL

#' Generic method for checking for infinite estimates
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call).
#' @param ... other options to be passed to the method.
#'
#' @seealso check_infinite_estimates.glm
#' @export
check_infinite_estimates <- function(object, ...) {
    UseMethod("check_infinite_estimates")
}
