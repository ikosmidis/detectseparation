#' detectseparation: Methods for Detecting and Checking for Separation
#' and Infinite Maximum Likelihood Estimates
#'
#' \pkg{detectseparation} provides pre-fit and post-fit methods for
#' the detection of separation and of infinite maximum likelihood
#' estimates in binomial response generalized linear models.
#'
#' The key methods are [detect_separation()] and
#' [check_infinite_estimates()].
#'
#' @seealso [detect_separation()], [check_infinite_estimates()]
#'
#' @name detectseparation
#' @import ROI.plugin.lpsolve
#' @importFrom stats coef coefficients gaussian update vcov sd binomial nobs
#' @importFrom graphics matplot
#'
"_PACKAGE"

#' Generic method for checking for infinite estimates
#' @param object a fitted model object (e.g. the result of a
#'     [glm()] call).
#' @param ... other options to be passed to the method.
#'
#' @seealso check_infinite_estimates.glm
#' @export
check_infinite_estimates <- function(object, ...) {
    UseMethod("check_infinite_estimates")
}
