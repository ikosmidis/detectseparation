#' detectseparation: Methods for Detecting and Checking for Separation
#' and Infinite Maximum Likelihood Estimates
#'
#' @seealso \code{\link{detect_separation}}, \code{\link{check_infinite_estimates}}
#'
#' @docType package
#' @name detectseparation
#' @import ROI.plugin.lpsolve
#' @importFrom stats coef coefficients gaussian update vcov sd
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
