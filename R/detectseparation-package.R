#' detectseparation: Methods for Detecting and Checking for Separation
#' and Infinite Maximum Likelihood Estimates
#'
#' \pkg{detectseparation} provides pre-fit and post-fit methods for
#' the detection of separation and of infinite maximum likelihood
#' estimates in binomial response generalized linear models.
#'
#' The key methods are [detect_separation()], [detect_infinite_estimates()] and
#' [check_infinite_estimates()].
#'
#' #' @details
#' # Declaration of AI-Assisted Development
#'
#' Generative artificial intelligence tools were used to assist with
#' selected aspects of software implementation, testing and
#' documentation.
#'
#' The use of these tools was limited to supporting the development
#' process. The scientific methodology, computational objectives, and
#' final decisions concerning the package remain the responsibility of
#' its human authors. AI-assisted material is subject to the same
#' review, testing, and quality standards as other contributions to
#' the package.
#'
#' This declaration is provided in the interest of transparency and
#' reproducibility. It does not modify the licensing terms of the
#' software.
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
