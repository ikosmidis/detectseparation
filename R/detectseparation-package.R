#' detectseparation: Detect or check for separation and infinite maximum likelihood esitmates in generalized linear models with categorical responses
#'
#' [SHORT DESCRIPTION HERE]
#'
#' @docType package
#' @name detectseparation
#' @importFrom stats coef coefficients gaussian update vcov

NULL

#' Generic method for checking for infinite estimates
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call).
#' @param ... other options to be passed to the method.
#' @export
check_infinite_estimates <- function(object, ...) {
    UseMethod("check_infinite_estimates")
}
