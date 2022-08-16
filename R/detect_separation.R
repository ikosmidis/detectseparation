# Copyright (C) 2017- Ioannis Kosmidis, Dirk Schumacher

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title Detect Separation
#'
#' @description
#'
#' Method for \code{\link{glm}} that tests for data separation and
#' finds which parameters have infinite maximum likelihood estimates
#' in generalized linear models with binomial responses
#'
#' \code{\link{detect_separation}()} is a method for \code{\link{glm}}
#' that tests for the occurrence of complete or quasi-complete
#' separation in datasets for binomial response generalized linear
#' models, and finds which of the parameters will have infinite
#' maximum likelihood estimates. \code{\link{detect_separation}()}
#' relies on the linear programming methods developed in Konis (2007).
#'
#' @inheritParams stats::glm.fit
#'
#' @aliases detectSeparation print.detect_separation
#'
#' @param x \code{x} is a design matrix of dimension \code{n * p}.
#' @param y \code{y} is a vector of observations of length \code{n}.
#' @param control a list of parameters controlling separation
#'     detection. See \code{\link{detect_separation_control}()} for
#'     details.
#' @param start currently not used.
#' @param mustart currently not used.
#' @param etastart currently not used.
#' @param singular.ok logical. If \code{FALSE}, a singular model is an
#'     error.
#'
#' @details
#'
#' Following the definitions in Albert and Anderson (1984), the data
#' for a binomial-response generalized linear model with logistic link
#' exhibit quasi-complete separation if there exists a non-zero
#' parameter vector \eqn{\beta} such that \eqn{X^0 \beta \le 0} and
#' \eqn{X^1 \beta \ge 0}, where \eqn{X^0} and \eqn{X^1} are the
#' matrices formed by the rows of the model matrix $X$ corresponding
#' to zero and non-zero responses, respectively. The data exhibits
#' complete separation if there exists a parameter vector \eqn{\beta} such
#' that the aforementioned conditions are satisfied with strict
#' inequalities. If there are no vectors \eqn{\beta} that can satisfy the
#' conditions, then the data points are said to overlap.
#'
#' If the inverse link function \eqn{G(t)} of a generalized linear
#' model with binomial responses is such that \eqn{\log G(t)} and
#' \eqn{\log (1 - G(t))} are concave and the model has an intercept
#' parameter, then overlap is a necessary and sufficient condition for
#' the maximum likelihood estimates to be finite (see Silvapulle, 1981
#' for a proof). Such link functions are, for example, the logit,
#' probit and complementary log-log.
#'
#' \code{\link{detect_separation}()} determines whether or not the
#' data exhibits (quasi-)complete separation. Then, if separation is
#' detected and the link function \eqn{G(t)} is such that \eqn{\log
#' G(t)} and \eqn{\log (1 - G(t))} are concave, the maximum likelihood
#' estimates has infinite components.
#'
#' \code{\link{detect_separation}()} is a wrapper to the
#' \code{\link{detect_infinite_estimates}()} method. Separation
#' detection, as separation is defined above, takes place using the
#' linear programming methods in Konis (2007) regardless of the link
#' function. The output of those methods is also used to determine
#' which estimates are infinite, unless the link is "log". In the
#' latter case the linear programming methods in Schwendinger et
#' al. (2021) are called to establish if and which estimates are
#' infinite. If the link function is not one of `"logit"`, `"log"`,
#' `"probit"`, `"cauchit"`, `"cloglog"` then a warning is issued.
#'
#' The \code{\link{coefficients}} method extracts a vector of values
#' for each of the model parameters under the following convention:
#' \code{0} if the maximum likelihood estimate of the parameter is
#' finite, and \code{Inf} or \code{-Inf} if the maximum likelihood
#' estimate of the parameter if plus or minus infinity. This
#' convention makes it easy to adjust the maximum likelihood estimates
#' to their actual values by element-wise addition.
#'
#' \code{\link{detect_separation}()} can be passed directly as
#' a method to the \code{\link{glm}} function. See, examples.
#'
#' \code{detectSeparation}() is an alias for \code{detect_separation}().
#'
#' @note
#'
#' For the definition of complete and quasi-complete separation, see
#' Albert and Anderson (1984). Kosmidis and Firth (2021) prove that
#' the reduced-bias estimator that results by the penalization of the
#' logistic regression log-likelihood by Jeffreys prior takes always
#' finite values, even when some of the maximum likelihood estimates
#' are infinite. The reduced-bias estimates can be computed using the
#' \pkg{brglm2} R package.
#'
#' \code{\link{detect_separation}} was designed in 2017 by Ioannis
#' Kosmidis for the **brglm2** R package, after correspondence with
#' Kjell Konis, and a port of the \code{separator} function had been
#' included in **brglm2** under the permission of Kjell Konis. In
#' 2020, \code{\link{detect_separation}} and
#' \code{\link{check_infinite_estimates}} were moved outside
#' **brglm2** into the dedicated **detectseparation** package. Dirk
#' Schumacher authored the \code{separator_ROI} function, which
#' depends on the **ROI** R package and is now the default
#' implementation used for detecting separation. In 2022, Florian
#' Schwendinger authored the \code{dielb_ROI} function for detecting
#' infinite estimates in log-binomial regression, and, with Ioannis
#' Kosmidis, they refactored the codebase to properly accommodate for
#' the support of log-binomial regression.
#'
#' @return
#'
#' A list that inherits from class \code{detect_separation},
#' \code{glm} and \code{lm}. A \code{print} method is provided for
#' \code{detect_separation} objects.
#'
#'
#' @author Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}, Dirk Schumacher [aut] \email{mail@dirk-schumacher.net}, Florian Schwendinger [aut] \email{FlorianSchwendinger@gmx.at}, Kjell Konis [ctb] \email{kjell.konis@me.com}
#'
#' @seealso \code{\link{glm.fit}} and \code{\link{glm}}, \code{\link{detect_infinite_estimates}}, \code{\link{check_infinite_estimates}}, \code{\link[brglm2]{brglm_fit}}
#'
#' @references
#'
#' Konis K. (2007). *Linear Programming Algorithms for Detecting
#' Separated Data in Binary Logistic Regression
#' Models*. DPhil. University of Oxford.
#' \url{https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}
#'
#' Konis K. (2013). safeBinaryRegression: Safe Binary Regression. R
#' package version 0.1-3.
#' \url{https://CRAN.R-project.org/package=safeBinaryRegression}
#'
#' Kosmidis I. and Firth D. (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71–82. \doi{10.1093/biomet/asaa052}
#'
#' Silvapulle, M. J. (1981).  On the Existence of Maximum Likelihood
#' Estimators for the Binomial Response Models.  *Journal of the Royal
#' Statistical Society. Series B (Methodological)*, **43**, 310–313.
#' \url{https://www.jstor.org/stable/2984941}
#'
#' Schwendinger, F., Grün, B. & Hornik, K. (2021). A comparison of
#' optimization solvers for log binomial regression including conic
#' programming.  *Computational Statistics*, **36**,
#' 1721–1754. \doi{10.1007/s00180-021-01084-5}
#'
#'
#' @examples
#'
#' # endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "detectseparation")
#' endometrial_sep <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                        family = binomial("logit"),
#'                        method = "detect_separation")
#' endometrial_sep
#' # The maximum likelihood estimate for NV is infinite
#' summary(update(endometrial_sep, method = "glm.fit"))
#'
#' \donttest{
#' # Example inspired by unpublished microeconometrics lecture notes by
#' # Achim Zeileis https://eeecon.uibk.ac.at/~zeileis/
#' # The maximum likelihood estimate of sourhernyes is infinite
#' if (requireNamespace("AER", quietly = TRUE)) {
#'     data("MurderRates", package = "AER")
#'     murder_sep <- glm(I(executions > 0) ~ time + income +
#'                       noncauc + lfp + southern, data = MurderRates,
#'                       family = binomial(), method = "detect_separation")
#'     murder_sep
#'     # which is also evident by the large estimated standard error for NV
#'     murder_glm <- update(murder_sep, method = "glm.fit")
#'     summary(murder_glm)
#'     # and is also revealed by the divergence of the NV column of the
#'     # result from the more computationally intensive check
#'     plot(check_infinite_estimates(murder_glm))
#'     # Mean bias reduction via adjusted scores results in finite estimates
#'     if (requireNamespace("brglm2", quietly = TRUE))
#'         update(murder_glm, method = brglm2::brglm_fit)
#' }
#' }
#' @export
detect_separation <- function(x, y, weights = NULL,
                              start = NULL, etastart = NULL,  mustart = NULL,
                              offset = NULL, family = gaussian(),
                              control = list(), intercept = TRUE, singular.ok = TRUE) {
    log_link <- isTRUE(family$link == "log")
    if (isTRUE(family$family == "binomial")) {
        reliable_links <- c("logit", "log", "probit", "cauchit", "cloglog")
        if (!isTRUE(family$link %in% reliable_links)) {
            warning("`detect_separation` results may be unreliable for binomial-response GLMs",
                    " with links other than ", paste(shQuote(reliable_links), collapse = ", "))
        }
        if (log_link) {
            message("Data separation in log-binomial models does not necessarily result in infinite estimates")
        }
    }
    else {
        warning("`detect_separation` has been developed for use with binomial-response GLMs")
    }
    out <- .detect_infinite_estimates(x = x, y = y, weights = weights, start = start,
                                      etastart = etastart,  mustart = mustart,
                                      offset = offset, family = family, control = control,
                                      intercept = control, singular.ok = singular.ok,
                                      log_link = FALSE)
    if (log_link) {
        # test for existence using the linear program in Schwendinger et al (2021)
        out$coefficients <- .detect_infinite_estimates(x = x, y = y, weights = weights, start = start,
                                                       etastart = etastart,  mustart = mustart,
                                                       offset = offset, family = family, control = control,
                                                       intercept = control, singular.ok = singular.ok,
                                                       log_link = TRUE)$coefficients
    }
    class(out) <- out$class <- "detect_separation"
    out
}

#' @method print detect_separation
#' @export
print.detect_separation <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
    cat("Implementation:", x$control$implementation, "| ")
    if (identical(x$control$implementation, "ROI")) {
        cat("Solver:", x$control$solver, "\n")
    }
    else {
        cat("Linear program:", x$control$linear_program, "| Purpose:", x$control$purpose, "\n")
    }
    cat("Separation:", x$outcome, "\n")
    if (!is.null(x$coefficients)) {
        cat("Existence of maximum likelihood estimates\n")
        print(coefficients(x))
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}

#' Auxiliary function for the \code{\link{glm}} interface when
#' \code{method} is \code{\link{detect_separation}}.
#'
#' Typically only used internally by \code{\link{detect_separation}}
#' but may be used to construct a \code{control} argument.
#'
#' @aliases detectSeparationControl
#' @param implementation should the implementation using \code{ROI} or
#'     the implementation using \code{lpSolveAPI} be used? Default is
#'     \code{ROI}.
#' @param tolerance maximum absolute variable value from the linear
#'     program, before separation is declared. Default is
#'     \code{1e-04}.
#' @param linear_program should \code{\link{detect_separation}} solve
#'     the \code{"primal"} (default) or \code{"dual"} linear program
#'     for separation detection? Only relevant if \code{implementation
#'     = "lpSolveAPI"}.
#' @param purpose should \code{\link{detect_separation}} simply
#'     \code{"test"} for separation or also \code{"find"} (default)
#'     which parameters are infinite? Only relevant if
#'     \code{implementation = "lpSolveAPI"}.
#' @param solver should the linear program be solved using the
#'     \code{"lpsolve"} (using the \pkg{ROI.plugin.lpsolve} package;
#'     default) or another solver? Alternative solvers are
#'     \code{"glpk"}, \code{"cbc"}, \code{"clp"}, \code{"cplex"},
#'     \code{"ecos"}, \code{"gurobi"}, \code{"scs"},
#'     \code{"symphony"}. If \pkg{ROI.plugin.[solver]} is not
#'     installed then the user will be prompted to install it before
#'     continuing.
#' @param solver_control a list with additional control parameters for
#'     the \code{"solver"}. This is solver specific, so consult the
#'     corresponding documentation. Default is \code{list()} unless
#'     \code{solver} is \code{"alabama"} when the default is \code{list(start
#'     = rep(0, p))}, where p is the number of parameters.
#'
#' @return
#'
#' A list with the supplied \code{linear_program}, \code{solver},
#' \code{solver_control}, \code{purpose}, \code{tolerance},
#' \code{implementation}, and the matched \code{separator} function
#' (according to the value of \code{implementation}).
#'
#' @export
detect_separation_control <- function(implementation = c("ROI", "lpSolveAPI"),
                                      solver = "lpsolve",
                                      linear_program = c("primal", "dual"),
                                      purpose = c("find", "test"),
                                      tolerance = 1e-04,
                                      solver_control = list()) {
    implementation <- match.arg(implementation)
    purpose <- match.arg(purpose)
    linear_program <- match.arg(linear_program)
    separator <- getNamespace("detectseparation")[[paste("separator", implementation, sep = "_")]]
    check_ROI_solver(solver)
    list(linear_program = linear_program,
         solver = solver,
         solver_control = solver_control,
         purpose = purpose,
         tolerance = tolerance,
         separator = separator,
         implementation = implementation)
}
