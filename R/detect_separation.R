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
#' Method for [glm()] that tests for data separation and
#' finds which parameters have infinite maximum likelihood estimates
#' in generalized linear models with binomial responses
#'
#' [detect_separation()] is a method for [glm()]
#' that tests for the occurrence of complete or quasi-complete
#' separation in datasets for binomial response generalized linear
#' models, and finds which of the parameters will have infinite
#' maximum likelihood estimates. [detect_separation()]
#' relies on the linear programming methods developed in Konis (2007).
#'
#' @inheritParams stats::glm.fit
#'
#' @aliases detectSeparation print.detect_separation
#'
#' @param x `x` is a design matrix of dimension `n * p`.
#' @param y `y` is a vector of observations of length `n`.
#' @param control a list of parameters controlling separation
#'     detection. See `detect_separation_control()` for
#'     details.
#' @param start currently not used.
#' @param mustart currently not used.
#' @param etastart currently not used.
#' @param singular.ok logical. If `FALSE`, a singular model is an
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
#' [detect_separation()] determines whether or not the
#' data exhibits (quasi-)complete separation. Then, if separation is
#' detected and the link function \eqn{G(t)} is such that \eqn{\log
#' G(t)} and \eqn{\log (1 - G(t))} are concave, the maximum likelihood
#' estimates has infinite components.
#'
#' [detect_separation()] is a wrapper to the
#' [detect_infinite_estimates()] method. Separation
#' detection, as separation is defined above, takes place using the
#' linear programming methods in Konis (2007) regardless of the link
#' function. The output of those methods is also used to determine
#' which estimates are infinite, unless the link is `"log"`. In the
#' latter case the linear programming methods in Schwendinger et
#' al. (2021) are called to establish if and which estimates are
#' infinite. If the link function is not one of `"logit"`, `"log"`,
#' `"probit"`, `"cauchit"`, `"cloglog"` then a warning is issued.
#'
#' If `separation_type = TRUE` in
#' `detect_separation_control()`, then, whenever
#' separation is detected, [detect_separation()] attempts
#' to distinguish between complete and quasi-complete separation by
#' solving an additional linear program that maximizes the minimum
#' transformed margin. If \eqn{x_i} is the \eqn{i}th row of the model
#' matrix and \eqn{\tilde{y}_i = -1 + 2 y_i}, where \eqn{y_i} is the
#' \eqn{i}th Bernoulli response (the representation to which any
#' binomial data supplied by the user are transformed internally),
#' then \eqn{\bar{x}_i = \tilde{y}_i x_i}. The additional linear
#' program that is solved maximizes \eqn{t} subject to
#' \eqn{\bar{X}\beta \ge t {\bf 1}} and \eqn{-1 \le \beta_j \le
#' 1}. This is equivalent to solving \eqn{\max_\beta \min_i
#' (\bar{x}_i^\top \beta)} subject to \eqn{-1 \le \beta_j \le 1} for
#' all \eqn{j}. A positive optimal value for \eqn{t} implies complete
#' separation, while an optimal value of zero implies
#' quasi-complete separation provided that separation has already been
#' detected. See Konis (2007, Section 1.3 and Chapter 4) for the
#' definitions of separation and the transformed linear programming
#' formulation.
#'
#' The [coefficients()] method extracts a vector of values
#' for each of the model parameters under the following convention:
#' `0` if the maximum likelihood estimate of the parameter is
#' finite, and `Inf` or `-Inf` if the maximum likelihood
#' estimate of the parameter if plus or minus infinity. This
#' convention makes it easy to adjust the maximum likelihood estimates
#' to their actual values by element-wise addition.
#'
#' [detect_separation()] can be passed directly as
#' a method to the [glm()] function. See, examples.
#'
#' [detectSeparation()] is an alias for [detect_separation()].
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
#' [detect_separation()] was designed in 2017 by Ioannis
#' Kosmidis for the \pkg{brglm2} R package, after correspondence with
#' Kjell Konis, and a port of the `separator()` function had been
#' included in \pkg{brglm2} under the permission of Kjell Konis. In
#' 2020, [detect_separation()] and
#' [check_infinite_estimates()] were moved outside
#' **brglm2** into the dedicated **detectseparation** package. Dirk
#' Schumacher authored the `separator_ROI()` function, which
#' depends on the \pkg{ROI} R package and is now the default
#' implementation used for detecting separation. In 2022, Florian
#' Schwendinger authored the `dielb_ROI()` function for detecting
#' infinite estimates in log-binomial regression, and, with Ioannis
#' Kosmidis, they refactored the codebase to properly accommodate for
#' the support of log-binomial regression.
#'
#' @return
#'
#' A list that inherits from class [`detect_separation`],
#' [`glm`] and [`lm`]. A `print` method is provided for
#' [`detect_separation`] objects. If
#' `detect_separation_control(separation_type = TRUE)` is used and
#' separation is detected, then the returned object has a `complete`
#' component, which is `TRUE` for complete separation and `FALSE` for
#' quasi-complete separation. Otherwise, the `complete` component is
#' `NULL`.
#'
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}, Dirk Schumacher `[aut]` \email{mail@dirk-schumacher.net}, Florian Schwendinger `[aut]` \email{FlorianSchwendinger@gmx.at}, Kjell Konis `[ctb]` \email{kjell.konis@me.com}
#'
#' @seealso [glm()], [detect_infinite_estimates()], [check_infinite_estimates()], [brglm2::brglm_fit()]
#'
#' @references
#'
#' Konis K. (2007). *Linear Programming Algorithms for Detecting
#' Separated Data in Binary Logistic Regression
#' Models*. DPhil. University of Oxford.
#' <https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a>
#'
#' Konis K. (2013). safeBinaryRegression: Safe Binary Regression. R
#' package version 0.1-3.
#' <https://CRAN.R-project.org/package=safeBinaryRegression>
#'
#' Kosmidis I. and Firth D. (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71–82. \doi{10.1093/biomet/asaa052}
#'
#' Silvapulle, M. J. (1981).  On the Existence of Maximum Likelihood
#' Estimators for the Binomial Response Models.  *Journal of the Royal
#' Statistical Society. Series B (Methodological)*, **43**, 310–313.
#' <https://www.jstor.org/stable/2984941>
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
#' # If we want to further check for the type of separation (complete
#' # or quasi-complete) we can do
#' endometrial_sep_type <- update(endometrial_sep, separation_type = TRUE)
#' endometrial_sep
#' endometrial_sep_type
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
            warning("`detect_separation()` results may be unreliable for binomial-response GLMs",
                    " with links other than ", paste(shQuote(reliable_links), collapse = ", "))
        }
        if (log_link) {
            message("Data separation in log-binomial models does not necessarily result in infinite estimates")
        }
    }
    else {
        warning("`detect_separation()` has been developed for use with binomial-response GLMs")
    }
    out <- .detect_infinite_estimates(x = x, y = y, weights = weights, start = start,
                                      etastart = etastart,  mustart = mustart,
                                      offset = offset, family = family, control = control,
                                      intercept = control, singular.ok = singular.ok,
                                      log_link = FALSE)
    out$complete <- NULL
    if (isTRUE(out$control$separation_type) && isTRUE(out$outcome)) {
        complete <- separation_type_ROI(x = out$x, y = out$y,
                                        solver = out$control$solver,
                                        tolerance = out$control$tolerance,
                                        solver_control = out$control$solver_control)
        out$complete <- complete$outcome
        if (is.na(complete$outcome)) {
            warning("unexpected result from `separation_type_ROI` with solver: ", out$control$solver, "\n")
        }
    }
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
    if (is.null(x$complete)) {
        str <- "\n"
    } else {
        if (is.na(x$complete)) {
            str <- "(type undetermined)\n"
        }
        else {
            str <- ifelse(isTRUE(x$complete), "(complete)\n", "(quasi-complete)\n")
        }
    }
    cat("Separation:", x$outcome, str)
    if (!is.null(x$coefficients)) {
        cat("Existence of maximum likelihood estimates\n")
        print(coefficients(x))
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}

#' Auxiliary function for the [glm()] interface when
#' `method` is [detect_separation()]
#'
#' Typically only used internally by [detect_separation()]
#' but may be used to construct a `control` argument.
#'
#' @aliases detectSeparationControl
#' @param implementation should the implementation using `ROI` or the
#'     implementation using `lpSolveAPI` be used? Default is `ROI`.
#' @param tolerance maximum absolute variable value from the linear
#'     program, before separation is declared. Default is `1e-04`.
#' @param linear_program should [detect_separation()] solve the
#'     `"primal"` (default) or `"dual"` linear program for separation
#'     detection? Only relevant if `implementation = "lpSolveAPI"`.
#' @param purpose should [detect_separation()] simply `"test"` for
#'     separation or also `"find"` (default) which parameters are
#'     infinite? Only relevant if `implementation = "lpSolveAPI"`.
#' @param solver should the linear program be solved using the
#'     `"lpsolve"` (using the `ROI.plugin.lpsolve` package; default)
#'     or another solver? Alternative solvers are `"glpk"`, `"cbc"`,
#'     `"clp"`, `"cplex"`, `"ecos"`, `"gurobi"`, `"scs"`,
#'     `"symphony"`. If the corresponding package
#'     `ROI.plugin.[solver]` is not installed then the user will be
#'     prompted to install it before continuing.
#' @param solver_control a list with additional control parameters for
#'     the `"solver"`. This is solver specific, so consult the
#'     corresponding documentation. Default is `list()` unless
#'     `solver` is `"alabama"` when the default is `list(start =
#'     rep(0, p))`, where `p` is the number of parameters.
#' @param separation_type logical. Should [detect_separation()]
#'     attempt to distinguish complete from quasi-complete separation
#'     after separation has been detected? If `TRUE` and separation is
#'     detected, then an additional linear program is solved. Default
#'     is `FALSE`.
#'
#'
#' @return
#'
#' A list with the supplied `linear_program`, `solver`,
#' `solver_control`, `purpose`, `tolerance`, `separation_type`,
#' and `implementation`. The returned list is intended to be passed to
#' the `control` argument of [detect_separation()] and
#' [detect_infinite_estimates()].
#'
#' @references
#'
#' Konis K. (2007). *Linear Programming Algorithms for Detecting
#' Separated Data in Binary Logistic Regression Models*. DPhil.
#' University of Oxford.
#' <https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a>
#'
#' @examples
#' data("endometrial", package = "detectseparation")
#' ctrl <- detect_separation_control(separation_type = TRUE)
#' glm(HG ~ NV + PI + EH, data = endometrial,
#'     family = binomial("logit"),
#'     method = "detect_separation",
#'     control = ctrl)
#'
#' @export
detect_separation_control <- function(implementation = c("ROI", "lpSolveAPI"),
                                      solver = "lpsolve",
                                      linear_program = c("primal", "dual"),
                                      purpose = c("find", "test"),
                                      tolerance = 1e-04,
                                      solver_control = list(),
                                      separation_type = FALSE) {
    implementation <- match.arg(implementation)
    purpose <- match.arg(purpose)
    linear_program <- match.arg(linear_program)
    check_ROI_solver(solver)
    list(linear_program = linear_program,
         solver = solver,
         solver_control = solver_control,
         purpose = purpose,
         tolerance = tolerance,
         separation_type = separation_type,
         implementation = implementation)
}
