# Copyright (C) 2017-2021 Ioannis Kosmidis, Dirk Schumacher

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
#' According to the definitions in Albert and Anderson (1984), the
#' data exhibits quasi-complete separation if there exists a parameter
#' vector $\beta \ne 0$ such that $X^0 \beta \le 0$ and $X^1 \beta \ge
#' 0$, where $X^0$ and $X^1$ are the matrices formed by the rows of
#' the model matrix $X$ corresponding to zero and non-zero responses,
#' respectively. The data exhibits complete separation if there exists
#' a parameter vector $\beta$ such that the aforementioned conditions
#' are satisfied with strict inequalities. If there are no vectors
#' $\beta$ that can satisfy the conditions, then the data points are
#' said to overlap.
#'
#' If the inverse link function $G(t)$ of a generalized linear models
#' with binomial responses is such that $\log G(t)$ and $\log (1 -
#' G(t))$ are convave and the model has an intercept parameter, then
#' overlap is a necessary and sufficient condition for the maximum
#' likelihood estimates to be finite (see Silvapulle, 1981 for a
#' proof). Such link functions are, for example, the logit, probit and
#' complementary log-log.
#'
#' \code{\link{detect_separation}()} determines whether or not the
#' data exhibits (quasi-)complete separation. Then, if separation is
#' detected and the link function $G(t)$ is such that $\log G(t)$ and
#' $\log (1 - G(t))$ are concave, the maximum likelihood estimates
#' has infinite components.
#'
#' \code{\link{detect_separation}()} is a wrapper to the
#' \code{separator_ROI} function and \code{separator_lpSolveAPI}
#' function (a modified version of the \code{separator} function from
#' the **safeBinaryRegression** R
#' package). \code{\link{detect_separation}()} can be passed directly as
#' a method to the \code{\link{glm}} function. See, examples.
#'
#' The \code{\link{coefficients}} method extracts a vector of values
#' for each of the model parameters under the following convention:
#' \code{0} if the maximum likelihood estimate of the parameter is
#' finite, and \code{Inf} or \code{-Inf} if the maximum likelihood
#' estimate of the parameter if plus or minus infinity. This
#' convention makes it easy to adjust the maximum likelihood estimates
#' to their actual values by element-wise addition.
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
#' included in **brglm2** under the permission of Kjell Konis.
#'
#' In 2020, \code{\link{detect_separation}} and
#' \code{\link{check_infinite_estimates}} were moved outside
#' **brglm2** into the dedicated **detectseparation** package. Dirk Schumacher
#' authored the \code{separator_ROI} function, which depends on the
#' **ROI** R package and is now the default implementation used for
#' detecting separation.
#'
#' @return
#'
#' A list that inherits from class \code{detect_separation},
#' \code{glm} and \code{lm}. A \code{print} method is provided for
#' \code{detect_separation} objects.
#'
#'
#' @author Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}, Dirk Schumacher [aut] \email{mail@dirk-schumacher.net}, Kjell Konis [ctb] \email{kjell.konis@me.com}
#'
#' @seealso \code{\link{glm.fit}} and \code{\link{glm}}, \code{\link{check_infinite_estimates}}, \code{\link[brglm2]{brglm_fit}},
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
#' models. *Biometrika*, **108**, 71–82
#'
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' \url{http://www.jstor.org/stable/2984941}

#' @examples
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "detectseparation")
#' endometrial_sep <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                        family = binomial("logit"),
#'                        method = "detect_separation")
#' endometrial_sep
#' ## The maximum likelihood estimate for NV is infinite
#' summary(update(endometrial_sep, method = "glm.fit"))
#'
#' \donttest{
#' ## Example inspired by unpublished microeconometrics lecture notes by
#' ## Achim Zeileis https://eeecon.uibk.ac.at/~zeileis/
#' ## The maximum likelihood estimate of sourhernyes is infinite
#' if (requireNamespace("AER", quietly = TRUE)) {
#'     data("MurderRates", package = "AER")
#'     murder_sep <- glm(I(executions > 0) ~ time + income +
#'                       noncauc + lfp + southern, data = MurderRates,
#'                       family = binomial(), method = "detect_separation")
#'     murder_sep
#'     ## which is also evident by the large estimated standard error for NV
#'     murder_glm <- update(murder_sep, method = "glm.fit")
#'     summary(murder_glm)
#'     ## and is also reveal by the divergence of the NV column of the
#'     ## result from the more computationally intensive check
#'     plot(check_infinite_estimates(murder_glm))
#'     ## Mean bias reduction via adjusted scores results in finite estimates
#'     if (requireNamespace("brglm2", quietly = TRUE))
#'         update(murder_glm, method = brglm2::brglm_fit)
#' }
#' }
#' @export
detect_separation <- function(x, y, weights = rep.int(1, nobs),
                              start = NULL, etastart = NULL,  mustart = NULL,
                              offset = rep.int(0, nobs), family = gaussian(),
                              control = list(), intercept = TRUE, singular.ok = TRUE) {
    if (isTRUE(family$family != "binomial")) {
        warning("`detect_separation` has been developed for use with binomial-response GLMs")
    }
    reliable_links <- c("logit", "log", "probit", "cauchit", "cloglog")
    if (!isTRUE(family$link %in% reliable_links)) {
        warning("`detect_separation` results may be unreliable for binomial-response GLMs",
                " with links other than ", paste(shQuote(reliable_links), collapse = ", "))
    }
    if (isTRUE(family$link == "log")) {
        warning("Data separation in log-binomial models does not necessarily result in infinite estimates")
    }
    control <- do.call("detect_separation_control", control)
    separator <- control$separator
    ## ensure x is a matrix
    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
    ##
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missingOffset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    ## Initialize as prescribed in family
    eval(family$initialize)
    if (EMPTY) {
        out <- list(separation = FALSE)
    }
    else {
        if (control$solver == "alabama" & is.null(control$solver_control$start)) {
            control$solver_control$start <- rep(0, nvars)
        }
        ## as in brglmFit
        boundary <- converged <- FALSE
        ## Detect aliasing
        qrx <- qr(x)
        rank <- qrx$rank
        is_full_rank <- rank == nvars
        if (!singular.ok && !is_full_rank) {
            stop("singular fit encountered")
        }
        if (!isTRUE(is_full_rank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            X_all <- x
            x <- x[, -aliased]
            nvars_all <- nvars
            nvars <- ncol(x)
            betas_names_all <- betas_names
            betas_names <- betas_names[-aliased]
        }
        else {
            nvars_all <- nvars
            betas_names_all <- betas_names
        }
        betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
        ## Observations with zero weight do not enter calculations so ignore
        keep <- weights > 0
        x <- x[keep, , drop = FALSE]
        y <- y[keep]
        ## Reshape data set: keep 0 and 1, and replace anything in (0,
        ## 1) with one zero and one 1
        ones <- y == 1
        zeros <- y == 0
        non_boundary <- !(ones | zeros)
        x <- x[c(which(ones), which(zeros), rep(which(non_boundary), 2)), , drop = FALSE]
        y <- c(y[ones], y[zeros], rep(c(0., 1.), each = sum(non_boundary)))
        ## Run linear program
        out <- separator(x = x, y = y,
                         linear_program = control$linear_program,
                         purpose = control$purpose,
                         tolerance = control$tolerance,
                         solver = control$solver,
                         solver_control = control$solver_control)
        if (is.na(out$separation)) {
            if (identical(control$implementation, "ROI")) {
                warning("unexpected result from implementation ", control$implementation, " with solver: ", control$solver, "\n")
            }
            else {
                warning("unexpected result from implementation ", control$implementation, " with linear_program: ", control$linear_program, " and purpose: ", control$purpose, "\n")
            }
        }
        if (is.null(out$beta)) {
            betas_all <- NULL
        }
        else {
            betas <- out$beta
            names(betas) <- betas_names
            inds <- abs(betas) < control$tolerance
            betas <- Inf * betas
            betas[inds] <- 0
            betas_all[betas_names] <- betas
        }
        out <- list(x = x,
                    y = y,
                    coefficients = betas_all,
                    separation = out$separation)
    }
    out$control <- control
    out$class <- "detect_separation"
    class(out) <- "detect_separation"
    return(out)
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
    separator <- match.fun(paste("separator", implementation, sep = "_"))
    check_ROI_solver(solver)
    list(linear_program = linear_program,
         solver = solver,
         solver_control = solver_control,
         purpose = purpose,
         tolerance = tolerance,
         separator = separator,
         implementation = implementation)
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
    cat("Separation:", x$separation, "\n")
    if (!is.null(x$coefficients)) {
        cat("Existence of maximum likelihood estimates\n")
        print(coefficients(x))
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}
