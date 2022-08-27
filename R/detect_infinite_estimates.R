# Copyright (C) 2022- Florian Schwendinger, Ioannis Kosmidis

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

#' @title Detect Infinite Estimates
#'
#' @description
#'
#' Method for \code{\link{glm}} that detects infinite components in
#' the maximum likelihood estimates of generalized linear models
#' with binomial responses.
#'
#' @inheritParams stats::glm.fit
#'
#' @aliases detectInfiniteEstimates
#'
#' @param x \code{x} is a design matrix of dimension \code{n * p}.
#' @param y \code{y} is a vector of observations of length \code{n}.
#' @param control a list of parameters controlling separation
#'     detection. See \code{\link{detect_separation_control}} for
#'     details.
#' @param start currently not used.
#' @param mustart currently not used.
#' @param etastart currently not used.
#' @param singular.ok logical. If \code{FALSE}, a singular model is an
#'     error.
#'
#'
#' @details
#'
#' For binomial-response generalized linear models with \code{"log"}
#' link, separated data allocations do not necessarily lead to
#' infinite maximum likelihood estimates. For this reason, for models
#' with the \code{"log"} link
#' \code{\link{detect_infinite_estimates}()} relies on an alternative
#' linear optimization model developed in Schwendinger et al. (2021),
#' and for all the other supported links it relies on the linear
#' programming methods developed in Konis (2007). See
#' \code{\link{detect_separation}()} for definitions and details.
#'
#' \code{\link{detect_infinite_estimates}()} is a wrapper to the
#' functions \code{separator_ROI()}, \code{separator_lpSolveAPI()} (a
#' modified version of the \code{separator()} function from the
#' **safeBinaryRegression** R package), and \code{dielb_ROI()}.
#'
#' The \code{\link{coefficients}()} method extracts a vector of values
#' for each of the model parameters under the following convention:
#' \code{0} if the maximum likelihood estimate of the parameter is
#' finite, and \code{Inf} or \code{-Inf} if the maximum likelihood
#' estimate of the parameter if plus or minus infinity. This
#' convention makes it easy to adjust the maximum likelihood estimates
#' to their actual values by element-wise addition.
#'
#' \code{\link{detect_infinite_estimates}()} can be passed directly as
#' a method to the \code{\link{glm}} function. See, examples.
#'
#' \code{detectInfiniteEstimates}() is an alias for
#' \code{detect_infinite_estimates}().
#'
#' @author Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}, Florian Schwendinger [aut] \email{FlorianSchwendinger@gmx.at}, Dirk Schumacher [aut] \email{mail@dirk-schumacher.net}, Kjell Konis [ctb] \email{kjell.konis@me.com}
#'
#' @seealso \code{\link{glm.fit}} and \code{\link{glm}}, \code{\link{detect_separation}}, \code{\link{check_infinite_estimates}}, \code{\link[brglm2]{brglm_fit}}
#'
#' @references
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' \url{https://www.jstor.org/stable/2984941}
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
#' Schwendinger, F., Grün, B. & Hornik, K. (2021). A comparison of
#' optimization solvers for log binomial regression including conic
#' programming.  *Computational Statistics*, **36**,
#' 1721–1754. \doi{10.1007/s00180-021-01084-5}
#'
#'
#' @examples
#' # The classical example given in Silvapulle (1981) can be utilized
#' # to show that for the Log-Binomial model there exist data allocations
#' # which are separated but produce finite estimates.
#' data("silvapulle1981", package = "detectseparation")
#'
#' # Since the data is separated the MLE does not exist for the logit link.
#' glm(y ~ ghqs, data = silvapulle1981, family = binomial(),
#'     method = "detect_infinite_estimates")
#'
#' # However, for the log link all components of the MLE are finite.
#' glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
#'     method = "detect_infinite_estimates")
#' glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"), start = c(-1, 0))
#'
#' @export
detect_infinite_estimates <- function(x, y, weights = NULL,
                                      start = NULL, etastart = NULL,  mustart = NULL,
                                      offset = NULL, family = gaussian(),
                                      control = list(), intercept = TRUE, singular.ok = TRUE) {
    if (isTRUE(family$family != "binomial")) {
        warning("`detect_infinite_estimates` has been developed for use with binomial-response GLMs")
    }
    nobs <- NROW(y)
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    out <- .detect_infinite_estimates(x = x, y = y, weights = weights, start = start,
                                      etastart = etastart,  mustart = mustart,
                                      offset = offset, family = family, control = control,
                                      intercept = control, singular.ok = singular.ok,
                                      log_link = isTRUE(family$link == "log"))
    class(out) <- out$class <- "detect_infinite_estimates"
    out
}

#' @method print detect_infinite_estimates
#' @export
print.detect_infinite_estimates <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
    cat("Implementation:", x$control$implementation, "| ")
    if (identical(x$control$implementation, "ROI")) {
        cat("Solver:", x$control$solver, "\n")
    }
    else {
        cat("Linear program:", x$control$linear_program, "| Purpose:", x$control$purpose, "\n")
    }
    cat("Infinite estimates:", x$outcome, "\n")
    if (!is.null(x$coefficients)) {
        cat("Existence of maximum likelihood estimates\n")
        print(coefficients(x))
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}

## Workhorse function
.detect_infinite_estimates <- function(x, y, weights = NULL,
                                       start = NULL, etastart = NULL,  mustart = NULL,
                                       offset = NULL, family = gaussian(),
                                       control = list(), intercept = TRUE, singular.ok = TRUE,
                                       log_link = FALSE) {
    control <- do.call("detect_separation_control", control)
    lp <- if (log_link) dielb_ROI else control$separator
    # ensure x is a matrix
    x <- as.matrix(x)
    betas_names_all <- betas_names <- if (is.null(colnames(x))) make.names(seq_len(NCOL(x))) else colnames(x)
    #
    nobs <- NROW(y)
    nvars <- ncol(x)
    if (nvars == 0) {
        return(list(outcome = FALSE, control = control))
    }
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missingOffset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    # Initialize as prescribed in family
    eval(family$initialize)
    if (control$solver == "alabama" & is.null(control$solver_control$start)) {
        control$solver_control$start <- rep(0, nvars)
    }
    # Detect aliasing
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
        betas_names <- betas_names[-aliased]
    }
    betas_all <- structure(rep(NA_real_, length(betas_names_all)), .Names = betas_names_all)
    # Observations with zero weight do not enter calculations so ignore
    keep <- weights > 0
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    # Reshape data set: keep 0 and 1, and replace anything in (0,
    # 1) with one zero and one 1
    ones <- y == 1
    zeros <- y == 0
    non_boundary <- !(ones | zeros)
    x <- x[c(which(ones), which(zeros), rep(which(non_boundary), 2)), , drop = FALSE]
    y <- c(y[ones], y[zeros], rep(c(0., 1.), each = sum(non_boundary)))
    # Run linear program
    out <- lp(x = x, y = y,
              linear_program = control$linear_program,
              purpose = control$purpose,
              tolerance = control$tolerance,
              solver = control$solver,
              solver_control = control$solver_control)
    if (is.na(out$outcome)) {
        if (identical(control$implementation, "ROI")) {
            warning("unexpected result from implementation ", control$implementation, " with solver: ", control$solver, "\n")
        }
        if (identical(control$implementation, "lsSolveAPI") & !log_link) {
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
        betas_all[betas_names] <- betas[betas_names]
    }
    out <- list(x = x,
                y = y,
                coefficients = betas_all,
                outcome = out$outcome,
                control = control)
    return(out)
}
