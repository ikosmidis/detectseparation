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
#' Method for [glm()] that detects infinite components in
#' the maximum likelihood estimates of generalized linear models
#' with binomial responses.
#'
#' @inheritParams stats::glm.fit
#'
#' @aliases detectInfiniteEstimates
#'
#' @param x `x` is a design matrix of dimension `n * p`.
#' @param y `y` is a vector of observations of length `n`.
#' @param control a list of parameters controlling separation
#'     detection. See [detect_separation_control()] for
#'     details.
#' @param start currently not used.
#' @param mustart currently not used.
#' @param etastart currently not used.
#' @param singular.ok logical. If `FALSE`, a singular model is an
#'     error.
#'
#'
#' @details
#'
#' For binomial-response generalized linear models with `"log"` link,
#' separated data allocations do not necessarily lead to infinite
#' maximum likelihood estimates. For this reason, for models with the
#' `"log"` link [detect_infinite_estimates()] relies on an alternative
#' linear optimization model developed in Schwendinger et
#' al. (2021). For all the other supported links
#' [detect_infinite_estimates()] relies on the linear programming
#' methods developed in Konis (2007). Log-binomial calculations use ROI;
#' non-log models use the selected implementation. See [detect_separation()]
#' for definitions, the coordinate-wise LP calculation, coefficient flag
#' meanings, numerical failure handling, and the log-binomial feasibility
#' assumption.
#'
#' [detect_infinite_estimates()] can be passed directly as
#' a method to the [glm()] function. See, examples.
#'
#' [detectInfiniteEstimates()] is an alias for
#' [detect_infinite_estimates()].
#'
#' @return A list of class `"detect_infinite_estimates"`, also inheriting from
#'     `"glm"` and `"lm"` when used as a [glm()] method. The common output
#'     components are described in the Value section of [detect_separation()].
#'     Here `outcome` is `TRUE` if infinite components are detected, `FALSE`
#'     if finiteness is established, and `NA` if the decision remains unresolved.
#'
#' @author Ioannis Kosmidis \[aut, cre\] \email{ioannis.kosmidis@warwick.ac.uk}, Florian Schwendinger \[aut\] \email{FlorianSchwendinger@gmx.at}, Dirk Schumacher \[aut\] \email{mail@dirk-schumacher.net}, Kjell Konis \[ctb\] \email{kjell.konis@me.com}
#'
#' @seealso [glm()], [detect_separation()], [check_infinite_estimates()], [brglm2::brglm_fit()]
#'
#' @references
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' <https://www.jstor.org/stable/2984941>
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
#' Schwendinger, F., Grün, B. & Hornik, K. (2021). A comparison of
#' optimization solvers for log binomial regression including conic
#' programming.  *Computational Statistics*, **36**,
#' 1721–1754. \doi{10.1007/s00180-021-01084-5}
#'
#'
#' @examples
#' # The classical example given in Silvapulle (1981) can be utilized
#' # to show that for the log-binomial model there exist data allocations
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
#' # The maximum likelihood fit is
#' glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"), start = c(-1, 0))
#'
#' @export
detect_infinite_estimates <- function(x, y, weights = NULL,
                                      start = NULL, etastart = NULL,  mustart = NULL,
                                      offset = NULL, family = gaussian(),
                                      control = list(), intercept = TRUE, singular.ok = TRUE) {
    if (isTRUE(family$family != "binomial")) {
        warning("`detect_infinite_estimates()` has been developed for use with binomial-response GLMs")
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

        coefs <- coefficients(x)
        display <- format(coefs, digits = digits, trim = TRUE)

        labels <- c(aliased = "NA(A)",
                    undetermined = "NA(UD)",
                    unresolved = "NA(UR)")
        missing <- is.na(coefs)
        display[missing] <- unname(labels[x$coefficient_status[missing]])

        print(display, quote = FALSE, right = TRUE)
        cat("0: finite value, Inf: infinity, -Inf: -infinity, ",
            "NA(A): aliased, NA(UD): undetermined, NA(UR): unresolved\n",
            sep = "")
        ## cat("Existence of maximum likelihood estimates\n")
        ## print(coefficients(x))
        ## cat("0: finite value, Inf: infinity, -Inf: -infinity, NA: undetermined, aliased or unresolved\n")
    }
}

## Workhorse function
.detect_infinite_estimates <- function(x, y, weights = NULL,
                                       start = NULL, etastart = NULL,  mustart = NULL,
                                       offset = NULL, family = gaussian(),
                                       control = list(), intercept = TRUE, singular.ok = TRUE,
                                       log_link = FALSE, classify = TRUE) {
    control <- do.call("detect_separation_control", control)
    find_coefficients <- classify && (log_link || control$implementation == "ROI" ||
                                      control$purpose == "find")
    lp <- if (log_link) {
        dielb_ROI
    } else {
        getNamespace("detectseparation")[[paste("separator", control$implementation, sep = "_")]]
    }
    # ensure x is a matrix
    x <- as.matrix(x)
    betas_names_all <- if (is.null(colnames(x))) make.names(seq_len(NCOL(x))) else colnames(x)
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
    # Initialization can change weights, for example for zero binomial totals.
    # Only observations with positive effective weight determine the rank.
    keep <- weights > 0
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    # Detect aliasing
    estimable <- seq_len(nvars)
    qrx <- qr(x)
    rank <- qrx$rank
    is_full_rank <- rank == nvars
    if (!singular.ok && !is_full_rank) {
        stop("singular fit encountered")
    }
    if (!isTRUE(is_full_rank)) {
        aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
        x <- x[, -aliased, drop = FALSE]
        estimable <- estimable[-aliased]
    }
    betas_all <- structure(rep(NA_real_, length(betas_names_all)), .Names = betas_names_all)
    coefficient_status <- structure(rep("aliased", nvars), .Names = betas_names_all)
    coefficient_status[estimable] <- "unresolved"
    if (control$solver == "alabama" & is.null(control$solver_control$start)) {
        control$solver_control$start <- rep(0, rank)
    }
    # Reshape data set: keep 0 and 1, and replace anything in (0,
    # 1) with one zero and one 1
    ones <- y == 1
    zeros <- y == 0
    non_boundary <- !(ones | zeros)
    x <- x[c(which(ones), which(zeros), rep(which(non_boundary), 2)), , drop = FALSE]
    y <- c(y[ones], y[zeros], rep(c(0., 1.), each = sum(non_boundary)))
    # Run linear program
    if (rank == 0) {
        # No nonzero fitted direction remains after removing aliased columns.
        out <- list(outcome = FALSE, beta = numeric(0), status = 0L)
    } else {
        out <- lp(x = x, y = y,
                  linear_program = control$linear_program,
                  purpose = control$purpose,
                  tolerance = control$tolerance,
                  solver = control$solver,
                  solver_control = control$solver_control)
    }
    A <- binomial_cone_constraints(x, y, log_link = log_link)
    if (!log_link && !is.null(out$beta) && rank > 0L && !is.na(out$outcome)) {
        # Use the strict-margin objective for the global decision, and reject
        # an infeasible returned direction rather than interpreting its signs.
        d <- out$beta
        if (length(d) != rank || any(!is.finite(d)) ||
            any(abs(d) > 1 + control$tolerance) ||
            any(drop(A %*% d) < -control$tolerance * pmax(1, rowSums(abs(A))))) {
            out$outcome <- NA
        } else {
            out$outcome <- sum(A %*% d) > control$tolerance
        }
    }
    coordinate_lp <- NULL
    if (!find_coefficients) {
        betas_all <- NULL
        coefficient_status <- NULL
    } else if (identical(out$outcome, FALSE)) {
        # A successful global finiteness check needs no coordinate LPs.
        betas_all[estimable] <- 0
        coefficient_status[estimable] <- "finite"
    } else {
        coordinate_lp <- if (log_link || control$implementation == "ROI") {
            coordinate_ranges_ROI(A, solver = control$solver,
                                  tolerance = control$tolerance,
                                  solver_control = control$solver_control)
        } else {
            coordinate_ranges_lpSolveAPI(A, tolerance = control$tolerance)
        }
        ranges <- coordinate_lp$ranges
        valid <- is.finite(ranges) & !is.na(coordinate_lp$status) & coordinate_lp$status == 0L
        resolved <- rowSums(valid) == 2L
        positive <- ranges[, "upper"] > control$tolerance
        negative <- ranges[, "lower"] < -control$tolerance
        finite <- which(resolved & !positive & !negative)
        pos_inf <- which(resolved & positive & !negative)
        neg_inf <- which(resolved & negative & !positive)
        undetermined <- which(resolved & positive & negative)
        betas_all[estimable[finite]] <- 0
        betas_all[estimable[pos_inf]] <- Inf
        betas_all[estimable[neg_inf]] <- -Inf
        coefficient_status[estimable[finite]] <- "finite"
        coefficient_status[estimable[pos_inf]] <- "positive_infinity"
        coefficient_status[estimable[neg_inf]] <- "negative_infinity"
        coefficient_status[estimable[undetermined]] <- "undetermined"
        if (is.na(out$outcome)) {
            # Any successfully established nonzero range certifies an improving
            # direction. A negative global decision needs every bound resolved.
            if (any(valid[, "upper"] & positive, na.rm = TRUE) ||
                any(valid[, "lower"] & negative, na.rm = TRUE)) {
                out$outcome <- TRUE
            } else if (all(resolved)) {
                out$outcome <- FALSE
            }
        }
    }
    out <- list(x = x,
                y = y,
                coefficients = betas_all,
                coefficient_status = coefficient_status,
                coordinate_lp = coordinate_lp,
                outcome = out$outcome,
                lp_status = out$status,
                control = control)
    return(out)
}
