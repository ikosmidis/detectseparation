# detect infinite estimates log binomial (dielb)
dielb_ROI <- function(x, y,
                      solver = "lpsolve",
                      tolerance = 1e-03,
                      solver_control = list(),
                      ...) {
  ## Use linear program from
  ##
  ## Schwendinger, F., Grün, B. & Hornik, K. A comparison of
  ## optimization solvers for log binomial regression including
  ## conic programming. Comput Stat 36, 1721–1754
  ## (2021). https://doi.org/10.1007/s00180-021-01084-5

  direction <- ifelse(y == 0, leq(1), eq(1))
  op <- OP(objective = -colSums(x), maximum = TRUE,
           constraints = L_constraint(x, direction, double(nrow(x))),
           bounds = V_bound(ld = -Inf, ud = Inf, nobj = ncol(x)))

  if (isTRUE(solver == "lpsolve")) {
      control <- list(pivoting = "firstindex", simplextype = c("primal", "primal"))
      solver_control <- modifyList(control, solver_control)
  }

  s <- ROI_solve(op, solver, control = solver_control)
  sol <- ROI::solution(s)
  names(sol) <- colnames(x)
  has_infinite_estimates <- !isTRUE(solution(s, "status_code") == 0L)
  list(infinite_estimates = has_infinite_estimates, beta = sol)
}


detect_infinite_estimates_log_binomial <- function(x, y, weights = rep(1, nobs),
                                                   start = NULL, etastart = NULL,  mustart = NULL,
                                                   offset = rep.int(0, nobs), family = gaussian(),
                                                   control = list(), intercept = TRUE, singular.ok = TRUE) {
    control <- do.call("detect_infinite_estimates_control", control)
    separator <- dielb_ROI
    ## ensure x is a matrix
    x <- as.matrix(x)
    betas_names_all <- betas_names <- if (is.null(colnames(x))) make.names(seq_len(NCOL(x))) else colnames(x)
    ##
    nobs <- NROW(y)
    nvars <- ncol(x)
    if (nvars == 0) {
        return(list(separation = FALSE, control = control))
    }
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missingOffset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    ## Initialize as prescribed in family
    eval(family$initialize)
    if (control$solver == "alabama" & is.null(control$solver_control$start)) {
        control$solver_control$start <- rep(0, nvars)
    }
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
        betas_names <- betas_names[-aliased]
    }

    betas_all <- structure(rep(NA_real_, length(betas_names_all)), .Names = betas_names_all)
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
    if (is.na(out$infinite_estimates)) {
        warning("unexpected result from implementation ", control$implementation, " with solver: ", control$solver, "\n")
    }
    if (is.null(out$beta)) {
        betas_all <- NULL
    } else {
        betas <- out$beta
        names(betas) <- betas_names
        inds <- abs(betas) < control$tolerance
        betas <- Inf * betas
        betas[inds] <- 0
        # Why betas_all is indexed by the manes
        betas_all[betas_names] <- betas[betas_names]
    }
    out <- list(x = x,
                y = y,
                coefficients = betas_all,
                infinite_estimates = out$infinite_estimates,
                control = control)
    return(out)
}



#' @title Detect Infinite Estimates
#'
#' @description
#' Method for \code{\link{glm}} that detects infinite components in
#' the MLE estimates of generalized linear models with binomial responses.
#'
#' For all links except the \code{"log"} link, function
#' \code{detect_infinite_estimates} is a wrapper around
#' \code{\link{detect_separation}}. For the \code{"log"}
#' link separated data allocations not necessarily lead to
#' infinite components in the MLE estimates, therefore another method is used.
#'
#' \code{\link{detect_infinite_estimates}} for the \code{"log"} link
#' relies on the linear optimization model developed in
#' Schwendinger et al. (2021) and for all the other
#' supported links it relies on the linear programming methods
#' developed in Konis (2007).
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
#' @references
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' \url{http://www.jstor.org/stable/2984941}
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
#' Schwendinger, F., Grün, B. & Hornik, K. A comparison of optimization solvers
#' for log binomial regression including conic programming.
#' Comput Stat 36, 1721–1754 (2021). \url{https://doi.org/10.1007/s00180-021-01084-5}
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
detect_infinite_estimates <- function(x, y, weights = rep(1, nobs),
                                      start = NULL, etastart = NULL,  mustart = NULL,
                                      offset = rep.int(0, nobs), family = gaussian(),
                                      control = list(), intercept = TRUE, singular.ok = TRUE) {
    if (isTRUE(family$family != "binomial")) {
        stop("`detect_infinite_estimates` has been developed for use with binomial-response GLMs")
    }
    if (isTRUE(family$link == "log")) {
        .detect_infinite_estimates <- detect_infinite_estimates_log_binomial
    } else {
        .detect_infinite_estimates <- detect_separation
    }
    out <- .detect_infinite_estimates(x = x, y = y, weights = weights, start = start,
                                      etastart = etastart,  mustart = mustart,
                                      offset = offset, family = family, control = control,
                                      intercept = control, singular.ok = singular.ok)
    names(out)[names(out) == "separation"] <- "infinite_estimates"
    class(out) <- out$class <- "detect_infinite_estimates"
    out
}


#' @method print detect_infinite_estimates
#' @export
print.detect_infinite_estimates <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
    cat("Implementation:", x$control$implementation, "| ")
    if (identical(x$control$implementation, "ROI")) {
        cat("Solver:", x$control$solver, "\n")
    } else {
        cat("Linear program:", x$control$linear_program, "| Purpose:", x$control$purpose, "\n")
    }
    cat("Infinite estimates:", x$infinite_estimates, "\n")
    if (!is.null(x$coefficients)) {
        cat("Existence of maximum likelihood estimates\n")
        print(coefficients(x))
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}


#' Auxiliary function for the \code{\link{glm}} interface when
#' \code{method} is \code{\link{detect_infinite_estimates}}.
#'
#' Typically only used internally by \code{\link{detect_infinite_estimates}}
#' but may be used to construct a \code{control} argument.
#'
#' @aliases detectInfiniteEstimatesControl
#' @param solver should the linear program be solved using the
#'     \code{"lpsolve"} (using the \pkg{ROI.plugin.lpsolve} package;
#'     default) or another solver? Alternative solvers are
#'     \code{"glpk"}, \code{"cbc"}, \code{"clp"}, \code{"cplex"},
#'     \code{"ecos"}, \code{"gurobi"}, \code{"scs"},
#'     \code{"symphony"}. If \pkg{ROI.plugin.[solver]} is not
#'     installed then the user will be prompted to install it before
#'     continuing.
#' @param tolerance maximum absolute variable value from the linear
#'     program, before separation is declared. Default is
#'     \code{1e-04}.
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
detect_infinite_estimates_control <- function(solver = "lpsolve",
                                              tolerance = 1e-04,
                                              solver_control = list()) {
    ## ensure the solver is loaded using the ROI plugin mechanism
    if (solver != "lpsolve") {
        if (!solver %in% names(ROI_registered_solvers())) {
            plugin_name <- sprintf("ROI.plugin.%s", gsub("\\..*", "", solver))
            if (plugin_name %in% ROI_installed_solvers()) {
                requireNamespace(plugin_name, quietly = TRUE)
            } else {
                stop(sprintf("'%s' can not be found among the installed solvers ", plugin_name),
                     "(in `ROI_installed_solvers()`) please make sure that is installed.")
            }
        }
    }
    list(solver = solver,
         solver_control = solver_control,
         tolerance = tolerance,
         implementation = "ROI")
}

check_solver <- function(solver) {

}
