#  separation_type_ROI:
#  Copyright (C) 2026- Ioannis Kosmidis
#
#  separator_ROI:
#  Copyright (C) 2021- Dirk Schumacher, Ioannis Kosmidis
#
#  separator_lpSolveAPI: Port of the separator function from the safeBinaryRegression (version 0.1-3) R package (see safeBinaryRegression/R/separator.R)
#  Copyright (C) 2009-2013 Kjell Konis; Copyright (C) 2017- Ioannis Kosmidis
#
#  dielb_ROI:
#  Copyright (C) 2022- Florian Schwendinger, Ioannis Kosmidis
#
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

separator_ROI <- function(x, y,
                          solver = "lpsolve",
                          tolerance = 1e-03,
                          solver_control = list(),
                          ...) {
    # The model here is based on Konis (2007), chapter 4. In particular
    # sections 4.2, 4.4.3 and 4.4.4
    #
    # build the ROI model
    # transform the model matrix so that all constraints are >=
    # that should also work with doubles
    y[y == 0] <- -1
    X_bar <- x * y
    m <- ncol(X_bar)
    n <- nrow(X_bar)
    constraints <- ROI::L_constraint(X_bar, rep.int(">=", n), rep.int(0, n))
    bounds <- ROI::V_bound(li = seq_len(m),
                           lb = rep.int(-1, m),
                           ui = seq_len(m),
                           ub = rep.int(1, m))
    # max t(rep.int(1, n)) %*% X_bar %*% beta = colSums(X_bar) %*% beta
    # subject to X_bar >= 0
    # beta between -1 and 1
    opt_model <- ROI::OP(objective = colSums(X_bar),
                         constraints = constraints,
                         types = rep.int("C", m),
                         bounds = bounds,
                         maximum = TRUE)

    # if the LP is unbounded, seperation exists. Otherwise an optiomal solution
    # with obj. value 0 exists.
    result <- ROI::ROI_solve(opt_model, solver = solver, control = solver_control)
    # compare to 0 zero with tolerance
    sol <- ROI::solution(result, "primal")
    names(sol) <- colnames(x)
    has_separation <- any(abs(sol) > tolerance, na.rm = TRUE)
    # an optimal solution should always exist
    if (!isTRUE(ROI::solution(result, "status_code") == 0L)) {
        has_separation <- NA
    }
    list(outcome = has_separation,
         beta = sol,
         status = ROI::solution(result, "status_code"))
}

# detect infinite estimates log binomial (dielb)
dielb_ROI <- function(x, y,
                      solver = "lpsolve",
                      tolerance = 1e-03,
                      solver_control = list(),
                      ...) {
    # Use linear program from
    #
    # Schwendinger, F., Grün, B. & Hornik, K. A comparison of
    # optimization solvers for log binomial regression including
    # conic programming. Comput Stat 36, 1721–1754
    # (2021). https://doi.org/10.1007/s00180-021-01084-5

    direction <- ifelse(y == 0, ROI::leq(1), ROI::eq(1))
    op <- ROI::OP(objective = -colSums(x), maximum = TRUE,
                  constraints = ROI::L_constraint(x, direction, double(nrow(x))),
                  bounds = ROI::V_bound(ld = -Inf, ud = Inf, nobj = ncol(x)))

    if (isTRUE(solver == "lpsolve")) {
        control <- list(pivoting = "firstindex", simplextype = c("primal", "primal"))
        solver_control <- utils::modifyList(control, solver_control)
    }

    result <- ROI::ROI_solve(op, solver, control = solver_control)
    sol <- ROI::solution(result)
    names(sol) <- colnames(x)
    status <- ROI::solution(result, "status_code")
    # A nonoptimal status need not mean unboundedness. Other statuses leave
    # the global decision unresolved; coordinate LPs can still supply evidence.
    symbol <- as.character(result$status$msg$symbol)
    unbounded <- length(symbol) == 1L &&
        symbol %in% c("LPS_UNBOUNDED", "GLP_UNBND", "ECOS_DINF")
    has_infinite_estimates <- if (isTRUE(status == 0L)) FALSE else if (unbounded) TRUE else NA
    if (identical(has_infinite_estimates, FALSE) &&
        (length(sol) != ncol(x) || any(!is.finite(sol)) || any(abs(sol) > tolerance))) {
        has_infinite_estimates <- NA
    }
    list(outcome = has_infinite_estimates,
         beta = sol,
         status = status)
}

#  Port of the separator function from the safeBinaryRegression
#  (version 0.1-3) R package (see safeBinaryRegression/R/separator.R)
separator_lpSolveAPI <- function(x, y,
                                 linear_program = c("primal", "dual"),
                                 purpose = c("test", "find"),
                                 tolerance = 1e-03,
                                 ...) {
    n <- dim(x)[1L]
    p <- dim(x)[2L]
    p_seq <- seq.int(p)
    zeros <- rep.int(0, n)
    dimnames(x) <- NULL
    y.bar <- -sign(y - 0.5)
    x.bar <- y.bar * x
    ans <- list()
    linear_program <- match.arg(linear_program)
    purpose <- match.arg(purpose)
    if (linear_program == "primal" && purpose == "test") {
        lp <- lpSolveAPI::make.lp(n, p)
        for(j in p_seq) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
        }
        status <- lpSolveAPI::set.rhs(lp,  zeros)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(1L, n))
        status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.bounds(lp, lower = rep(-Inf, p), upper = rep(Inf, p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                                          simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        if (status == 0) {
            ans$outcome <- FALSE
        }
        else {
            if (status == 3) {
                ans$outcome <- TRUE
            }
            else {
                ans$outcome <- NA
            }
        }
    }
    if (linear_program == "primal" && purpose == "find") {
        lp <- lpSolveAPI::make.lp(n, p)
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
        }
        status <- lpSolveAPI::set.rhs(lp, zeros)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(1, n))
        status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.bounds(lp, lower = rep.int(-1, p), upper = rep.int(1, p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                              simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        beta <- if (status == 0) lpSolveAPI::get.variables(lp) else NULL
        if (status != 0) {
            ans$outcome <- NA
        }
        else if (any(abs(beta) > tolerance)) {
            ans$outcome <- TRUE
        }
        else {
            ans$outcome <- FALSE
        }
        ans$beta <- beta
    }
    if (linear_program == "dual" && purpose == "test") {
        lp <- lpSolveAPI::make.lp(p, n)
        for (j in 1:n) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
        }
        status <- lpSolveAPI::set.rhs(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.constr.type(lp, rep.int(3, p))
        status <- lpSolveAPI::set.objfn(lp, zeros)
        status <- lpSolveAPI::set.bounds(lp, lower = zeros, upper = rep(Inf, n))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                                          simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        if (status == 0) {
            ans$outcome <- FALSE
        }
        else {
            if (status == 2) {
                ans$outcome <- TRUE
            }
            else {
                ans$outcome <- NA
            }
        }
    }
    if (linear_program == "dual" && purpose == "find") {
        lp <- lpSolveAPI::make.lp(p, n + 2*p)
        for (j in 1:n) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
        }
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, n+j, -1.0, j)
        }
        # IK, 12 April 2017: p_seq instead 1:n below;
        # safeBinaryRegression:::separator (version 0.1-3) has 1:n
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, n+p+j, 1.0, j)
        }
        b <- -colSums(x.bar)
        status <- lpSolveAPI::set.rhs(lp, b)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(3, p))
        status <- lpSolveAPI::set.objfn(lp, rep.int(c(0.0, 1.0), c(n, 2*p)))
        status <- lpSolveAPI::set.bounds(lp, lower = rep.int(0, n + 2*p), upper = rep(Inf, n + 2*p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                                          simplextype = c("primal", "primal"))
        basis <- p_seq
        basis[b >= 0.0] <- basis[b >= 0.0] + p
        status <- lpSolveAPI::set.basis(lp, -(n + p + basis))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        beta <- if (status == 0) lpSolveAPI::get.dual.solution(lp)[2:(p+1)] else NULL
        if (status != 0) {
            ans$outcome <- NA
        }
        else if (any(abs(beta) > tolerance)) {
            ans$outcome <- TRUE
        }
        else {
            ans$outcome <- FALSE
        }
        ans$beta <- beta
    }
    ans$status <- status
    ans
}

# Cone constraints for initialized responses, after removing zero weights.
# Mixed responses contribute both signs of a row. For the log link, success
# predictors must remain constant and failure predictors must not increase.
binomial_cone_constraints <- function(x, y, log_link = FALSE) {
    x <- as.matrix(x)
    if (length(y) != nrow(x) || any(!is.finite(y)) || any(y < 0 | y > 1)) {
        stop("y must contain one initialized binomial response per row of x")
    }
    success <- x[y > 0, , drop = FALSE]
    failure <- x[y < 1, , drop = FALSE]
    if (log_link) {
        rbind(-failure, success, -success)
    } else {
        rbind(success, -failure)
    }
}

# Directional minima and maxima over {d: A d >= 0, -1 <= d <= 1}.
# These are ranges of directions, not intervals for coefficient estimates.
# Each backend supplies an optimality status, a primal solution, and its
# objective. An NA range denotes a nonoptimal or numerically unusable solve;
# status codes retain their backend-specific meanings (zero means optimal).
.coordinate_ranges <- function(A, solve, tolerance, make_solver = NULL) {
    A <- as.matrix(A)
    if (!is.numeric(A) || any(!is.finite(A))) {
        stop("A must be a finite numeric matrix")
    }
    if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
        stop("tolerance must be a finite nonnegative number")
    }
    q <- ncol(A)
    labels <- list(colnames(A), c("lower", "upper"))
    ranges <- matrix(NA_real_, q, 2L, dimnames = labels)
    status <- matrix(NA_integer_, q, 2L, dimnames = labels)
    if (q == 0L) {
        return(list(ranges = ranges, status = status))
    }
    if (nrow(A) == 0L) {
        ranges[, "lower"] <- -1
        ranges[, "upper"] <- 1
        status[] <- 0L
        return(list(ranges = ranges, status = status))
    }
    # Construct reusable solver objects only after validation and empty cases.
    if (!is.null(make_solver)) {
        solve <- make_solver(A)
    }
    # Check primal feasibility with row-relative residual tolerances.
    margin_tolerance <- tolerance * pmax(1, rowSums(abs(A)))
    for (j in seq_len(q)) {
        objective <- numeric(q)
        objective[j] <- 1
        for (k in seq_len(2L)) {
            result <- solve(objective, maximum = k == 2L)
            status[j, k] <- result$status
            if (!isTRUE(result$status == 0L)) {
                next
            }
            d <- result$solution
            value <- result$objective
            if (length(d) != q || any(!is.finite(d)) ||
                length(value) != 1L || !is.finite(value)) {
                next
            }
            if (any(abs(d) > 1 + tolerance) ||
                any(drop(A %*% d) < -margin_tolerance) ||
                abs(value - d[j]) > tolerance * max(1, abs(value), abs(d[j]))) {
                next
            }
            # Zero is feasible, so minima cannot be positive and maxima cannot
            # be negative, apart from the allowed numerical residual.
            if ((k == 1L && value > tolerance) ||
                (k == 2L && value < -tolerance)) {
                next
            }
            ranges[j, k] <- value
        }
    }
    failed <- which(is.na(ranges), arr.ind = TRUE)
    if (nrow(failed) > 0L) {
        coefficient_names <- colnames(A)
        if (is.null(coefficient_names)) {
            coefficient_names <- paste("coefficient", seq_len(q))
        }
        details <- vapply(seq_len(nrow(failed)), function(i) {
            j <- failed[i, 1L]
            k <- failed[i, 2L]
            reason <- if (isTRUE(status[j, k] == 0L)) {
                "solution failed numerical checks"
            } else {
                paste("solver status", status[j, k])
            }
            paste0(coefficient_names[j], " (", c("minimum", "maximum")[k],
                   ": ", reason, ")")
        }, character(1L))
        warning("Coordinate LP calculations failed for: ",
                paste(details, collapse = "; "),
                ". The affected coefficient classifications are unresolved.",
                call. = FALSE)
    }
    list(ranges = ranges, status = status)
}

coordinate_ranges_ROI <- function(A, solver = "lpsolve", tolerance = 1e-04,
                                  solver_control = list()) {
    make_solver <- function(A) {
        q <- ncol(A)
        constraints <- ROI::L_constraint(A, rep.int(">=", nrow(A)),
                                          numeric(nrow(A)))
        bounds <- ROI::V_bound(li = seq_len(q), lb = rep.int(-1, q),
                               ui = seq_len(q), ub = rep.int(1, q))
        function(objective, maximum) {
            op <- ROI::OP(objective = objective, constraints = constraints,
                           bounds = bounds, maximum = maximum)
            result <- ROI::ROI_solve(op, solver = solver, control = solver_control)
            status <- ROI::solution(result, "status_code")
            if (!isTRUE(status == 0L)) {
                return(list(status = status))
            }
            list(status = status, solution = ROI::solution(result, "primal"),
                 objective = ROI::solution(result, "objval"))
        }
    }
    .coordinate_ranges(A, solve = NULL, tolerance = tolerance, make_solver = make_solver)
}

coordinate_ranges_lpSolveAPI <- function(A, tolerance = 1e-04,
                                         solver_control = list()) {
    make_solver <- function(A) {
        q <- ncol(A)
        lp <- lpSolveAPI::make.lp(nrow(A), q)
        for (j in seq_len(q)) {
            lpSolveAPI::set.column(lp, j, A[, j])
        }
        lpSolveAPI::set.rhs(lp, numeric(nrow(A)))
        lpSolveAPI::set.constr.type(lp, rep.int(">=", nrow(A)))
        lpSolveAPI::set.bounds(lp, lower = rep.int(-1, q),
                               upper = rep.int(1, q))
        control <- utils::modifyList(list(pivoting = "firstindex",
                                         simplextype = c("primal", "primal")),
                                     solver_control)
        function(objective, maximum) {
            lpSolveAPI::set.objfn(lp, objective)
            control$sense <- if (maximum) "max" else "min"
            do.call(lpSolveAPI::lp.control, c(list(lprec = lp), control))
            # A basis from a different objective or sense need not be suitable.
            lpSolveAPI::set.basis(lp, default = TRUE)
            status <- lpSolveAPI::solve.lpExtPtr(lp)
            if (!isTRUE(status == 0L)) {
                return(list(status = status))
            }
            list(status = status, solution = lpSolveAPI::get.variables(lp),
                 objective = lpSolveAPI::get.objective(lp))
        }
    }
    .coordinate_ranges(A, solve = NULL, tolerance = tolerance, make_solver = make_solver)
}

separation_type_ROI <- function(x, y,
                         solver = "lpsolve",
                         tolerance = 1e-03,
                         solver_control = list(),
                         ...) {
    # The model here is based on the transformed design matrix X_bar and
    # solves a max-min-margin problem. In particular, it maximizes the
    # minimum margin t subject to X_bar %*% beta >= t. A strictly positive
    # optimal value for t implies complete separation.
    #
    # build the ROI model
    # transform the model matrix so that all constraints are >=
    # that should also work with doubles
    y[y == 0] <- -1
    X_bar <- x * y
    m <- ncol(X_bar)
    n <- nrow(X_bar)
    constraints <- ROI::L_constraint(cbind(X_bar, -1), rep.int(">=", n), rep.int(0, n))
    bounds <- ROI::V_bound(li = seq_len(m),
                           lb = rep.int(-1, m),
                           ui = seq_len(m),
                           ub = rep.int(1, m),
                           nobj = m + 1)
    # max t
    # subject to X_bar %*% beta >= t
    # beta between -1 and 1
    opt_model <- ROI::OP(objective = c(rep.int(0, m), 1),
                         constraints = constraints,
                         types = rep.int("C", m + 1),
                         bounds = bounds,
                         maximum = TRUE)
    result <- ROI::ROI_solve(opt_model, solver = solver, control = solver_control)
    status <- ROI::solution(result, "status_code")
    sol <- ROI::solution(result, "primal")
    if (length(sol) != m + 1L) {
        sol <- rep(NA_real_, m + 1L)
    }
    beta <- sol[seq_len(m)]
    names(beta) <- colnames(x)
    minimum_margin <- unname(sol[m + 1])
    has_complete_separation <- isTRUE(minimum_margin > tolerance)
    if (!isTRUE(status == 0L)) {
        has_complete_separation <- NA
    }
    list(outcome = has_complete_separation,
         beta = beta,
         minimum_margin = minimum_margin)
}
