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
    has_seperation <- any(abs(sol) > tolerance, na.rm = TRUE)
    # an optimal solution should always exist
    if (!isTRUE(ROI::solution(result, "status_code") == 0L)) {
        has_separation <- NA
    }
    list(outcome = has_seperation,
         beta = sol)
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
    has_infinite_estimates <-  !isTRUE(ROI::solution(result, "status_code") == 0L)
    list(outcome = has_infinite_estimates,
         beta = sol)
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
        if (status != 0) {
            ans$outcome <- NA
        }
        beta <- lpSolveAPI::get.variables(lp)
        if (any(abs(beta) > tolerance)) {
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
        beta <- lpSolveAPI::get.dual.solution(lp)[2:(p+1)]
        if (any(abs(beta) > tolerance)) {
            ans$outcome <- TRUE
        }
        else {
            ans$outcome <- FALSE
        }
        ans$beta <- beta
    }
    ans
}
