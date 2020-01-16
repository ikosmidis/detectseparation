#  Copyright (C) 2020 Dirk Schumacher, Ioannis Kosmidis
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

separator_ROI <- function(x, y,
                          solver = "lpsolve",
                          tolerance = 1e-03,
                          ...) {
  betas_names <- dimnames(x)[[2L]]

  ## IK: 13/01/2020: 
  ## detect_separation ensures that y is always 1 or zero so below checks is not necessary
  # just to be safe
  ## unique_classes <- unique(y)
  ## stopifnot(
  ##   is.numeric(unique_classes),
  ##   length(unique_classes) == 2L,
  ##   all(unique_classes %in% c(0, 1))
  ## )

  # the model here is based on Konis (2007), chapter 4. In particular
  # sections 4.2, 4.4.3 and 4.4.4

  # build the ROI model

  # transform the model matrix so that all constraints are >=
  # that should also work with doubles
  y[y == 0] <- -1
  X_bar <- x * y

  m <- ncol(X_bar)
  n <- nrow(X_bar)

  constraints <- ROI::L_constraint(X_bar, rep.int(">=", n), rep.int(0, n))

  bounds <- ROI::V_bound(
    li = seq_len(m), lb = rep.int(-1, m),
    ui = seq_len(m), ub = rep.int(1, m)
  )

  # max t(rep.int(1, n)) %*% X_bar %*% beta = colSums(X_bar) %*% beta
  # subject to X_bar >= 0
  # beta between -1 and 1
  opt_model <- ROI::OP(
    objective = colSums(X_bar),
    constraints = constraints,
    types = rep.int("C", m),
    bounds = bounds,
    maximum = TRUE
  )

  # This will be now checked by match.arg in detect_separation_control
  # ensure the solver is loaded using the ROI plugin mechanism
  ## require_solver(solver)

  # if the LP is unbounded, seperation exists. Otherwise an optiomal solution
  # with obj. value 0 exists.
  result <- ROI::ROI_solve(opt_model, solver = solver)
  
  ## an optimal solution should always exists
  ## stopifnot(identical(as.integer(result$status$code), 0L))
  ## IK: 13/01/2020
  ## stoping with a message instead
  if (!isTRUE(identical(result$status$code, 0L))) 
      stop("unexpected result from ", solver)
  
  # compare to 0 zero with tolerance
  solution <- ROI::solution(result, "primal")
  non_zero <- abs(solution) > tolerance
  names(solution) <- betas_names
  has_seperation <- any(non_zero, na.rm = TRUE)

  list(
    separation = has_seperation,
    beta = solution
  )
}

## require_solver <- function(solver_name) {
##   solver_name <- match.arg(solver_name, choices = c("lpsolve", "glpk"))
##   roi_plugin_name <- paste0("ROI.plugin.", solver_name)
##   if (!requireNamespace(roi_plugin_name, quietly = TRUE)) {
##     stop(
##       "No ROI solver plugin loaded for linear programs. ",
##       "Please install the package ", roi_plugin_name,
##       " or use a different solver.",
##       call. = FALSE
##     )
##   }
## }
