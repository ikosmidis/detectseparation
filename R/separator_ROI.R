#  Copyright (C) 2020 Dirk Schumacher; Ioannis Kosmidis
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
                          solver_control = list(),
                          ...) {
  betas_names <- dimnames(x)[[2L]]
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
  bounds <- ROI::V_bound(li = seq_len(m),
                         lb = rep.int(-1, m),
                         ui = seq_len(m),
                         ub = rep.int(1, m))
  ## max t(rep.int(1, n)) %*% X_bar %*% beta = colSums(X_bar) %*% beta
  ## subject to X_bar >= 0
  ## beta between -1 and 1
  opt_model <- ROI::OP(objective = colSums(X_bar),
                       constraints = constraints,
                       types = rep.int("C", m),
                       bounds = bounds,
                       maximum = TRUE)
  
  ## if the LP is unbounded, seperation exists. Otherwise an optiomal solution
  ## with obj. value 0 exists.  
  result <- ROI::ROI_solve(opt_model, solver = solver, control = solver_control)  
  ## compare to 0 zero with tolerance
  solution <- ROI::solution(result, "primal")
  non_zero <- abs(solution) > tolerance
  names(solution) <- betas_names
  has_seperation <- any(non_zero, na.rm = TRUE) 
  ## an optimal solution should always exists
  if (!isTRUE(identical(result$status$code, 0L))) {
      has_separation <- NA
  }  
  list(separation = has_seperation,
       beta = solution)
}
