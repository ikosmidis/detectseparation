separator <- function(x, y,
                      linear_program = c("primal", "dual"), purpose = c("test", "find"), # this will most probably not be necessary
                      tolerance = 1e-03,
                      solver = "lpsolve") {
  betas_names <- dimnames(x)[[2L]]

  # just to be safe
  unique_classes <- unique(y)
  stopifnot(
    is.numeric(unique_classes),
    length(unique_classes) == 2L,
    all(unique_classes %in% c(0, 1))
  )

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

  # ensure the solver is loaded using the ROI plugin mechanism
  require_solver(solver)

  # if the LP is unbounded, seperation exists. Otherwise an optiomal solution
  # with obj. value 0 exists.
  result <- ROI::ROI_solve(opt_model, solver = solver)

  # an optimal solution should always exists
  stopifnot(identical(as.integer(result$status$code), 0L))

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

require_solver <- function(solver_name) {
  solver_name <- match.arg(solver_name, choices = c("lpsolve", "glpk"))
  roi_plugin_name <- paste0("ROI.plugin.", solver_name)
  if (!requireNamespace(roi_plugin_name, quietly = TRUE)) {
    stop(
      "No ROI solver plugin loaded for linear programs. ",
      "Please install the package ", roi_plugin_name,
      " or use a different solver.",
      call. = FALSE
    )
  }
}
