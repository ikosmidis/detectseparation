cone <- detectseparation:::binomial_cone_constraints
roi <- detectseparation:::coordinate_ranges_ROI
api <- detectseparation:::coordinate_ranges_lpSolveAPI
engine <- detectseparation:::.coordinate_ranges

capture_warnings <- function(code) {
    captured <- new.env(parent = emptyenv())
    captured$messages <- character()
    value <- withCallingHandlers(code, warning = function(w) {
        captured$messages <- c(captured$messages, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = captured$messages)
}

expected_ranges <- function(lower, upper, coefficient_names) {
    matrix(c(lower, upper), ncol = 2L,
           dimnames = list(coefficient_names, c("lower", "upper")))
}

## Hand-derived ranges cover finite, either forced sign, and ambiguous limits.
X_two <- cbind("(Intercept)" = 1, x = c(-1, 1))
X_missed <- cbind("(Intercept)" = 1,
                  x1 = c(0, 0, 1, 0, -1),
                  x2 = c(0, 0, 0, 1, -1),
                  x3 = c(0, 0, 0, 0, 1))
X_overlap <- cbind("(Intercept)" = 1, x = c(-1, 0, 1, 2))
X_log <- cbind("(Intercept)" = 1, x = c(0, 1))
X_flat <- cbind("(Intercept)" = 1, x = c(-1, 0, 1))

cases <- list(
    ## d_b >= |d_a|: intercept has both signs, slope only a positive sign.
    ambiguous = list(A = cone(X_two, c(0, 1)),
                     expected = expected_ranges(c(-1, 0), c(1, 1), colnames(X_two))),
    ## d_a = 0, d_1,d_2 >= 0, d_3 >= d_1 + d_2: all three slopes diverge.
    missed_infinities = list(A = cone(X_missed, c(0, 1, 1, 1, 1)),
                             expected = expected_ranges(rep(0, 4), c(0, 1, 1, 1),
                                                        colnames(X_missed))),
    negative_infinities = list(A = cone(X_missed, c(1, 0, 0, 0, 0)),
                               expected = expected_ranges(c(0, -1, -1, -1), rep(0, 4),
                                                          colnames(X_missed))),
    overlap = list(A = cone(X_overlap, c(0, 1, 0, 1)),
                   expected = expected_ranges(c(0, 0), c(0, 0), colnames(X_overlap))),
    ## Log link: the success at x = 0 fixes d_a; the failure allows d_b <= 0.
    log_infinite = list(A = cone(X_log, c(1, 0), log_link = TRUE),
                        expected = expected_ranges(c(0, -1), c(0, 0), colnames(X_log))),
    log_mixed = list(A = cone(X_log, c(0.5, 0), log_link = TRUE),
                     expected = expected_ranges(c(0, -1), c(0, 0), colnames(X_log))),
    ## All failures: d_a <= -|d_b|, so the slope has both signs.
    log_ambiguous = list(A = cone(X_two, c(0, 0), log_link = TRUE),
                         expected = expected_ranges(c(-1, -1), c(0, 1), colnames(X_two))),
    ## Zero cone despite nonunique finite log-binomial maximizers.
    log_finite = list(A = cone(X_flat, c(1, 0, 1), log_link = TRUE),
                      expected = expected_ranges(c(0, 0), c(0, 0), colnames(X_flat))),
    ## d_2 >= 2*d_1 >= 0: the box restricts d_1 to [0, 1/2].
    fractional_range = list(A = rbind(c(a = 1, b = 0), c(a = -2, b = 1)),
                            expected = expected_ranges(c(0, 0), c(0.5, 1), c("a", "b"))),
    one_column = list(A = matrix(c(1, -1), ncol = 1L, dimnames = list(NULL, "a")),
                      expected = expected_ranges(0, 0, "a")),
    no_constraints = list(A = matrix(numeric(), 0L, 2L,
                                     dimnames = list(NULL, c("a", "b"))),
                          expected = expected_ranges(c(-1, -1), c(1, 1), c("a", "b"))),
    no_coordinates = list(A = matrix(numeric(), 2L, 0L),
                          expected = expected_ranges(numeric(), numeric(), NULL)))

## Required backends always run. Optional ROI plugins are checked separately.
backends <- list(ROI_lpsolve = function(A) roi(A), lpSolveAPI = api)
for (solver in c("glpk", "ecos")) {
    if (requireNamespace(paste0("ROI.plugin.", solver), quietly = TRUE)) {
        backends[[paste0("ROI_", solver)]] <- local({
            selected_solver <- solver
            function(A) roi(A, solver = selected_solver)
        })
    }
}
for (backend in names(backends)) {
    for (case_name in names(cases)) {
        case <- cases[[case_name]]
        out <- backends[[backend]](case$A)
        info <- paste(backend, case_name)
        expect_equal(out$ranges, case$expected, tolerance = 1e-04, info = info)
        expect_equal(out$status, matrix(0L, nrow(case$expected), 2L,
                                        dimnames = dimnames(case$expected)), info = info)
    }
}

## Mixed responses contribute both signs, including with one-column designs.
X_mixed <- matrix(c(1, 2), ncol = 1L)
expect_equal(cone(X_mixed, c(0.5, 0)), matrix(c(1, -1, -2), ncol = 1L))
expect_equal(cone(X_mixed, c(0.5, 0), log_link = TRUE),
             matrix(c(-1, -2, 1, -1), ncol = 1L))

## A failed coordinate solve must remain unresolved, without erasing other bounds.
partial_result <- capture_warnings(engine(diag(2), function(objective, maximum) {
    if (!maximum) return(list(status = 2L))
    list(status = 0L, solution = objective, objective = 1)
}, tolerance = 1e-04))
partial <- partial_result$value
expect_equal(partial_result$warnings,
             paste0("Coordinate LP calculations failed for: ",
                    "coefficient 1 (minimum: solver status 2); ",
                    "coefficient 2 (minimum: solver status 2). ",
                    "The affected coefficient classifications are unresolved."))
expect_true(all(is.na(partial$ranges[, "lower"])))
expect_equal(unname(partial$ranges[, "upper"]), c(1, 1))
expect_equal(unname(partial$status), matrix(c(2L, 2L, 0L, 0L), 2L, 2L))

failed_result <- capture_warnings(engine(diag(2), function(...) list(status = 3L),
                                        tolerance = 1e-04))
failed <- failed_result$value
expect_equal(length(failed_result$warnings), 1L)
expect_true(grepl("maximum: solver status 3", failed_result$warnings, fixed = TRUE))
expect_true(all(is.na(failed$ranges)))
expect_true(all(failed$status == 3L))

## Optimality status alone does not make an unusable primal solution trustworthy.
invalid_results <- list(
    missing_solution = list(status = 0L, solution = c(NA_real_, 0), objective = NA_real_),
    outside_box = list(status = 0L, solution = c(2, 2), objective = 2),
    violated_cone = list(status = 0L, solution = c(-1, -1), objective = -1),
    inconsistent_objective = list(status = 0L, solution = c(0, 0), objective = 0.5),
    wrong_dimension = list(status = 0L, solution = 0, objective = 0))
for (reason in names(invalid_results)) {
    invalid_result <- capture_warnings(engine(diag(2), function(...) invalid_results[[reason]],
                                             tolerance = 1e-04))
    invalid <- invalid_result$value
    expect_equal(length(invalid_result$warnings), 1L, info = reason)
    expect_true(grepl("solution failed numerical checks", invalid_result$warnings,
                     fixed = TRUE), info = reason)
    expect_true(all(is.na(invalid$ranges)), info = reason)
    expect_true(all(invalid$status == 0L), info = reason)
}

## A small equality residual is rejected when it exceeds the chosen tolerance.
equality <- matrix(c(1, -1), ncol = 1L)
residual <- function(...) list(status = 0L, solution = 1e-10, objective = 1e-10)
loose <- engine(equality, residual, tolerance = 1e-04)
tight_result <- capture_warnings(engine(equality, residual, tolerance = 1e-12))
tight <- tight_result$value
expect_equal(length(tight_result$warnings), 1L)
expect_true(grepl("solution failed numerical checks", tight_result$warnings, fixed = TRUE))
expect_equal(unname(loose$ranges), matrix(1e-10, 1L, 2L), tolerance = 0)
expect_true(all(is.na(tight$ranges)))
expect_true(all(tight$status == 0L))

## Different failure reasons and bounds are collected into one named warning.
named_A <- diag(3)
colnames(named_A) <- c("a", "b", "c")
mixed <- capture_warnings(engine(named_A, function(objective, maximum) {
    if (objective[1L] == 1 && !maximum) return(list(status = 2L))
    if (objective[3L] == 1 && !maximum) return(list(status = NA_integer_))
    if (objective[2L] == 1 && maximum) {
        return(list(status = 0L, solution = 2 * objective, objective = 2))
    }
    d <- if (maximum) objective else numeric(3)
    list(status = 0L, solution = d, objective = as.numeric(maximum))
}, tolerance = 1e-04))
expect_equal(mixed$warnings,
             paste0("Coordinate LP calculations failed for: ",
                    "a (minimum: solver status 2); c (minimum: solver status NA); ",
                    "b (maximum: solution failed numerical checks). ",
                    "The affected coefficient classifications are unresolved."))
expect_equal(unname(mixed$value$ranges), matrix(c(NA, 0, NA, 1, NA, 1), 3L, 2L))
expect_equal(unname(mixed$value$status), matrix(c(2L, 0L, NA, 0L, 0L, 0L), 3L, 2L))

## Successful and empty calculations issue no warnings.
for (backend in backends) {
    expect_silent(backend(cases$ambiguous$A))
    expect_silent(backend(cases$no_constraints$A))
    expect_silent(backend(cases$no_coordinates$A))
}

## Strict input checks in both backends.
for (bad in list(matrix(NA_real_, 1L, 1L), matrix(Inf, 1L, 1L), matrix("x", 1L, 1L))) {
    expect_error(roi(bad), "finite numeric matrix")
    expect_error(api(bad), "finite numeric matrix")
}
for (bad in list(NA_real_, -1, Inf, c(1, 2))) {
    expect_error(roi(diag(2), tolerance = bad), "finite nonnegative number")
    expect_error(api(diag(2), tolerance = bad), "finite nonnegative number")
}
expect_error(cone(X_two, 1), "initialized binomial response")
expect_error(cone(X_two, c(-1, 1)), "initialized binomial response")

## Solver controls must not override whether each coordinate LP minimizes or maximizes.
controlled <- api(cases$ambiguous$A, solver_control = list(sense = "max", presolve = "rows"))
expect_equal(controlled$ranges, cases$ambiguous$expected)
expect_true(all(controlled$status == 0L))
