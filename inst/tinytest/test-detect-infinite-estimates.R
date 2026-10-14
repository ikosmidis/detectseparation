data("silvapulle1981", package = "detectseparation")

m_inf_log <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
                 method = "detect_infinite_estimates")
m_inf_logit <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("logit"),
                   method = "detect_infinite_estimates")

expect_equal(unname(coef(m_inf_log)), c(0, 0))
expect_equal(unname(coef(m_inf_logit)), c(-Inf, Inf))




## Empty model
mempty <- update(m_inf_log, y ~ -1)
expect_null(coef(mempty))
expect_false(mempty$outcome)

## Coordinate classifications through both public methods and implementations.
missed <- data.frame(y = c(0, 1, 1, 1, 1),
                     x1 = c(0, 0, 1, 0, -1),
                     x2 = c(0, 0, 0, 1, -1),
                     x3 = c(0, 0, 0, 0, 1))
two <- data.frame(y = c(0, 1), x = c(-1, 1))
overlap <- data.frame(y = c(0, 1, 0, 1), x = c(-1, 0, 1, 2))
controls <- list(list(implementation = "ROI"),
                 list(implementation = "lpSolveAPI", linear_program = "primal"),
                 list(implementation = "lpSolveAPI", linear_program = "dual"))
for (method in c("detect_infinite_estimates", "detect_separation")) {
    for (control in controls) {
        fit <- glm(y ~ ., data = missed, family = binomial(),
                   method = method, control = control)
        expect_equal(unname(coef(fit)), c(0, Inf, Inf, Inf))
        expect_equal(unname(fit$coefficient_status),
                     c("finite", rep("positive_infinity", 3)))
        expect_equal(names(fit$coefficient_status), names(coef(fit)))
        expect_true(fit$outcome)
        expect_true(all(fit$coordinate_lp$status == 0L))

        negative <- missed
        negative$y <- 1 - negative$y
        fit_negative <- glm(y ~ ., data = negative, family = binomial(),
                            method = method, control = control)
        expect_equal(unname(coef(fit_negative)), c(0, -Inf, -Inf, -Inf))
        expect_equal(unname(fit_negative$coefficient_status),
                     c("finite", rep("negative_infinity", 3)))

        fit_two <- glm(y ~ x, data = two, family = binomial(),
                       method = method, control = control)
        expect_equal(coef(fit_two), c("(Intercept)" = NA_real_, x = Inf))
        expect_equal(unname(fit_two$coefficient_status), c("undetermined", "positive_infinity"))
        expect_true(fit_two$outcome)

        fit_overlap <- glm(y ~ x, data = overlap, family = binomial(),
                           method = method, control = control)
        expect_equal(unname(coef(fit_overlap)), c(0, 0))
        expect_equal(unname(fit_overlap$coefficient_status), c("finite", "finite"))
        expect_false(fit_overlap$outcome)
        expect_null(fit_overlap$coordinate_lp)

        fit_alias <- glm(y ~ x1 + I(2 * x1) + x2 + x3, data = missed,
                         family = binomial(), method = method, control = control)
        expect_equal(unname(coef(fit_alias)), c(0, Inf, NA_real_, Inf, Inf))
        expect_equal(unname(fit_alias$coefficient_status),
                     c("finite", "positive_infinity", "aliased",
                       "positive_infinity", "positive_infinity"))
    }
}

## Same directional classification for the other supported full-support links.
for (link in c("probit", "cloglog", "cauchit")) {
    fit_link <- glm(y ~ x, data = two, family = binomial(link),
                    method = "detect_infinite_estimates")
    expect_equal(coef(fit_link), c("(Intercept)" = NA_real_, x = Inf))
    expect_equal(unname(fit_link$coefficient_status), c("undetermined", "positive_infinity"))
}

## Effective weights, grouped responses, and offsets preserve the cone signs.
mixed <- data.frame(success = c(1, 1, 1, 1), failure = c(1, 0, 0, 0),
                    x1 = c(0, 1, 0, -1), x2 = c(0, 0, 1, -1), x3 = c(0, 0, 0, 1),
                    off = c(2, -3, 4, -5))
fit_grouped <- glm(cbind(success, failure) ~ x1 + x2 + x3,
                   data = mixed, family = binomial(), method = "detect_infinite_estimates")
fit_fractional <- glm(success / (success + failure) ~ x1 + x2 + x3,
                      weights = success + failure, data = mixed, family = binomial(),
                      method = "detect_infinite_estimates")
expect_equal(unname(coef(fit_grouped)), c(0, Inf, Inf, Inf))
expect_equal(coef(fit_fractional), coef(fit_grouped))
expect_equal(fit_fractional$coefficient_status, fit_grouped$coefficient_status)
fit_offset <- detect_infinite_estimates(x = model.matrix(~ x1 + x2 + x3, mixed),
                                        y = mixed$success / (mixed$success + mixed$failure),
                                        weights = mixed$success + mixed$failure, offset = mixed$off,
                                        family = binomial())
expect_equal(coef(fit_offset), coef(fit_grouped))
expect_equal(fit_offset$coefficient_status, fit_grouped$coefficient_status)

## Distinguish the log-binomial likelihood cone from geometric separation.
log_data <- data.frame(y = c(1, 0), x = c(0, 1))
log_inf <- glm(y ~ x, data = log_data, family = binomial("log"),
                method = "detect_infinite_estimates")
expect_message(log_sep <- update(log_inf, method = "detect_separation"), "not necessarily")
expect_equal(coef(log_inf), c("(Intercept)" = 0, x = -Inf))
expect_equal(unname(log_inf$coefficient_status), c("finite", "negative_infinity"))
expect_true(log_inf$outcome)
expect_equal(coef(log_sep), coef(log_inf))
expect_equal(log_sep$coefficient_status, log_inf$coefficient_status)
expect_equal(log_sep$coordinate_lp, log_inf$coordinate_lp)

log_flat <- data.frame(y = c(1, 0, 1), x = c(-1, 0, 1))
flat <- glm(y ~ x, data = log_flat, family = binomial("log"),
             method = "detect_infinite_estimates")
expect_equal(unname(coef(flat)), c(0, 0))
expect_equal(unname(flat$coefficient_status), c("finite", "finite"))
expect_false(flat$outcome)
expect_null(flat$coordinate_lp)

log_all_failure <- glm(y ~ x, data = transform(two, y = 0), family = binomial("log"),
                        method = "detect_infinite_estimates")
expect_equal(unname(coef(log_all_failure)), c(-Inf, NA_real_))
expect_equal(unname(log_all_failure$coefficient_status), c("negative_infinity", "undetermined"))

## Aliasing and empty models keep their separate interpretation.
weighted <- data.frame(y = c(0, 1, 1), x = c(0, 0, 1), w = c(1, 1, 0))
weighted_fit <- glm(y ~ x, weights = w, data = weighted, family = binomial(),
                     method = "detect_infinite_estimates")
expect_equal(unname(weighted_fit$coefficient_status), c("finite", "aliased"))
zero_fit <- glm(y ~ x, weights = rep(0, 3), data = weighted, family = binomial(),
                 method = "detect_infinite_estimates")
expect_equal(unname(zero_fit$coefficient_status), c("aliased", "aliased"))
expect_false(zero_fit$outcome)
expect_null(mempty$coefficient_status)

## Optional solvers agree on classifications despite selecting different directions.
for (solver in c("glpk", "ecos")) {
    if (requireNamespace(paste0("ROI.plugin.", solver), quietly = TRUE)) {
        fit_solver <- glm(y ~ ., data = missed, family = binomial(),
                          method = "detect_infinite_estimates", solver = solver)
        expect_equal(unname(coef(fit_solver)), c(0, Inf, Inf, Inf))
        expect_equal(unname(fit_solver$coefficient_status), c("finite", rep("positive_infinity", 3)))
        log_solver <- update(log_inf, solver = solver)
        expect_equal(coef(log_solver), coef(log_inf))
        expect_equal(log_solver$coefficient_status, log_inf$coefficient_status)
        expect_true(log_solver$outcome)
    }
}

## Restore namespace bindings even if an expectation or a fit raises an error.
with_mock_lp <- function(name, value, code, package = "detectseparation") {
    original <- get(name, asNamespace(package))
    utils::assignInNamespace(name, value, ns = package)
    on.exit(utils::assignInNamespace(name, original, ns = package))
    force(code)
}

## No coordinate solves after overlap, or in the lpSolveAPI test-only mode.
unexpected <- function(...) stop("unexpected coordinate LP")
with_mock_lp("coordinate_ranges_ROI", unexpected, {
    fit_shortcut <- glm(y ~ x, data = overlap, family = binomial(),
                        method = "detect_infinite_estimates")
    expect_false(fit_shortcut$outcome)
    expect_equal(unname(coef(fit_shortcut)), c(0, 0))
})
with_mock_lp("coordinate_ranges_lpSolveAPI", unexpected, {
    for (linear_program in c("primal", "dual")) {
        fit_test <- glm(y ~ x, data = two, family = binomial(),
                       method = "detect_infinite_estimates", implementation = "lpSolveAPI",
                       linear_program = linear_program, purpose = "test")
        expect_true(fit_test$outcome)
        expect_null(coef(fit_test))
        expect_null(fit_test$coefficient_status)
        expect_null(fit_test$coordinate_lp)
    }
})

## A coordinate warning leaves successful classifications and separation intact.
partial_lp <- function(A, ...) {
    detectseparation:::.coordinate_ranges(A, function(objective, maximum) {
        if (objective[1L] == 1 && !maximum) return(list(status = 2L))
        d <- numeric(ncol(A))
        if (objective[2L] == 1 && maximum) d[2L] <- 1
        list(status = 0L, solution = d, objective = sum(objective * d))
    }, tolerance = 1e-04)
}
with_mock_lp("coordinate_ranges_ROI", partial_lp, {
    for (method in c("detect_infinite_estimates", "detect_separation")) {
        expect_warning(partial_fit <- glm(y ~ x, data = two, family = binomial(),
                                         method = method),
                       pattern = "minimum: solver status 2")
        expect_equal(unname(coef(partial_fit)), c(NA_real_, Inf))
        expect_equal(unname(partial_fit$coefficient_status), c("unresolved", "positive_infinity"))
        expect_true(partial_fit$outcome)
        expect_equal(partial_fit$coordinate_lp$status[1, "lower"], 2L)
    }
})

## A failed global LP can be resolved by coordinate evidence, not by its status.
failed_global <- function(...) list(outcome = NA, beta = NULL, status = 1L)
with_mock_lp("separator_ROI", failed_global, {
    recovered <- glm(y ~ x, data = two, family = binomial(), method = "detect_infinite_estimates")
    expect_true(recovered$outcome)
    expect_equal(unname(recovered$coefficient_status), c("undetermined", "positive_infinity"))
    recovered_overlap <- glm(y ~ x, data = overlap, family = binomial(),
                             method = "detect_infinite_estimates")
    expect_false(recovered_overlap$outcome)
    expect_equal(unname(coef(recovered_overlap)), c(0, 0))
})
failed_coordinates <- function(A, ...) {
    detectseparation:::.coordinate_ranges(A, function(...) list(status = 2L),
                                          tolerance = 1e-04)
}
with_mock_lp("separator_ROI", failed_global, {
    with_mock_lp("coordinate_ranges_ROI", failed_coordinates, {
        expect_warning(failed_fit <- glm(y ~ x, data = two, family = binomial(),
                                        method = "detect_infinite_estimates"),
                       pattern = "affected coefficient classifications are unresolved")
        expect_true(is.na(failed_fit$outcome))
        expect_true(all(is.na(coef(failed_fit))))
        expect_equal(unname(failed_fit$coefficient_status), c("unresolved", "unresolved"))
    })
})
unusable_coordinates <- function(A, ...) {
    detectseparation:::.coordinate_ranges(A, function(...) {
        list(status = 0L, solution = rep(NA_real_, ncol(A)), objective = NA_real_)
    }, tolerance = 1e-04)
}
with_mock_lp("coordinate_ranges_ROI", unusable_coordinates, {
    expect_warning(unusable_fit <- glm(y ~ x, data = two, family = binomial(),
                                      method = "detect_infinite_estimates"),
                   pattern = "solution failed numerical checks")
    expect_true(unusable_fit$outcome)
    expect_true(all(is.na(coef(unusable_fit))))
    expect_equal(unname(unusable_fit$coefficient_status), c("unresolved", "unresolved"))
})

## A nonoptimal log-binomial status other than unboundedness is unresolved.
op <- ROI::OP(objective = c(0, 1), constraints = ROI::L_constraint(diag(2), c(">=", ">="), c(0, 0)),
               bounds = ROI::V_bound(ld = -Inf, ud = Inf, nobj = 2L), maximum = TRUE)
timeout_result <- ROI::ROI_solve(op, solver = "lpsolve")
timeout_result$status$msg$symbol <- "LPS_TIMEOUT"
with_mock_lp("ROI_solve", function(...) timeout_result, {
    raw_log <- detectseparation:::dielb_ROI(model.matrix(y ~ x, log_data), log_data$y)
    expect_true(is.na(raw_log$outcome))
}, package = "ROI")
