## Rank must be assessed on observations with positive effective weights.
weighted <- data.frame(y = c(0, 1, 1), x = c(0, 0, 1), w = c(1, 1, 0))
grouped <- data.frame(success = c(0, 1, 0), failure = c(1, 0, 0),
                      x = c(0, 0, 1))
controls <- list(list(implementation = "ROI"),
                 list(implementation = "lpSolveAPI", linear_program = "primal"),
                 list(implementation = "lpSolveAPI", linear_program = "dual"))

for (method in c("detect_separation", "detect_infinite_estimates")) {
    for (control in controls) {
        fit <- glm(y ~ x, data = weighted, weights = w, family = binomial(),
                   method = method, control = control)
        expect_equal(coef(fit), c("(Intercept)" = 0, x = NA_real_))
        expect_false(fit$outcome)
        expect_equal(dim(fit$x), c(2L, 1L))
        expect_error(glm(y ~ x, data = weighted, weights = w,
                         family = binomial(), method = method,
                         control = control, singular.ok = FALSE),
                     "singular fit encountered")

        ## Initialization changes effective weights for zero trial totals.
        fit_grouped <- glm(cbind(success, failure) ~ x, data = grouped,
                           family = binomial(), method = method, control = control)
        expect_equal(coef(fit_grouped), coef(fit))
        expect_false(fit_grouped$outcome)
        expect_error(glm(cbind(success, failure) ~ x, data = grouped,
                         family = binomial(), method = method,
                         control = control, singular.ok = FALSE),
                     "singular fit encountered")

        ## One retained column must still be passed to the LP as a matrix.
        fit_alias <- glm(y ~ x, data = weighted[1:2, ], family = binomial(),
                         method = method, control = control)
        expect_equal(coef(fit_alias), coef(fit))
        expect_false(fit_alias$outcome)

        ## An aliased column must not hide separation in the retained column.
        separated <- weighted
        separated$y <- 1
        fit_separated <- glm(y ~ x, data = separated, weights = w,
                             family = binomial(), method = method, control = control)
        expect_equal(coef(fit_separated), c("(Intercept)" = Inf, x = NA_real_))
        expect_true(fit_separated$outcome)

        ## No estimable columns: all coefficients are aliased, with no separator.
        fit_zero <- glm(y ~ 0 + x, data = weighted, weights = w,
                        family = binomial(), method = method, control = control)
        expect_equal(coef(fit_zero), c(x = NA_real_))
        expect_false(fit_zero$outcome)
        expect_equal(dim(fit_zero$x), c(2L, 0L))

        fit_no_weight <- glm(y ~ x, data = weighted, weights = rep(0, 3),
                             family = binomial(), method = method, control = control)
        expect_equal(coef(fit_no_weight), c("(Intercept)" = NA_real_, x = NA_real_))
        expect_false(fit_no_weight$outcome)
        expect_equal(dim(fit_no_weight$x), c(0L, 0L))
    }
}

## Preserve the coefficient-free lpSolveAPI test mode when rank is zero.
for (linear_program in c("primal", "dual")) {
    fit_test <- glm(y ~ 0 + x, data = weighted, weights = w,
                    family = binomial(), method = "detect_separation",
                    implementation = "lpSolveAPI", linear_program = linear_program,
                    purpose = "test")
    expect_null(coef(fit_test))
    expect_false(fit_test$outcome)
}

## The log-binomial LP uses the same effective design and alias handling.
fit_log <- glm(cbind(success, failure) ~ x, data = grouped,
               family = binomial("log"), method = "detect_infinite_estimates")
expect_equal(coef(fit_log), c("(Intercept)" = 0, x = NA_real_))
expect_false(fit_log$outcome)
expect_error(glm(cbind(success, failure) ~ x, data = grouped,
                 family = binomial("log"), method = "detect_infinite_estimates",
                 singular.ok = FALSE), "singular fit encountered")
