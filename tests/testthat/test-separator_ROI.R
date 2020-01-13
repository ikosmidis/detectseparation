context("ROI implementation")


if (requireNamespace("AER", quietly = TRUE)) {
    data("MurderRates", package = "AER")
    murder_formula <- I(executions > 0) ~ time + income + noncauc + lfp + southern
    murder_sep_lpsolve <- glm(murder_formula, data = MurderRates,
                              family = binomial(),
                              method = "detect_separation",
                              implementation = "ROI",
                              solver = "lpsolve")

    murder_sep_glpk <- update(murder_sep_lpsolve,
                              implementation = "ROI",
                              solver = "glpk")

    
    test_that("ROI implementation returns the same result with solver lpsolve and solve glpk", {
        expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_glpk))    
    })    

}

