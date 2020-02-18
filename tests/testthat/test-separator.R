context("ROI and lpSolveAPI implementations agree")

data("endometrial", package = "detectseparation")
endo_sep_konis <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                      family = binomial("cloglog"),
                      method = "detect_separation",
                      implementation = "lpSolveAPI")

endo_sep_new <- update(endo_sep_konis, implementation = "ROI")

test_that("ROI implementation returns the same result as lpSolveAPI implementation [1]", {
    expect_identical(endo_sep_konis$separation, endo_sep_new$separation)
    expect_equal(coef(endo_sep_konis), coef(endo_sep_new))    
})


if (requireNamespace("AER", quietly = TRUE)) {
    data("MurderRates", package = "AER")
    murder_formula <- I(executions > 0) ~ time + income + noncauc + lfp + southern
    murder_sep_konis <- glm(murder_formula, data = MurderRates,
                            family = binomial(),
                            method = "detect_separation",
                            implementation = "lpSolveAPI")
    murder_sep_lpsolve <- update(murder_sep_konis,
                                 implementation = "ROI",
                                 solver = "lpsolve")
    
    test_that("ROI implementation returns the same result as lpSolveAPI implementation [2]", {
        expect_identical(murder_sep_konis$separation, murder_sep_lpsolve$separation)
        expect_equal(coef(murder_sep_konis), coef(murder_sep_lpsolve))    
    })

}

