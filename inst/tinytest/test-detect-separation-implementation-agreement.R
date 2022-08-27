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

    murder_sep_alabama <- update(murder_sep_lpsolve,
                                 implementation = "ROI",
                                 solver = "alabama")

    murder_sep_ecos <- update(murder_sep_lpsolve,
                                 implementation = "ROI",
                                 solver = "ecos")
    
    murder_sep_alabama1 <- update(murder_sep_lpsolve,
                                  implementation = "ROI",
                                  solver = "alabama",
                                  solver_control = list(start = rep(1, 6)))

    murder_sep_lpSolveAPI <- update(murder_sep_lpsolve,
                                    implementation = "lpSolveAPI")


    expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_glpk))
    expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_alabama))
    expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_ecos))
    expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_alabama1))
    expect_equal(coef(murder_sep_lpsolve), coef(murder_sep_lpSolveAPI))

}


data("endometrial", package = "detectseparation")
endo_sep_konis <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                      family = binomial("cloglog"),
                      method = "detect_separation",
                      implementation = "lpSolveAPI")

endo_sep_new <- update(endo_sep_konis, implementation = "ROI")


expect_identical(endo_sep_konis$separation, endo_sep_new$separation)
expect_equal(coef(endo_sep_konis), coef(endo_sep_new))    
