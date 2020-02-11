context("detect_separation output")

## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "detectseparation")
endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                              family = binomial("cloglog"),
                              method = "detect_separation")

## The lizards example from ?brglm::brglm
data("lizards", package = "detectseparation")
lizards_separation <- glm(cbind(grahami, opalinus) ~ height + diameter +
                              light + time, family = binomial(logit), data = lizards,
                          method = "detect_separation")


test_that("infinte estimates have been found as expected", {
    expect_equal(coef(endometrial_separation), c(0L, -Inf, 0L, 0L), check.attributes = FALSE)
    expect_equal(coef(lizards_separation), rep(0L, 6), check.attributes = FALSE)
})


endometrial_separation_lpsolve <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                      family = binomial("cloglog"),
                                      method = "detect_separation",
                                      implementation = "lpSolveAPI")

endometrial_separation_lpsolve2 <- update(endometrial_separation_lpsolve,
                                          linear_program = "dual",
                                          purpose = "test")


test_that("output is as expected", {
    expect_output(print(endometrial_separation), "ROI \\| Solver: lpsolve")
    expect_output(print(endometrial_separation), "0: finite value, Inf: infinity, -Inf: -infinity")
    expect_output(print(endometrial_separation), "Separation: TRUE")
    expect_output(print(lizards_separation), "ROI \\| Solver: lpsolve")
    expect_output(print(lizards_separation), "0: finite value, Inf: infinity, -Inf: -infinity")
    expect_output(print(endometrial_separation_lpsolve), "Implementation: lpSolveAPI \\| Linear program: primal \\| Purpose: find")
    expect_output(print(endometrial_separation_lpsolve), "Separation: TRUE")
    expect_output(print(endometrial_separation_lpsolve2), "Implementation: lpSolveAPI \\| Linear program: dual \\| Purpose: test \\nSeparation: TRUE")
})


## ## hepatitis
## data("hepatitis", package = "pmlr")
## hepat <- hepatitis
## hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
## hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
## y <- rowSums(hepat$counts*nnet::class.ind(hepat$type)[,c(1, 3)])
## glm(y/counts ~ group * time, data = hepat, weights = counts, family = binomial(),
##     method = "detect_separation")
## dd <- data.frame(y = c(1,1,1,0,0), x = c(1,2,4,4,5), off = c(1,2, 2,4,3))
## summary(glm(y ~ x + offset(off), family = binomial(), data = dd))
