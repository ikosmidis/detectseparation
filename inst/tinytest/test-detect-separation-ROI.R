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

## Check that the coefficients are as expected
expect_equal(coef(endometrial_separation), c(0L, -Inf, 0L, 0L), check.attributes = FALSE)
expect_equal(coef(lizards_separation), rep(0L, 6), check.attributes = FALSE)

## endometrial using lpSolveAPI
endometrial_separation_lpsolve <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                      family = binomial("cloglog"),
                                      method = "detect_separation",
                                      implementation = "lpSolveAPI")

## testing lpSolveAPI implementation options
endometrial_separation_lpsolve2 <- update(endometrial_separation_lpsolve,
                                          linear_program = "dual",
                                          purpose = "test")

## test output
expect_stdout(print(endometrial_separation), "ROI \\| Solver: lpsolve")
expect_stdout(print(endometrial_separation), "0: finite value, Inf: infinity, -Inf: -infinity")
expect_stdout(print(endometrial_separation), "Separation: TRUE")
expect_stdout(print(lizards_separation), "ROI \\| Solver: lpsolve")
expect_stdout(print(lizards_separation), "0: finite value, Inf: infinity, -Inf: -infinity")
expect_stdout(print(endometrial_separation_lpsolve), "Implementation: lpSolveAPI \\| Linear program: primal \\| Purpose: find")
expect_stdout(print(endometrial_separation_lpsolve), "Separation: TRUE")
expect_stdout(print(endometrial_separation_lpsolve2), "Implementation: lpSolveAPI \\| Linear program: dual \\| Purpose: test \\nSeparation: TRUE")

## Test separation type detection
counterexample_x <- rbind(c(1, -1, 2),
                          c(1, -1, -2),
                          c(1, -2, -2),
                          c(1, 2, 2))
counterexample_y <- c(0, 1, 1, 1)

counterexample_separation <- detect_separation(x = counterexample_x,
                                               y = counterexample_y,
                                               family = binomial(),
                                               control = list(separation_type = TRUE))
expect_true(counterexample_separation$outcome)
expect_true(counterexample_separation$complete)
expect_stdout(print(counterexample_separation), "Separation: TRUE \\(complete\\)")

counterexample_separation_lpsolve <- detect_separation(x = counterexample_x,
                                                       y = counterexample_y,
                                                       family = binomial(),
                                                       control = list(implementation = "lpSolveAPI",
                                                                      linear_program = "primal",
                                                                      purpose = "test",
                                                                      separation_type = TRUE))
expect_true(counterexample_separation_lpsolve$outcome)
expect_true(counterexample_separation_lpsolve$complete)

overlap_x <- cbind(1, c(-1, 0, 1, 2))
overlap_y <- c(0, 1, 0, 1)
overlap_separation <- detect_separation(x = overlap_x,
                                        y = overlap_y,
                                        family = binomial(),
                                        control = list(separation_type = TRUE))
expect_false(overlap_separation$outcome)
expect_null(overlap_separation$complete)

## Test aliasing
endometrial_separation_a <- glm(HG ~ NV + I(2 * NV) + PI + EH, data = endometrial,
                                family = binomial("cloglog"),
                                method = "detect_separation")
expect_equal(coef(endometrial_separation_a), c(0L, Inf, NA, 0L, 0L), check.attributes = FALSE)


## Test EMPTY fits
endometrial_separation_n <- glm(HG ~ 0, data = endometrial,
                                family = binomial("cloglog"),
                                method = "detect_separation")
expect_null(coef(endometrial_separation_n))
expect_false(endometrial_separation_n$outcome)


## Test log-binomial models
data("silvapulle1981", package = "detectseparation")

expect_message(m1sep <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
                            method = "detect_separation"), pattern = "necessarily result in")

m1inf <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
             method = "detect_infinite_estimates")

expect_message(m1sep_type <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
                                 method = "detect_separation", separation_type = TRUE),
               pattern = "necessarily result in")

expect_equal(unname(coef(m1sep)), unname(coef(m1inf)))
expect_equal(unname(coef(m1sep)), c(0, 0))
expect_true(m1sep_type$outcome)
expect_false(m1sep_type$complete)





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
