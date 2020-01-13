context("separator and separator_konis results agree")

data("endometrial", package = "detectseparation")
endo_sep_konis <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                      family = binomial("cloglog"),
                      method = "detect_separation",
                      separator = "separator_konis")

endo_sep_new <- update(endo_sep_konis, separator = "separator")

test_that("separator returns the same result as separator_konis [1]", {
    expect_identical(endo_sep_konis$separation, endo_sep_new$separation)
    expect_equal(coef(endo_sep_konis), coef(endo_sep_new))    
})


if (requireNamespace("AER", quietly = TRUE)) {
    data("MurderRates", package = "AER")
    murder_formula <- I(executions > 0) ~ time + income + noncauc + lfp + southern
    murder_sep_konis <- glm(murder_formula, data = MurderRates,
                            family = binomial(),
                            method = "detect_separation", separator = "separator_konis")
    murder_sep_new <- update(murder_sep_konis, separator = "separator")

    test_that("separator returns the same result as separator_konis [2]", {
        expect_identical(murder_sep_konis$separation, murder_sep_new$separation)
        expect_equal(coef(murder_sep_konis), coef(murder_sep_new))    
    })    
}
