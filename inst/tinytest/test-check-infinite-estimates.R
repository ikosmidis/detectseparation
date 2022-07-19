data("endometrial", package = "detectseparation")

endometrial_ml <- glm(HG ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
inf_check <- check_infinite_estimates(endometrial_ml)

## NV is infinite
expect_true(max(diff(inf_check[, "NV"])) > 1e+06)

## Make sure a plot method is available for inf_check objects
expect_true("plot.inf_check" %in% methods(class = "inf_check"))
