settings <- expand.grid(lp = c("dual", "primal"),
                        purpose = c("find", "test"),
                        stringsAsFactors = FALSE)


if (requireNamespace("AER", quietly = TRUE)) {

    data("MurderRates", package = "AER")
    murder_formula <- I(executions > 0) ~ time + income + noncauc + lfp + southern
    murder_sep <- as.list(numeric(4))
    for (j in seq.int(nrow(settings))) {
        murder_sep[[j]] <- glm(murder_formula, data = MurderRates,
                               family = binomial(),
                               method = "detect_separation",
                               implementation = "lpSolveAPI",
                               linear_program = settings[j, "lp"],
                               purpose = settings[j, "purpose"])
    }
    names(murder_sep) <- apply(settings, 1, paste, collapse = "-")
    ## test that lpSolveAPI implementation returns the same result across linear_program and purpose settings [MurderRates]"
    expect_equal(coef(murder_sep[["dual-find"]]), coef(murder_sep[["primal-find"]]))
    expect_equal(murder_sep[["dual-find"]]$separation, murder_sep[["primal-find"]]$separation)
    expect_equal(murder_sep[["dual-find"]]$separation, murder_sep[["primal-test"]]$separation)
    expect_equal(murder_sep[["dual-find"]]$separation, murder_sep[["dual-test"]]$separation)
}

data("endometrial", package = "detectseparation")
endo_sep <- as.list(numeric(4))
for (j in seq.int(nrow(settings))) {
    endo_sep[[j]] <- glm(HG ~ I(-NV) + PI + EH,
                         data = endometrial,
                         family = binomial("cloglog"),
                         method = "detect_separation",
                         implementation = "lpSolveAPI",
                         linear_program = settings[j, "lp"],
                         purpose = settings[j, "purpose"])
}
names(endo_sep) <- apply(settings, 1, paste, collapse = "-")
## test that lpSolveAPI implementation returns the same result across linear_program and purpose settings [endometrial]
expect_equal(coef(endo_sep[["dual-find"]]), coef(endo_sep[["primal-find"]]))
expect_equal(endo_sep[["dual-find"]]$separation, endo_sep[["primal-find"]]$separation)
expect_equal(endo_sep[["dual-find"]]$separation, endo_sep[["primal-test"]]$separation)
expect_equal(endo_sep[["dual-find"]]$separation, endo_sep[["dual-test"]]$separation)

data("lizards", package = "detectseparation")
lizo_sep <- as.list(numeric(4))
for (j in seq.int(nrow(settings))) {
    lizo_sep[[j]] <- glm(cbind(grahami, opalinus) ~ height + diameter +
                             light + time,
                         data = lizards,
                         family = binomial("logit"),
                         method = "detect_separation",
                         implementation = "lpSolveAPI",
                         linear_program = settings[j, "lp"],
                         purpose = settings[j, "purpose"])
}
names(lizo_sep) <- apply(settings, 1, paste, collapse = "-")
## test that lpSolveAPI implementation returns the same result across linear_program and purpose settings [lizards]
expect_equal(coef(lizo_sep[["dual-find"]]), coef(lizo_sep[["primal-find"]]))
expect_equal(lizo_sep[["dual-find"]]$separation, lizo_sep[["primal-find"]]$separation)
expect_equal(lizo_sep[["dual-find"]]$separation, lizo_sep[["primal-test"]]$separation)
expect_equal(lizo_sep[["dual-find"]]$separation, lizo_sep[["dual-test"]]$separation)

