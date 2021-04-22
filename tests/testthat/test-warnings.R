context("detect_separation warnings")

## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "detectseparation")

expect_warning(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = binomial("log"),
                                             method = "detect_separation"), regexp = "reliable")


expect_warning(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = binomial("identity"),
                                             method = "detect_separation"), regexp = "reliable")


expect_warning(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = poisson("log"),
                                             method = "detect_separation"), regexp = "binomial-response")
