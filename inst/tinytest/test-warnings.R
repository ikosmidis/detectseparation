## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "detectseparation")

expect_message(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = binomial("log"),
                                             method = "detect_separation"), pattern = "not necessarily")


expect_warning(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = binomial("identity"),
                                             method = "detect_separation"), pattern = "reliable")


expect_warning(endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                                             family = poisson("log"),
                                             method = "detect_separation"), pattern = "for use with binomial-response")
