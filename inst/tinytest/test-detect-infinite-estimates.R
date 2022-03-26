data("silvapulle1981", package = "detectseparation")

m_inf_log <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("log"),
                 method = "detect_infinite_estimates")
m_inf_logit <- glm(y ~ ghqs, data = silvapulle1981, family = binomial("logit"),
                   method = "detect_infinite_estimates")

expect_equal(unname(coef(m_inf_log)), c(0, 0))
expect_equal(unname(coef(m_inf_logit)), c(-Inf, Inf))




## Empty model
mempty <- update(m_inf_log, y ~ -1)
expect_null(coef(mempty))
expect_false(mempty$outcome)
