## test defaults
out0 <- detect_separation_control()

expect_equal(out0$linear_program, "primal")
expect_equal(out0$solver, "lpsolve")
expect_equal(out0$purpose, "find")
expect_equal(out0$tolerance, 1e-04)
expect_equal(out0$separator, detectseparation:::separator_ROI)
expect_equal(out0$implementation, "ROI")

## test arguments
out1 <- detect_separation_control(implementation = "lpSolveAPI", solver = "alabama", linear_program = "dual", solver_control = list(abs = 2), tolerance = 1e-07)

expect_equal(out1$linear_program, "dual")
expect_equal(out1$solver, "alabama")
expect_equal(out1$purpose, "find")
expect_equal(out1$tolerance, 1e-07)
expect_equal(out1$separator, detectseparation:::separator_lpSolveAPI)
expect_equal(out1$implementation, "lpSolveAPI")

