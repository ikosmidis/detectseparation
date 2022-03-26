## test defaults
out0 <- detect_separation_control()

expect_equal(out0$linear_program, "primal")
expect_equal(out0$solver, "lpsolve")
expect_equal(out0$purpose, "find")
expect_equal(out0$tolerance, 1e-04)
expect_equal(out0$separator, detectseparation::separator_ROI)
expect_equal(out0$implementation, "ROI")
