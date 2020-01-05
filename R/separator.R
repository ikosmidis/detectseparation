separator <- function(x, y,
                      linear_program = c("primal", "dual"), purpose = c("test", "find"), # this will most probably not be necessary
                      tolerance = 1e-03) {
    n <- dim(x)[1L]
    p <- dim(x)[2L]

    ## Some dummy code with garbage output

    ans <- list()

    ans$separation <- TRUE
    ans$beta <- rep.int(0, p)

    ans
}
