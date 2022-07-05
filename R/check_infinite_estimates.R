# Copyright (C) 2016- Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


#' A simple diagnostic of whether the maximum likelihood estimates are
#' infinite
#'
#'
#' @param object the result of a \code{\link{glm}} call.
#' @param nsteps starting from \code{maxit = 1}, the GLM is refitted
#'     for \code{maxit = 2}, \code{maxit = 3}, \ldots, \code{maxit =
#'     nsteps}. Default value is 30.
#' @param ... currently not used.
#'
#' @return
#'
#' An object of class \code{inf_check} that has a \code{plot} method.
#'
#' @details
#'
#' \code{check_infinite_estimates}() attempts to identify the occurrence
#' of infinite estimates in GLMs with binomial responses by
#' successively refitting the model. At each iteration the maximum
#' number of allowed IWLS iterations is fixed starting from 1 to
#' \code{nsteps} (by setting \code{control = glm.control(maxit = j)},
#' where \code{j} takes values 1, \ldots, nsteps in
#' \code{\link{glm}}). For each value of \code{maxit}, the estimated
#' asymptotic standard errors are divided to the corresponding ones
#' from \code{control = glm.control(maxit = 1)}. Then, based on the
#' results in Lesaffre & Albert (1989), if the sequence of ratios in
#' any column of the resultant matrix diverges, then complete or
#' quasi-complete separation occurs and the maximum likelihood
#' estimate for the corresponding parameter has value minus or plus
#' infinity.
#'
#' \code{check_infinite_estimates}() can also be used to identify the
#' occurrence of infinite estimates in baseline category logit models
#' for nominal responses (see \code{\link[brglm2]{brmultinom}()} from
#' the \pkg{brglm2} R package), and adjacent category logit models for
#' ordinal responses (see \code{\link[brglm2]{bracl}()} from the
#' \pkg{brglm2} R package).
#'
#' @return
#'
#' A matrix inheriting from class \code{inf_check}, with \code{nsteps}
#' rows and \code{p} columns, where \code{p} is the number of model
#' parameters. A \code{plot} method is provided for \code{inf_check}
#' objects for the easy inspection of the ratios of the standard
#' errors.
#'
#' @note
#'
#' For the definition of complete and quasi-complete separation, see
#' Albert and Anderson (1984). Kosmidis and Firth (2021) prove that
#' the reduced-bias estimator that results by the penalization of the
#' logistic regression log-likelihood by Jeffreys prior takes always
#' finite values, even when some of the maximum likelihood estimates
#' are infinite. The reduced-bias estimates can be computed using the
#' \pkg{brglm2} R package.
#'
#' @seealso \code{\link[nnet]{multinom}},
#'     \code{\link{detect_separation}},
#'     \code{\link[brglm2]{brmultinom}},
#'     \code{\link[brglm2]{bracl}}
#'
#' @references
#'
#' Lesaffre, E., & Albert, A. (1989). Partial Separation in Logistic
#' Discrimination. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **51**, 109-116
#'
#' Kosmidis I. and Firth D. (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71–82
#'
#' @examples
#'
#' # endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "detectseparation")
#' endometrial_ml <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                       family = binomial("probit"))
#' # clearly the maximum likelihood estimate for the coefficient of
#' # NV is infinite
#' (estimates <- check_infinite_estimates(endometrial_ml))
#' plot(estimates)
#'
#'
#' \donttest{
#' # Aligator data (Agresti, 2002, Table~7.1)
#' if (requireNamespace("brglm2", quietly = TRUE)) {
#'     data("alligators", package = "brglm2")
#'     all_ml <- brglm2::brmultinom(foodchoice ~ size + lake , weights = round(freq/3),
#'                          data = alligators, type = "ML", ref = 1)
#'     # Clearly some estimated standard errors diverge as the number of
#'     # Fisher scoring iterations increases
#'     plot(check_infinite_estimates(all_ml))
#'     # Bias reduction the brglm2 R packages can be used to get finite estimates
#'     all_br <- brglm2::brmultinom(foodchoice ~ size + lake , weights = round(freq/3),
#'                          data = alligators, ref = 1)
#'     plot(check_infinite_estimates(all_br))
#' }
#' }
#' @export
check_infinite_estimates.glm <- function(object, nsteps = 20, ...) {
    valid_classes <- c("glm", "brglmFit", "brmultinom")
    is_brmultinom <- inherits(object, "brmultinom")
    if (!inherits(object, valid_classes)) {
        warning("check_infinite_estimates has been designed for objects of class 'glm', 'brglmFit', 'brmultinom'")
    }
    if ((object$family$family != "binomial") & (!is_brmultinom)) {
        warning("check_infinite_estimates has been designed for binomial- or multinomial-response models")
    }
    if (is_brmultinom) {
        betas <- coef(object)
        dims <- dim(betas)
        betasNames <- paste(rep(colnames(betas), dims[1]), rownames(betas), sep = ":")
        betas <- c(betas)
        names(betas) <- betasNames
    }
    else {
        betas <- coef(object)
        betasNames <- names(betas)
    }
    eps <- .Machine$double.eps
    noNA <- !is.na(betas)
    stdErrors <- matrix(0, nsteps, length(betas))
    start <- NULL
    for (i in 1:nsteps) {
        if (is_brmultinom) {
            suppressWarnings(temp.object <- update(object, control = list(maxit = i, epsilon = eps, type = object$type), start = start))
            stdErrors[i, noNA] <- sqrt(diag(vcov(temp.object))[noNA])
        }
        else {
            suppressWarnings(temp.object <- update(object, control = list(maxit = i, epsilon = eps), start = start))
            stdErrors[i, noNA] <- summary(temp.object)$coef[betasNames[noNA], "Std. Error"]
        }
        start <- c(coef(temp.object))
    }
    res <- sweep(stdErrors, 2, stdErrors[1, ], "/")
    colnames(res) <- betasNames
    class(res) <- "inf_check"
    res
}

#' @export
plot.inf_check <- function(x, tol = 1e+2, ...) {
    # heuristic for determining ploting ranges
    sds <- apply(x, 2, sd)
    matplot(x, type = "l", lty = 1, ylim = range(x[, sds < tol]) * c(1, 1.5),
            ylab = "estimate", xlab = "number of iterations")
}
