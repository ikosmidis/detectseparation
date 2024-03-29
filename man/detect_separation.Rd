% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_separation.R, R/zzz_conventions.R
\name{detect_separation}
\alias{detect_separation}
\alias{detectSeparation}
\alias{print.detect_separation}
\title{Detect Separation}
\usage{
detect_separation(
  x,
  y,
  weights = NULL,
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = NULL,
  family = gaussian(),
  control = list(),
  intercept = TRUE,
  singular.ok = TRUE
)

detectSeparation(
  x,
  y,
  weights = NULL,
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = NULL,
  family = gaussian(),
  control = list(),
  intercept = TRUE,
  singular.ok = TRUE
)
}
\arguments{
\item{x}{\code{x} is a design matrix of dimension \code{n * p}.}

\item{y}{\code{y} is a vector of observations of length \code{n}.}

\item{weights}{an optional vector of \sQuote{prior weights} to be used
    in the fitting process.  Should be \code{NULL} or a numeric vector.}

\item{start}{currently not used.}

\item{etastart}{currently not used.}

\item{mustart}{currently not used.}

\item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link[stats]{offset}} terms can be
    included in the formula instead or as well, and if more than one is
    specified their sum is used.  See \code{\link[stats]{model.offset}}.}

\item{family}{a description of the error distribution and link
    function to be used in the model.  For \code{glm} this can be a
    character string naming a family function, a family function or the
    result of a call to a family function.  For \code{glm.fit} only the
    third option is supported.  (See \code{\link[stats]{family}} for details of
    family functions.)}

\item{control}{a list of parameters controlling separation
detection. See \code{\link{detect_separation_control}()} for
details.}

\item{intercept}{logical. Should an intercept be included in the
    \emph{null} model?}

\item{singular.ok}{logical. If \code{FALSE}, a singular model is an
error.}
}
\value{
A list that inherits from class \code{detect_separation},
\code{glm} and \code{lm}. A \code{print} method is provided for
\code{detect_separation} objects.
}
\description{
Method for \code{\link{glm}} that tests for data separation and
finds which parameters have infinite maximum likelihood estimates
in generalized linear models with binomial responses

\code{\link{detect_separation}()} is a method for \code{\link{glm}}
that tests for the occurrence of complete or quasi-complete
separation in datasets for binomial response generalized linear
models, and finds which of the parameters will have infinite
maximum likelihood estimates. \code{\link{detect_separation}()}
relies on the linear programming methods developed in Konis (2007).
}
\details{
Following the definitions in Albert and Anderson (1984), the data
for a binomial-response generalized linear model with logistic link
exhibit quasi-complete separation if there exists a non-zero
parameter vector \eqn{\beta} such that \eqn{X^0 \beta \le 0} and
\eqn{X^1 \beta \ge 0}, where \eqn{X^0} and \eqn{X^1} are the
matrices formed by the rows of the model matrix $X$ corresponding
to zero and non-zero responses, respectively. The data exhibits
complete separation if there exists a parameter vector \eqn{\beta} such
that the aforementioned conditions are satisfied with strict
inequalities. If there are no vectors \eqn{\beta} that can satisfy the
conditions, then the data points are said to overlap.

If the inverse link function \eqn{G(t)} of a generalized linear
model with binomial responses is such that \eqn{\log G(t)} and
\eqn{\log (1 - G(t))} are concave and the model has an intercept
parameter, then overlap is a necessary and sufficient condition for
the maximum likelihood estimates to be finite (see Silvapulle, 1981
for a proof). Such link functions are, for example, the logit,
probit and complementary log-log.

\code{\link{detect_separation}()} determines whether or not the
data exhibits (quasi-)complete separation. Then, if separation is
detected and the link function \eqn{G(t)} is such that \eqn{\log
G(t)} and \eqn{\log (1 - G(t))} are concave, the maximum likelihood
estimates has infinite components.

\code{\link{detect_separation}()} is a wrapper to the
\code{\link{detect_infinite_estimates}()} method. Separation
detection, as separation is defined above, takes place using the
linear programming methods in Konis (2007) regardless of the link
function. The output of those methods is also used to determine
which estimates are infinite, unless the link is "log". In the
latter case the linear programming methods in Schwendinger et
al. (2021) are called to establish if and which estimates are
infinite. If the link function is not one of `"logit"`, `"log"`,
`"probit"`, `"cauchit"`, `"cloglog"` then a warning is issued.

The \code{\link{coefficients}} method extracts a vector of values
for each of the model parameters under the following convention:
\code{0} if the maximum likelihood estimate of the parameter is
finite, and \code{Inf} or \code{-Inf} if the maximum likelihood
estimate of the parameter if plus or minus infinity. This
convention makes it easy to adjust the maximum likelihood estimates
to their actual values by element-wise addition.

\code{\link{detect_separation}()} can be passed directly as
a method to the \code{\link{glm}} function. See, examples.

\code{detectSeparation}() is an alias for \code{detect_separation}().
}
\note{
For the definition of complete and quasi-complete separation, see
Albert and Anderson (1984). Kosmidis and Firth (2021) prove that
the reduced-bias estimator that results by the penalization of the
logistic regression log-likelihood by Jeffreys prior takes always
finite values, even when some of the maximum likelihood estimates
are infinite. The reduced-bias estimates can be computed using the
\pkg{brglm2} R package.

\code{\link{detect_separation}} was designed in 2017 by Ioannis
Kosmidis for the **brglm2** R package, after correspondence with
Kjell Konis, and a port of the \code{separator} function had been
included in **brglm2** under the permission of Kjell Konis. In
2020, \code{\link{detect_separation}} and
\code{\link{check_infinite_estimates}} were moved outside
**brglm2** into the dedicated **detectseparation** package. Dirk
Schumacher authored the \code{separator_ROI} function, which
depends on the **ROI** R package and is now the default
implementation used for detecting separation. In 2022, Florian
Schwendinger authored the \code{dielb_ROI} function for detecting
infinite estimates in log-binomial regression, and, with Ioannis
Kosmidis, they refactored the codebase to properly accommodate for
the support of log-binomial regression.
}
\examples{

# endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "detectseparation")
endometrial_sep <- glm(HG ~ NV + PI + EH, data = endometrial,
                       family = binomial("logit"),
                       method = "detect_separation")
endometrial_sep
# The maximum likelihood estimate for NV is infinite
summary(update(endometrial_sep, method = "glm.fit"))

\donttest{
# Example inspired by unpublished microeconometrics lecture notes by
# Achim Zeileis https://eeecon.uibk.ac.at/~zeileis/
# The maximum likelihood estimate of sourhernyes is infinite
if (requireNamespace("AER", quietly = TRUE)) {
    data("MurderRates", package = "AER")
    murder_sep <- glm(I(executions > 0) ~ time + income +
                      noncauc + lfp + southern, data = MurderRates,
                      family = binomial(), method = "detect_separation")
    murder_sep
    # which is also evident by the large estimated standard error for NV
    murder_glm <- update(murder_sep, method = "glm.fit")
    summary(murder_glm)
    # and is also revealed by the divergence of the NV column of the
    # result from the more computationally intensive check
    plot(check_infinite_estimates(murder_glm))
    # Mean bias reduction via adjusted scores results in finite estimates
    if (requireNamespace("brglm2", quietly = TRUE))
        update(murder_glm, method = brglm2::brglm_fit)
}
}
}
\references{
Konis K. (2007). *Linear Programming Algorithms for Detecting
Separated Data in Binary Logistic Regression
Models*. DPhil. University of Oxford.
\url{https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}

Konis K. (2013). safeBinaryRegression: Safe Binary Regression. R
package version 0.1-3.
\url{https://CRAN.R-project.org/package=safeBinaryRegression}

Kosmidis I. and Firth D. (2021). Jeffreys-prior penalty, finiteness
and shrinkage in binomial-response generalized linear
models. *Biometrika*, **108**, 71–82. \doi{10.1093/biomet/asaa052}

Silvapulle, M. J. (1981).  On the Existence of Maximum Likelihood
Estimators for the Binomial Response Models.  *Journal of the Royal
Statistical Society. Series B (Methodological)*, **43**, 310–313.
\url{https://www.jstor.org/stable/2984941}

Schwendinger, F., Grün, B. & Hornik, K. (2021). A comparison of
optimization solvers for log binomial regression including conic
programming.  *Computational Statistics*, **36**,
1721–1754. \doi{10.1007/s00180-021-01084-5}
}
\seealso{
\code{\link{glm.fit}} and \code{\link{glm}}, \code{\link{detect_infinite_estimates}}, \code{\link{check_infinite_estimates}}, \code{\link[brglm2]{brglm_fit}}
}
\author{
Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}, Dirk Schumacher [aut] \email{mail@dirk-schumacher.net}, Florian Schwendinger [aut] \email{FlorianSchwendinger@gmx.at}, Kjell Konis [ctb] \email{kjell.konis@me.com}
}
