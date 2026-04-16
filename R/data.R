#' Habitat preferences of lizards
#'
#' @format A data frame with 23 rows and 6 columns:
#'
#' * `grahami`. count of grahami lizards
#' * `opalinus`. count of opalinus lizards
#' * `height`. a factor with levels `<5ft`, `>=5ft`
#' * `diameter`. a factor with levels `<=2in`, `>2in`
#' * `light`. a factor with levels `sunny`, `shady`
#' * `time`. a factor with levels `early`, `midday`, `late`
#'
#' The variables `grahami` and `opalinus` are counts of two lizard
#' species at two different perch heights, two different perch
#' diameters, in sun and in shade, at three times of day.
#'
#' @seealso
#'
#' [detect_separation()]
#'
#' @source
#'
#' McCullagh P, Nelder J A (1989) _Generalized Linear
#' Models_ (2nd Edition).  London: Chapman and Hall.
#'
#' Originally from
#'
#' Schoener T W (1970) Nonsynchronous spatial overlap of lizards
#' in patchy habitats.  _Ecology_ *51*, 408-418.
#'
"lizards"

#' Histology grade and risk factors for 79 cases of endometrial cancer
#'
#' @format
#'
#' A data frame with 79 rows and 4 variables:
#'
#' * `NV`: neovasculization with coding 0 for absent and 1 for present
#' * `PI`: pulsality index of arteria uterina
#' * `EH`: endometrium height
#' * `HG`: histology grade with coding 0 for low grade and 1 for high grade
#'
#' @source
#'
#' The packaged data set was downloaded in `.dat` format from
#' <https://users.stat.ufl.edu/~aa/glm/data/>. The latter link
#' provides the data sets used in Agresti (2015).
#'
#' The endometrial data set was first analyzed in Heinze and
#' Schemper (2002), and was originally provided by Dr
#' E. Asseryanis from the Medical University of Vienna.
#'
#' @seealso
#'
#' [detect_separation()]
#'
#'
#' @references
#'
#' Agresti A (2015). *Foundations of Linear and Generalized Linear
#' Models*.  Wiley Series in Probability and Statistics. Wiley.
#'
#' Heinze G, Schemper M (2002). A Solution to the Problem of
#' Separation in Logistic Regression. *Statistics in Medicine*,
#' **21**, 2409–2419. \doi{10.1002/sim.1047}.
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82. \doi{10.1093/biomet/asaa052}.
#'
#'
"endometrial"

#' Separation Example Presented in Silvapulle (1981)
#'
#' @description
#' Separation example presented in Silvapulle (1981).
#'
#' @format
#'
#' A data frame with 35 rows and 2 variables:
#'
#' * `y`: a factor with the levels `case` and `none-case`, giving the
#'        outcome of a standardized psychiatric interview
#' * `ghqs`: an integer giving the general health questionnaire score.
#'
#' @seealso
#'
#' [detect_infinite_estimates()]
#'
#' @references
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' <https://www.jstor.org/stable/2984941>
#'
"silvapulle1981"
