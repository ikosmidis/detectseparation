#' Habitat preferences of lizards
#'
#' The lizards data frame has 23 rows and 6 columns. Variables
#' \code{grahami} and \code{opalinus} are counts of two lizard species
#' at two different perch heights, two different perch diameters, in
#' sun and in shade, at three times of day.
#'
#' \itemize{
#'   \item grahami. count of grahami lizards
#'   \item opalinus. count of opalinus lizards
#'   \item height. a factor with levels \code{<5ft}, \code{>=5ft}
#'   \item diameter. a factor with levels \code{<=2in}, \code{>2in}
#'   \item light. a factor with levels \code{sunny}, \code{shady}
#'   \item time. a factor with levels \code{early}, \code{midday}, \code{late}
#' }
#'
#' @seealso
#'
#' \code{\link[brglm2]{brglm_fit}}

#'
#' @source
#'
#'   McCullagh, P. and Nelder, J. A. (1989) _Generalized Linear
#'   Models_ (2nd Edition).  London: Chapman and Hall.
#'
#' Originally from
#'
#'     Schoener, T. W. (1970) Nonsynchronous spatial overlap of lizards
#'     in patchy habitats.  _Ecology_ *51*, 408-418.
#'
"lizards"

#' Histology grade and risk factors for 79 cases of endometrial cancer
#'
#' @format A data frame with 79 rows and 4 variables:
#' \describe{
#'
#' \item{NV}{neovasculization with coding 0 for absent and 1 for present}
#'
#' \item{PI}{pulsality index of arteria uterina}
#'
#' \item{EH}{endometrium height}
#'
#' \item{HG}{histology grade with coding 0 for low grade and 1 for high grade}
#'
#' }
#'
#' @source The packaged data set was downloaded in \code{.dat} format
#'     from \url{https://users.stat.ufl.edu/~aa/glm/data/}. The latter
#'     link provides the data sets used in Agresti (2015).
#'
#'     The endometrial data set was first analyzed in Heinze and
#'     Schemper (2002), and was originally provided by Dr
#'     E. Asseryanis from the Medical University of Vienna.
#'
#' @seealso
#'
#' \code{\link[brglm2]{brglm_fit}}
#'
#'
#' @references
#'
#' Agresti, A. (2015). *Foundations of Linear and Generalized Linear
#' Models*.  Wiley Series in Probability and Statistics. Wiley
#'
#' Heinze, G., & Schemper, M. (2002). A Solution to the Problem of
#' Separation in Logistic Regression. *Statistics in Medicine*, **21**, 2409–2419
#'
"endometrial"


#' Separation Example Presented in Silvapulle (1981)
#'
#' @description
#' Separation example presented in Silvapulle (1981).
#'
#' @format A data frame with 35 rows and 2 variables:
#' \describe{
#'
#' \item{y}{a factor with the levels \code{case} and \code{none-case},
#'          giving the outcome of a standardized psychiatric interview}
#' \item{ghqs}{an integer giving the general health questionnaire score.}
#'
#' }
#'
#' @references
#'
#' Silvapulle, M. J. (1981).
#' On the Existence of Maximum Likelihood Estimators for the Binomial Response Models.
#' Journal of the Royal Statistical Society. Series B (Methodological), 43(3), 310–313.
#' \url{https://www.jstor.org/stable/2984941}
#'
"silvapulle1981"
