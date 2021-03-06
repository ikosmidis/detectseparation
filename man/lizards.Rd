% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lizards}
\alias{lizards}
\title{Habitat preferences of lizards}
\format{
An object of class \code{data.frame} with 23 rows and 6 columns.
}
\source{
McCullagh, P. and Nelder, J. A. (1989) _Generalized Linear
  Models_ (2nd Edition).  London: Chapman and Hall.

Originally from

    Schoener, T. W. (1970) Nonsynchronous spatial overlap of lizards
    in patchy habitats.  _Ecology_ *51*, 408-418.
}
\usage{
lizards
}
\description{
The lizards data frame has 23 rows and 6 columns. Variables
\code{grahami} and \code{opalinus} are counts of two lizard species
at two different perch heights, two different perch diameters, in
sun and in shade, at three times of day.
}
\details{
\itemize{
  \item grahami. count of grahami lizards
  \item opalinus. count of opalinus lizards
  \item height. a factor with levels \code{<5ft}, \code{>=5ft}
  \item diameter. a factor with levels \code{<=2in}, \code{>2in}
  \item light. a factor with levels \code{sunny}, \code{shady}
  \item time. a factor with levels \code{early}, \code{midday}, \code{late}
}
}
\seealso{
\code{\link[brglm2]{brglm_fit}}
}
\keyword{datasets}
