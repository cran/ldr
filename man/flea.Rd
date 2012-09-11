\name{flea}
\alias{flea}
\docType{data}
\title{
Flea-beetles data
}
\description{
Six measurements on each of three species of flea-beetles: concinna, heptapotamica, and heikertingeri.
}
\usage{data(flea)}
\format{
  A data frame with 74 observations on the following 7 variables.
  \describe{
    \item{\code{species}}{a factor with levels \code{Concinna}, \code{Heikert.}, and \code{Heptapot.}}
    \item{\code{tars1}}{width of the first joint of the first tarsus in microns 
        (the sum of measurements for both tarsi).}
    \item{\code{tars2}}{the same for the second joint.}
    \item{\code{head}}{the maximal width of the head between the external edges of the
        eyes in 0.01 mm.}
    \item{\code{aede1}}{the maximal width of the aedeagus in the fore-part in microns.}
    \item{\code{aede2}}{the front angle of the aedeagus (1 unit = 7.5 degrees).}
    \item{\code{aede3}}{the aedeagus width from the side in microns.}
  }
}
\source{
Lubischew, AA "On the Use of Discriminant Functions in Taxonomy", Biometrics, Dec. 1962, pp. 455-477.
}
\references{Dianne Cook and Deborah F. Swayne, Interactive and Dynamic Graphics for Data Analysis: 
With Examples Using R and GGobi. URL: http://www.ggobi.org/book/data/flea.xml
}
\examples{
data(flea)
}
\keyword{datasets}
