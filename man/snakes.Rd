\name{snakes}
\alias{snakes}
\docType{data}
\title{
Snakes data
}
\description{
Genetic covariance matrices for six genetic traits of two female garter snake populations, one from a coastal and the other from inland site in
northern California. The data set was initially studied by Phillips and Arnold (1999).
}
\usage{data(snakes)}
\format{
List format of 3 components.
\describe{
  \item{\code{snakes[[1]]}}{sample genetic covariance matrix for the inland population, obtained.}
  \item{\code{snakes[[2]]}}{sample genetic covariance matrix for the coastal population.}
  \item{\code{snakes[[3]]}}{vector of sample sizes, respectively for inland and coastal samples.}
}
}
\details{
Both genetic variance-covariances are obtained on six traits of the snakes.
}

\references{
Phillips P, Arnold S (1999). Hierarchical Comparison of Genetic variance-Covariance matrix
using the Flury Hierarchy." Evolution, 53, 1506--1515.
}
\examples{
data(snakes)
}
\keyword{datasets}
