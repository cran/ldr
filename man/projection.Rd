\name{projection}
\alias{projection}
\title{
Projection operator
}
\description{
Compute a projection operator given the necessary arguments.
}
\usage{
projection(alpha, Sigma)
}
\arguments{
  \item{alpha}{
A semi-orthogonal matrix with \code{p} rows and \code{d} columns, \eqn{d \leq p}.
}
  \item{Sigma}{
A full rank \eqn{p \times p} square matrix.
}
}
\details{
The projection operator is obtained as \eqn{\alpha (\alpha^T \Sigma \alpha)^{-1}\alpha^T \Sigma}.
}
\author{
Kofi Placid Adragni
}
\keyword{internal}
