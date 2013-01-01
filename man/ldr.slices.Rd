\name{ldr.slices}
\alias{ldr.slices}
\title{
Function to slice continuous response.
}
\description{
Divides a vector of length \code{n} into slices of approximately equal size. It is used to construct the piecewise bases, and
internally used in \code{lad} functions.
}
\usage{
ldr.slices(y, nslices = 3)
}
\arguments{
  \item{y}{
a vector of length \code{n}.
}
  \item{nslices}{
the number of slices, no larger than \code{n}.
}
}
\details{
The number of observations per slice \code{m} is computed as the largest integer less or equal to \code{n/nslices}.
The \code{n} observations of \code{y} are ordered in the increasing order.
The first set of first \code{m} observations is allocated to the first slice, the second set is allocated into the second slice, and so on.
}
\value{
Returns a named list with four elements as follows:
  \item{bins}{Slices with their observations}
  \item{nslices}{The actual number of slices produced.}
  \item{slice.size}{The number of observations in each slice.}
  \item{slice.indicator}{Vector of length \code{n} indicating the slice number of each observed response value.}
}
\references{
Cook, RD and Weisberg, S (1999), Applied Regression Including Computing and Graphics, New York: Wiley.
}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}
