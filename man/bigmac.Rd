\name{bigmac}
\alias{bigmac}
\docType{data}
\title{
bigmac data
}
\description{The data give average values in 1991 on several economic indicators for 45 world cities.  All prices are in US dollars, using currency conversion at the time of publication.
}
\usage{data(bigmac)}
\format{
  A data frame with 45 observations on the following 10 variables.
  \describe{
    \item{\code{BigMac}}{Minimum labor to buy a BigMac and fries}
    \item{\code{Bread}}{Minimum labor to buy 1 kg bread}
    \item{\code{BusFare}}{Lowest cost of 10k public transit}
    \item{\code{EngSal}}{Electrical engineer annual salary, 1000s}
    \item{\code{EngTax}}{Tax rate paid by engineer}
    \item{\code{Service}}{Annual cost of 19 services}
    \item{\code{TeachSal}}{Primary teacher salary, 1000s}
    \item{\code{TeachTax}}{Tax rate paid by primary teacher}
    \item{\code{VacDays}}{Average days vacation per year}
    \item{\code{WorkHrs}}{Average hours worked per year}
  }
}
\source{
Rudolf Enz, "Prices and Earnings Around the Globe",
1991 edition, Published by the Union Bank of Switzerland.
}
\references{
Cook, RD and Weisberg, S (2004). Applied Regression Including Computing and Graphics, New York: Wiley, http://www.stat.umn.edu/arc. 
}
\examples{
data(bigmac)
pairs(bigmac)
}
\keyword{datasets}
