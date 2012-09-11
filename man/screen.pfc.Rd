\name{screen.pfc}
\alias{screen.pfc}
\title{
Adaptive Screening of Predictors
}
\description{
Given a set of \eqn{p} predictors and a response, this function selects all predictors that are statistically related to the response at a specified significance level, using a flexible basis function.
}
\usage{
screen.pfc(X, fy)
}
\arguments{
  \item{X}{
Matrix or data frame with \code{n} rows of observations and \code{p} columns of predictors of continuous type. 
}
  \item{fy}{
Function of y. Basis function to be used to capture the dependency between individual predictors and the response. See \code{\link{bf}} for detail.
}

}
\details{
For each predictor \eqn{X_j}, write the equation
\deqn{X_j= \mu + \phi f_y + \epsilon}
where \eqn{f_y} is a flexible basis function provided by the user. The basis function is constructed using the function \code{bf}.
The screening procedure uses a test statistic on the null hypothesis \eqn{\phi=0} against the alternative \eqn{\phi \ne 0}.
}

\value{
Return a data frame object with \eqn{p} rows corresponding to the variables with the following columns
\item{F }{\code{F} statistic for testing the above hypotheses.}
\item{Pvalue }{The p-value of the test statistic. The \code{F} test has 1 and \code{n-2} degrees of freedom}
\item{Index }{Index of the variable, as its position \code{j}.}
}
\references{
Adragni, KP and Cook, RD (2008) Discussion on the Sure Independence Screening for Ultrahigh Dimensional Feature Space of Jianqing Fan and Jinchi Lv (2007) Journal of the Royal Statistical Society Series B, 70, Part5, pp1:35 
}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\examples{
data(OH)
X <- OH[, -c(1,295)]; y=OH[,295]

# Correlation screening
out <- screen.pfc(X, fy=bf(y, case="poly", degree=1))
head(out)

# Special basis function
out1 <- screen.pfc(X, fy=scale(cbind(y, sqrt(y)), center=TRUE, scale=FALSE))
head(out1)

# Piecewise constant basis with 10 slices
out2 <- screen.pfc(X, fy=bf(y, case="pdisc", degree=0, nslices=10))
head(out2)
}

