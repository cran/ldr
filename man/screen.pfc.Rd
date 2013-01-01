\name{screen.pfc}
\alias{screen.pfc}
\title{
Adaptive Screening of Predictors
}
\description{
Given a set of \eqn{p} predictors and a response, this function selects all predictors that are statistically related to the response at a specified significance level, using a flexible basis function.
}
\usage{
screen.pfc(X, fy, cutoff=0.1)
}
\arguments{
  \item{X}{
Matrix or data frame with \code{n} rows of observations and \code{p} columns of predictors of continuous type.
}
  \item{fy}{
Function of \code{y}. Basis function to be used to capture the dependency between individual predictors and the response. See \code{\link{bf}} for detail.
}
  \item{cutoff}{
The level of significance to be used for the cutoff, by default 0.1.
}
}
\details{
For each predictor \eqn{X_j}, write the equation
\deqn{X_j= \mu + \phi f_y + \epsilon}
where \eqn{f_y} is a flexible basis function provided by the user. The basis function is constructed using the function \code{bf}.
The screening procedure uses a test statistic on the null hypothesis \eqn{\phi=0} against the alternative \eqn{\phi \ne 0}. Given the \eqn{r} components of the basis function \eqn{f_y}, the above model is a linear model where \eqn{X_j} is the response and \eqn{f_y} constitutes the predictors. The hypothesis test on \eqn{\phi} is essentially an F-test. Specifically, given the data, let \eqn{\hat{\phi}} be the ordinary least squares estimator of \eqn{\phi}. We consider the usual test statistic
%
\deqn{F_j=\frac{n-r-1}{r}.\frac{\sum_{i=1}^n [(X_{ji}-\bar{X}_{j.})^2 -  (X_{ji}-\bar{X}_{j.} - \hat{\phi}_j \mathbf{f}_{y_i})^2]}{\sum_{i=1}^n (X_{ji}-\bar{X}_{j.} - \hat{\phi}_j \mathbf{f}_{y_i})^2}}
%
where \eqn{\bar{X}_{j.}=\sum_{i=1}^n X_{ji}/n}. The statistic \eqn{F_j} follows an \eqn{F} distribution with \eqn{(r, n-r-1)} degrees of freedom. The sample size \eqn{n} is expected to be larger than \eqn{r}.
}

\value{
Return a data frame object with \eqn{p} rows corresponding to the variables with the following columns
\item{F }{\code{F} statistic for testing the above hypotheses.}
\item{P-value }{The p-value of the test statistic. The \code{F} test has 1 and \code{n-2} degrees of freedom}
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

