\name{lad}
\alias{lad}
\title{
Likelihood Acquired Directions
}
\description{
Method to estimate the central subspace, using inverse conditional mean and conditional variance functions.
}
\usage{
lad(X, y, numdir = NULL, nslices = NULL, numdir.test = FALSE, ...)
}

\arguments{
  \item{X}{Data matrix with \code{n} rows of observations and \code{p} columns of predictors. The predictors are assumed to have a continuous distribution.}
  \item{y}{Response vector of \code{n} observations, possibly categorical or continuous. It is assumed categorical if \code{nslices=NULL}.}
  \item{numdir}{Integer between 1 and \eqn{p}. It is the number of directions of the reduction to estimate. If not provided then it will equal the number of distinct values of the categorical response.}
  \item{nslices}{Integer number of slices. It must be provided if \code{y} is continuous, and must be less than \eqn{n}. It is used to discretize the continuous response.}
  \item{numdir.test}{Boolean. If \code{FALSE}, \code{core} computes the reduction for the specific number of directions \code{numdir}. If \code{TRUE}, it does the computation of the reduction for the \code{numdir} directions, from 0 to \code{numdir}.}
  \item{\dots}{Other arguments to pass to \code{GrassmannOptim}.}
}

\details{
Consider a regression in which the response \eqn{Y} is discrete with support \eqn{S_Y=\{1,2,...,h\}}.
Following standard practice, continuous response can be sliced into finite categories to meet this condition.
Let \eqn{X_y \in R^p} denote a random vector of predictors distributed as \eqn{X|(Y=y)} and assume
that \eqn{X_y \sim N(\mu_y, \Delta_y)}, \eqn{y \in S_Y}. Let \eqn{\mu=E(X)} and \eqn{\Sigma=\mathrm{Var}(X)}
denote the marginal mean and variance of \eqn{X} and let \eqn{\Delta=E(\Delta_Y)} denote the average covariance matrix.
Given \eqn{n_y} independent observations of \eqn{X_y, y \in S_{Y}}, the goal is to obtain the maximum likelihood
estimate of the \eqn{d}-dimensional central subspace \eqn{\mathcal{S}_{Y|X}}, which is defined informally as the smallest
subspace such that \eqn{Y} is independent of \eqn{X} given its projection \eqn{P_{\mathcal{S}_{Y|X}}X}
onto \eqn{\mathcal{S}_{Y|X}}.

Let \eqn{\tilde{\Sigma}} denote the sample covariance matrix of \eqn{X}, let \eqn{\tilde{\Delta}_y} denote the sample
covariance matrix for the data with \eqn{Y=y}, and let \eqn{\tilde{\Delta}=\sum_{y=1}^{h} m_y \tilde{\Delta}_y} where \eqn{m_y}
is the fraction of cases observed with \eqn{Y=y}. The maximum likelihood estimator of \eqn{\mathcal{S}_{Y|X}} maximizes over
\eqn{\mathcal{S} \in \mathcal{G}_{(d,p)}} the log-likelihood function
%
\deqn{L(\mathcal{S})=\frac{n}{2}\log|P_{\mathcal{S}} \tilde{\Sigma} P_{\mathcal{S}}|_0 - \frac{n}{2}\log|\tilde{\Sigma}|
- \frac{1}{2}\sum_{y=1}^{h} n_y \log|P_{\mathcal{S}} \tilde{\Delta}_y P_{\mathcal{S}}|_0,}
%
where \eqn{|A|_0} indicates the product of the non-zero eigenvalues of a positive semi-definite symmetric
matrix \eqn{A}, \eqn{P_{\mathcal{S}}} indicates the projection onto the subspace \eqn{\mathcal{S}} in the usual inner product, and \eqn{\mathcal{G}_{(d,p)}} is the set of all \eqn{d}-dimensional subspaces in \eqn{R^p}, called Grassmann manifold. The desired reduction is then \eqn{\hat{\Gamma}^{T}X}. Once the dimension of the reduction subspace is estimated, the columns of \eqn{\hat{\Gamma}} are a basis for the maximum likelihood estimate of \eqn{\mathcal{S}_{Y|X}}.

The dimension \eqn{d} of the sufficient reduction is to be estimated. A sequential likelihood ratio test, and information criteria (AIC, BIC) are implemented, following Cook and Forzani (2009).
}
\value{
This command returns a list object of class \code{ldr}. The output depends on the argument \code{numdir.test}. If \code{numdir.test=TRUE}, a list of matrices is provided corresponding to the \code{numdir} values (1 through \code{numdir}) for each of the parameters \eqn{\Gamma}, \eqn{\Delta}, and \eqn{\Delta_y}; otherwise, a single list of matrices for a single value of \code{numdir}.
The output of \code{loglik}, \code{aic}, \code{bic}, \code{numpar} are vectors of \code{numdir} elements if \code{numdir.test=TRUE}, and scalars otherwise. Following are the components returned:
 \item{R}{The reduction data-matrix of \eqn{X} obtained using the centered data-matrix \eqn{X}. The centering of the data-matrix of \eqn{X} is such that each column vector is centered around its sample mean.}
 \item{Gammahat}{Estimate of \eqn{\Gamma}}
 \item{Deltahat}{Estimate of \eqn{\Delta}}
 \item{Deltahat_y}{Estimate of \eqn{\Delta_y}}
 \item{loglik}{Maximized value of the LAD log-likelihood.}
 \item{aic}{Akaike information criterion value.}
 \item{bic}{Bayesian information criterion value.}
 \item{numpar}{Number of parameters in the model.}
}
\references{
Cook RD, Forzani L (2009). Likelihood-based Sufficient Dimension Reduction, J. of the American Statistical Association, Vol. 104, No. 485, 197--208.}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}
\seealso{
\code{\link{core}, \link{pfc}}
}
\examples{
data(flea)
fit <- lad(X=flea[,-1], y=flea[,1], numdir=2, numdir.test=TRUE)
summary(fit)
plot(fit)
}

