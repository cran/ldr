\name{core}
\alias{core}
\title{Covariance Reduction}
\description{
Method to reduce sample covariance matrices to an informational core that is sufficient to characterize the variance 
heterogeneity among different populations. 
}
\usage{
core(X, y, Sigmas = NULL, ns = NULL, numdir = 2, verbose = FALSE, 
    short = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix with \code{n} rows of observations and \code{p} columns of predictors. The predictors are assumed to have a continuous distribution.
}
  \item{y}{
Vector of group labels. Observations with the same label are considered to be in the same group.
}
  \item{Sigmas}{
A list object of covariance matrices corresponding the different populations.
}
  \item{ns}{
A vector of number of observations of the samples corrsponding to the different populations.
}
  \item{numdir}{
Integer between 1 and \code{p}. It is the number of directions to estimate for the reduction.
}
  \item{verbose}{
Boolean. If \code{TRUE}, the screen output provides some additional logging.
}
  \item{short}{
Boolean. If \code{TRUE}, \code{core} computes the reduction for the specific number of directions \code{numdir}. 
If \code{FALSE}, it does the computation of the reduction for the \code{numdir} directions, from 0 to \code{numdir}.
}
  \item{\dots}{
Other arguments to pass to \code{GrassmannOptim}.
}
}

\details{
Consider the problem of characterizing the covariance matrices \eqn{\Sigma_y, y=1,...,h}, of a random vector \eqn{X} observed in each of \eqn{h}
normal populations. Let \eqn{S_y = (n_y-1)\tilde{\Sigma_y}} where \eqn{\tilde{\Sigma_y}} is the sample covariance matrix corresponding to 
\eqn{\Sigma_y}, and \eqn{n_y} is the number of observations corresponding to \eqn{y}. The goal is to find a semi-orthogonal matrix \eqn{\Gamma \in R^{p \times d}, d < p}, with the property that for any two populations \eqn{j} and 
\eqn{k}
\deqn{S_j|(\Gamma' S_j \Gamma=B, n_j=m) \sim S_k|(\Gamma' S_k \Gamma=B, n_k=m).}
That is, given \eqn{\Gamma' S_g \Gamma} and \eqn{n_g}, the conditional distribution of \eqn{S_g} must must depend on \eqn{g}. Thus 
\eqn{\Gamma' S_g \Gamma} is sufficient to account for the heterogeneity among the population covariance matrices.  
The central subspace \eqn{\mathcal{S}}, spanned by the columns of \eqn{\Gamma} is obtained by optimizing the following log-likelihood function
\deqn{L(\mathcal{S})= c-\frac{n}{2} \log|\tilde{\Sigma}| + \frac{n}{2} \log|P_{\mathcal{S}} \tilde{\Sigma} P_{\mathcal{S}}|-\sum_{y=1}^{h}\frac{n_y}{2} \log|P_{\mathcal{S}} \tilde{\Sigma}_y P_{\mathcal{S}}|,} 
where \eqn{c} is a constant depending only on \eqn{p} and \eqn{n_y}, \eqn{\tilde{\Sigma}_y, y=1,...,h,} denotes the sample covariance matrix
from population \eqn{y} computed with divisor \eqn{n_y}, and \eqn{\tilde{\Sigma}=\sum_{y=1}^{h} (n_y/n)\tilde{\Sigma}}.
}
\value{
This command returns a list object of class \code{ldr}. The output depends on the argument \code{short}. If \code{short=FALSE}, 
a list of matrices is provided corresponding to the \code{numdir} values (1 through \code{numdir}) for each of the parameters 
\eqn{\Gamma}, \eqn{\Sigma}, and \eqn{\Sigma_g}. Otherwise, a single list of matrices for a single value of \code{numdir}. 
The output of \code{loglik}, \code{aic}, \code{bic}, \code{numpar} are vectors with \code{numdir} elements if \code{short=FALSE}, 
and scalars otherwise. Following are the components returned:
 \item{Gammahat}{Estimate of \eqn{\Gamma}.}
 \item{Sigmahat}{Estimate of overall \eqn{\Sigma}.}
 \item{Sigmashat}{Estimate of group-specific \eqn{\Sigma_g}'s.}
 \item{loglik}{Maximized value of the CORE log-likelihood.}
 \item{aic}{Akaike information criterion value.}
 \item{bic}{Bayesian information criterion value.}
 \item{numpar}{Number of parameters in the model.}
}

\references{
Cook RD and Forzani L (2008). Covariance reducing models: An alternative to spectral modelling of covariance matrices. 
Biometrika, Vol. 95, No. 4, 799--812.
}
\author{
Andrew Raim, University of Maryland, Baltimore County
}
\note{
Currently \code{loglik}, \code{AIC}, and \code{BIC} are computed up to a constant. 
Therefore, these can be compared relatively (e.g. two \code{loglik}'s can be subtracted to compute a likelihood ratio test), 
but they should not be treated as absolute quantities. 
}

\seealso{
\code{\link{lad}, \link{pfc}}
}
\examples{
\dontrun{
data(flea)
X <- flea[,-1]
y<-as.integer(factor(flea[,1], levels=unique(flea[,1])))
fit <- core(X=as.matrix(X), y=y, numdir=1, verbose=TRUE, short=FALSE)
summary(fit)
}
}

