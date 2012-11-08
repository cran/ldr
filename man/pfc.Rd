\name{pfc}
\alias{pfc}
\title{Principal fitted components}

\description{
Principal fitted components model for sufficient dimension reduction. This function estimates all parameters in the model.
}
\usage{
pfc(X, fy=fy, numdir=2, structure = c("iso", "aniso", "unstr", 
"unstr2"), eps_aniso=0.01, verbose=FALSE, short=TRUE, ...)
}
\arguments{
  \item{X}{
Design matrix with \code{n} rows of observations and \code{p} columns of predictors. 
The predictors are assumed to have a continuous distribution. 
}
  \item{fy}{
Basis function to be obtained using \code{bf} or defined by the user. 
It is a function of \code{y} alone and has \code{r} independent column vectors. 
See \code{\link{bf}}, for detail.
}
  \item{numdir}{
The number of directions to be used in estimating the reduction subspace. 
The dimension must be less than or equal to the minimum of \code{r} and \code{p}.
}
  \item{structure}{
Structure of \code{var(X|Y)}. The following options are available: \code{"iso"} for isotropic
 (predictors, conditionally on the response, are independent and on the same measurement scale); 
 \code{"aniso"} for anisotropic (predictors, conditionally on the response, are independent and on different measurement scales); 
 \code{"unstr"} for unstructured variance. The fourth structure \code{"unstr2"} refers to an extended PFC model
 with an heterogenous error structure.
}
  \item{eps_aniso}{Precision term used in estimating \code{var(X|Y)} for the anisotropic structure.
}
  \item{verbose}{Boolean. If \code{TRUE}, the screen output provides some additional logging.}
  \item{short}{Boolean. If \code{TRUE}, \code{pfc} fits with the \code{numdir} provided only. 
  If \code{FALSE}, PFC models are fit for all dimensions less than or equal to \code{numdir}.
}
  \item{\dots}{
Additional arguments to \code{Grassmannoptim}.
}
}
\details{
Let \eqn{X} be a column vector of \eqn{p} predictors, and \eqn{Y} be a univariate response variable. 
Principal fitted components model is an inverse regression model for sufficient dimension reduction. It is an inverse regression model
given by \eqn{X|(Y=y) \sim N(\mu + \Gamma \beta f_y, \Delta)}. The term \eqn{\Delta} is assumed independent of \eqn{y}. 
Its simplest structure is the isotropic (\code{iso}) with \eqn{\Delta=\delta^2 I_p}, 
where, conditionally on the response, the predictors are independent and are on the same measurement scale. 
The sufficient reduction is \eqn{\Gamma^TX}. The anisotropic (\code{aniso}) PFC model assumes that 
\eqn{\Delta=}diag\eqn{(\delta_1^2, ..., \delta_p^2)}, where the conditional predictors are independent and on different measurement scales.
The unstructured (\code{unstr}) PFC model allows a general structure for \eqn{\Delta}. With the anisotropic and unstructured \eqn{\Delta}, the 
sufficient reduction is \eqn{\Gamma^T \Delta^{-1}X}. it should be noted that \eqn{X \in R^{p}} while the data-matrix to use is in \eqn{R^{n \times p}}.] 

The error structure of the extended structure has the following form 
\deqn{\Delta=\Gamma \Omega \Gamma^T + \Gamma_0 \Omega_0 \Gamma_0^T,}
%
where \eqn{\Gamma_0} is the orthogonal completion of \eqn{\Gamma} such that \eqn{(\Gamma, \Gamma_0)} is a 
\eqn{p \times p} orthogonal matrix. The matrices \eqn{\Omega \in R^{d \times d}} and \eqn{\Omega_0 \in 
R^{(p-d) \times (p-d)}} are assumed to be symmetric and full-rank. The sufficient reduction is \eqn{\Gamma^{T}X}. 
Let \eqn{\mathcal{S}_{\Gamma}} be the subspace spanned by the columns of \eqn{\Gamma}. The parameter space of \eqn{\mathcal{S}_{\Gamma}} 
is the set of all \eqn{d} dimensional subspaces in \eqn{R^p}, called Grassmann manifold
 and denoted by \eqn{\mathcal{G}_{(d,p)}}. 
Let \eqn{\hat{\Sigma}}, \eqn{\hat{\Sigma}_{\mathrm{fit}}} be the sample variance of \eqn{X} and 
the fitted covariance matrix, and let \eqn{\hat{\Sigma}_{\mathrm{res}}=\hat{\Sigma} - \hat{\Sigma}_{\mathrm{fit}}}. The 
 MLE of \eqn{\mathcal{S}_{\Gamma}} under \code{unstr2} setup is obtained by maximizing the log-likelihood 
\deqn{L(\mathcal{S}_U) =  - \log|U^T \hat{\Sigma}_{\mathrm{res}} U| - \log|V^T \hat{\Sigma}V|}
over \eqn{\mathcal{G}_{(d,p)}}, where \eqn{V} is an orthogonal completion of \eqn{U}.
}
\value{
This command returns a list object of class \code{ldr}. The output depends on the argument \code{short}. 
If \code{short=FALSE}, a list of matrices is provided corresponding to the \code{numdir} values (1 through \code{numdir})
 for each of the parameters \eqn{\mu}, \eqn{\beta}, \eqn{\Gamma}, \eqn{\Gamma_0}, \eqn{\Omega}, and \eqn{\Omega_0}. 
 Otherwise, a single list of matrices for a single value of \code{numdir}. The outputs of \code{loglik}, \code{aic}, \code{bic}, \code{numpar} are 
 vectors of \code{numdir} elements if \code{short=FALSE}, and scalars otherwise. Following are the components returned:
  \item{Muhat}{Estimate of \eqn{\mu}.}
  \item{Betahat}{Estimate of \eqn{\beta}.}
  \item{Deltahat}{The estimate of the covariance \eqn{\Delta}.}
  \item{Gammahat}{An estimated orthogonal basis representative of \eqn{\hat{\mathcal{S}}_{\Gamma}}, the subspace spanned by \eqn{\Gamma}.}
  \item{Gammahat0}{An estimated orthogonal basis representative of \eqn{\hat{\mathcal{S}}_{\Gamma_0}}, the subspace spanned by \eqn{\Gamma_0}.}
  \item{Omegahat}{The estimate of the covariance \eqn{\Omega} if an extended model is used.}
  \item{Omegahat0}{The estimate of the covariance \eqn{\Omega_0} if an extended model is used.}
  \item{loglik }{The value of the log-likelihood for the model.}
  \item{aic}{Akaike information criterion value.}
  \item{bic}{Bayesian information criterion value.}
  \item{numdir}{The number of directions to estimate.}
  \item{numpar}{The number of parameters in the model.}
  \item{evalues}{The first \code{numdir} largest eigenvalues of \eqn{\hat{\Sigma}_{\mathrm{fit}}}.}
}
\references{
Cook, RD (2007): Fisher Lecture - Dimension Reduction in Regression (with discussion).  Statistical Science, 22, 1--26.

Cook, RD and Forzani, L (2009): Principal fitted components for dimension reduction in regression. 
Statistical Science 23, 485--501.
}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\examples{
data(bigmac)
fit1 <- pfc(X=bigmac[,-1], fy=bf(y=bigmac[,1], case="poly", degree=3), 
numdir=3, structure="iso")
summary(fit1, y=bigmac[,1])

fit2 <- pfc(X=bigmac[,-1], fy=bf(y=bigmac[,1], case="poly", degree=3), 
numdir=3, structure="aniso", short=FALSE)
summary(fit2, y=bigmac[,1])

fit3 <- pfc(X=bigmac[,-1], fy=bf(y=bigmac[,1], case="poly", degree=3), 
numdir=3, structure="unstr", short=FALSE)
summary(fit3)
summary(fit3, y=bigmac[,1])	
}

