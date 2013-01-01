\name{structure.test}
\alias{structure.test}

\title{
Test of covariance structure for PFC models
}
\description{
Information criterion and likelihood ratio test for the structure of the covariance matrix of PFC models.
}
\usage{
structure.test(object1, object2)
}
\arguments{
  \item{object1}{
An object of class \code{pfc}
}
  \item{object2}{
A second object of class \code{pfc}, fitted exactly as for \code{object1} except for the covariance structure \eqn{\Delta}.
}
}
\details{
Consider two PFC models \eqn{M_1} and \eqn{M_2}, with the same parameters, except for the conditional covariance that is \eqn{\Delta_1} for \eqn{M_1} and \eqn{\Delta_2} for \eqn{M_2} such that model \eqn{M_1} is nested in model \eqn{M_2}. We implemented the likelihood ratio test for the hypotheses: \eqn{H_0: \Delta=\Delta_1} versus \eqn{H_a: \Delta=\Delta_2}. The test is implemented for the isotropic, anisotropic, and the unstructured PFC models. One may test isotropic against either anisotropic or unstructured, or test anisotropic against unstructured. The degrees of freedom are given by the difference in the number of parameters in the covariances. Information criterion AIC and BIC are also provided.

}

\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\examples{
data(bigmac)
fit1 <- pfc(X=bigmac[,-1], y=bigmac[,1], fy=bf(y=bigmac[,1], case="poly", 
        degree=3), numdir=3, structure="iso")
fit2 <- pfc(X=bigmac[,-1], y=bigmac[,1], fy=bf(y=bigmac[,1], case="poly", 
        degree=3), numdir=3, structure="aniso")
fit3 <- pfc(X=bigmac[,-1], y=bigmac[,1], fy=bf(y=bigmac[,1], case="poly", 
        degree=3), numdir=3, structure="unstr")
structure.test(fit1, fit3)
structure.test(fit2, fit3)
}
