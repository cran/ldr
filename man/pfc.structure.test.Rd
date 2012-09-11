\name{pfc.structure.test}
\alias{pfc.structure.test}
\title{
Likelihood ratio test for testing the structure of PFC models.
}
\description{
Given a PFC model, \code{pfc.structure.test} uses a likelihood ratio to test either isotropic against anisotropic or unstructured,  or anisotropic against unstructured.  
}
\usage{
pfc.structure.test(object1, object2, verbose = FALSE)
}
\arguments{
  \item{object1}{
A fitted object of \code{pfc} model.
}
  \item{object2}{
A fitted object of \code{pfc} model.
}
  \item{verbose}{
If TRUE, print steps on screen.
}
}
\value{
  \item{IC}{Akaike and Bayesian information criterion (AIC, and BIC) to select the appropriate covariance structure of the model.}
  \item{LRT}{Likelihood ratio test to test the covariance structure of the model.}
}
\references{
Adragni, K. P. and Cook, R. D. (2009). Sufficient dimension reduction and prediction in regression. Philosophical Transactions of the Royal Society, ser A., 367, 1906.
}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\examples{
data(bigmac)
fy=bf(y=bigmac[,1], case="poly", degree=3);
fit1 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="iso"); 
fit2 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="aniso"); 
fit3 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="unstr");
(pfc.structure.test(fit1, fit2))
(pfc.structure.test(fit2, fit3))

fit1 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="iso", short=FALSE); 
fit2 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="aniso", short=FALSE); 
fit3 <- pfc(X=bigmac[,-1], fy=fy,  numdir=2, structure="unstr", short=FALSE);
(pfc.structure.test(fit1, fit2))
(pfc.structure.test(fit2, fit3))
}
