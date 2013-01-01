\name{ldr}
\alias{ldr}
\title{
Likelihood-based Dimension Reduction
}
\description{
Main function of the package. It creates objects of one of classes \code{core}, \code{lad}, or \code{pfc} to estimate a sufficient dimension reduction subspace using covariance reducing models (CORE), likelihood acquired directions (LAD), or principal fitted components (PFC).
}

\usage{
ldr(X, y = NULL, fy = NULL, Sigmas = NULL, ns = NULL,
    numdir = NULL, nslices = NULL, model = c("core", "lad", "pfc"),
    numdir.test = FALSE, ...)
}

\arguments{
  \item{X}{Design matrix with \code{n} rows of observations and \code{p} columns of predictors. The predictors are assumed to have a continuous distribution.}
  \item{y}{The response vector of length \code{n}. It can be continuous or categorical.}
  \item{fy}{Basis function to be obtained using \code{bf} or defined by the user. It is a function of \code{y} alone and has independent column vectors. It is used exclusively with \code{pfc}.  See \code{\link{bf}} for detail.}
  \item{Sigmas}{A list object of sample covariance matrices corresponding to the different populations. It is used exclusively with \code{core}.}
  \item{ns}{A vector of number of observations of the samples corresponding to the different populations.}
  \item{numdir}{The number of directions to be used in estimating the reduction subspace. When calling \code{pfc}, the dimension \code{numdir} must be less than or equal to the minimum of \code{p} and \code{r}, where \code{r} is the number of columns of \code{fy}. When calling \code{lad} and \code{y} is continuous, \code{numdir} is the number of slices to use.}
  \item{nslices}{Number of slices for a continuous response. It is used exclusively with \code{lad}.}
  \item{model}{One of the following: \code{"pfc"}, \code{"lad"}, \code{"core"}.}
  \item{numdir.test}{Boolean. If \code{FALSE}, the chosen model fits with the provided \code{numdir}. If \code{TRUE}, the model is fit for all dimensions less or equal to \code{numdir}.}
  \item{\dots}{Additional arguments for specific models and/or Grassmannoptim.}
}

\details{
Likelihood-based methods to sufficient dimension reduction are model-based inverse regression approaches using the conditional distribution of the \eqn{p}-vector of predictors \eqn{X} given the response \eqn{Y=y}. Three methods are implemented in this package: covariance reduction (CORE), principal fitted components (PFC), and likelihood acquired directions (LAD). All three assume that \eqn{X|(Y=y) \sim N(\mu_y, \Delta_y)}.

For CORE, given a set of \eqn{h} covariance matrices, the goal is to find a sufficient reduction that accounts for the heterogeneity among the population covariance matrices. See the documentation of \code{"core"} for details.

For PFC, \eqn{\mu_y=\mu + \Gamma \beta f_y}, with various structures of \eqn{\Delta}. The simplest is the isotropic (\code{"iso"}) with \eqn{\Delta=\delta^2 I_p}. The anisotropic (\code{"aniso"}) PFC model assumes that \eqn{\Delta=\mathrm{diag}(\delta_1^2, ..., \delta_p^2)}, where the conditional predictors are independent and on different measurement scales. The unstructured (\code{"unstr"}) PFC model allows a general structure for \eqn{\Delta}. Extended structures are considered. See the help file of \code{pfc} for more detail.

LAD assumes that the response \eqn{Y} is discrete. A continuous response is sliced into finite categories to meet this condition. It estimates the central subspace \eqn{\mathcal{S}_{Y|X}} by modeling both \eqn{\mu_y} and \eqn{\Delta_y}. See \code{lad} for more detail.
}

\value{ An object of one of the classes \code{core}, \code{lad}, or \code{pfc} . The output depends on the model used. See \code{pfc}, \code{lad}, and \code{core} for further detail.}

\references{
Adragni, KP and Cook, RD (2009): Sufficient dimension reduction and prediction in regression. Phil. Trans. R. Soc. A 367, 4385-4405.

Cook, RD (2007): Fisher Lecture - Dimension Reduction in Regression (with discussion). Statistical Science, 22, 1--26.

Cook, R. D. and Forzani, L. (2008a).  Covariance reducing models: An alternative to spectral modelling of
covariance matrices. Biometrika 95, 799-812.

Cook, R. D. and Forzani, L. (2008b).  Principal fitted components for dimension reduction in regression.
Statistical Science 23, 485--501.

Cook, R. D. and Forzani, L. (2009).  Likelihood-based sufficient dimension reduction.
Journal of the American Statistical Association, Vol. 104, 485, pp 197--208.


}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\seealso{\code{\link{pfc}}, \code{\link{lad}}, \code{\link{core}}
}
\examples{
data(bigmac)
fit1 <- ldr(X=bigmac[,-1], y=bigmac[,1], fy=bf(y=bigmac[,1], case="pdisc",
        degree=0, nslices=5), numdir=3, structure="unstr", model="pfc")
summary(fit1)
plot(fit1)

fit2 <- ldr(X=bigmac[,-1], y=bigmac[,1], fy=bf(y=bigmac[,1], case="poly",
        degree=2), numdir=2, structure="aniso", model="pfc")
summary(fit2)
plot(fit2)

fit3 <- ldr(X=as.matrix(bigmac[,-1]), y=bigmac[,1], model="lad", nslices=5)
summary(fit3)
plot(fit3)
}

