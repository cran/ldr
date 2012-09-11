\name{ldr}
\alias{ldr}
\title{
Likelihood-based Dimension Reduction
}
\description{
Main function of the package. It creates objects of class \code{ldr} to estimate the central mean subspace. 
Implementation includes principal fitted components (PFC), likelihood acquired directions (LAD), 
and covariance reducing models (CORE). Adjunct functions can be used on the ldr object to estimate and 
test the dimension of the central mean subspace.
}
\usage{
ldr(X, y=NULL, fy=NULL, ycat=TRUE, numdir=2, 
model=c("core", "lad", "pfc"), verbose=FALSE, short=TRUE,...)
}
\arguments{
  \item{X}{Design matrix with \code{n} rows of observations and \code{p} columns of predictors. 
  The predictors are assumed to have a continuous distribution. 
}
  \item{y}{The response vector of length \code{n}. It must be numerical.}
  \item{fy}{Basis function to be obtained using \code{bf} or defined by the user. 
  It is a function of \code{y} alone and has indenpendent column vectors. 
  It is used exclusively with \code{pfc}. 
  See \code{\link{bf}} for detail.
}
  \item{numdir}{The number of directions to be used in estimating the reduction subspace. 
 When calling \code{pfc}, the dimension \code{numdir} must be less than or equal to the minimum of \code{p} and \code{r}, where 
 \code{r} is the number of columns of \code{fy}. When calling \code{lad} and \code{y} is continuous, \code{numdir} is the number of 
 slices to use.
}
  \item{ycat}{Boolean. Set to TRUE if \code{y} is categorical. It is used exclusively with \code{lad}.
  }
  \item{model}{One of the following: \code{pfc}, \code{lad}, \code{core}.
  }
    \item{verbose}{Boolean. If \code{TRUE}, the screen output provides some additional logging.
    }
  \item{short}{Boolean. If \code{TRUE}, the chosen model fits with the provided \code{numdir}. 
  If \code{FALSE}, the model is fit for all dimensions less or equal to \code{numdir}.
  }
  \item{\dots}{Additional arguments for specfic models and/or Grassmannoptim.
}
}
\details{
Likelihood-based methods to sufficient dimension reduction are model-based inverse regression approaches using the conditional distribution of 
the \eqn{p}-vector of predictors \eqn{X} given the response \eqn{Y=y}. Three methods are implemented in this package: covariance reduction (CORE),
principal fitted components (PFC), and likelihood acquired directions (LAD). All three assume that \eqn{X|(Y=y) \sim N(\mu_y, \Delta_y)}.

For CORE, given a set of \eqn{h} covariance matrices, the goal is to find a sufficient reduction that accounts for the heterogeneity 
among the population covariance matrices. See the documentation of \code{core} for details.

For PFC, \eqn{\mu_y=\mu + \Gamma \beta f_y}, with various structures of \eqn{\Delta}. The simplest is the isotropic 
(\code{iso}) with \eqn{\Delta=\delta^2 I_p}. 
The anisotropic (\code{aniso}) PFC model assumes that \eqn{\Delta=\mathrm{diag}(\delta_1^2, ..., \delta_p^2)}, 
where the conditional predictors are independent and on different measurement scales.
The unstructured (\code{unstr}) PFC model allows a general structure for \eqn{\Delta}. 
Extended structures are considered. See the help file of \code{pfc}
for more detail.

LAD assumes that the response \eqn{Y} is discrete. Continuous response are sliced into finite categories to meet this condition.  
It estimates the central subspace \eqn{\mathcal{S}_{Y|X}} by modeling both \eqn{\mu_y} and \eqn{\Delta_y}. See \code{lad} for more detail.
}
\value{ An object of class \code{ldr}. The output depends on the model used. See \code{pfc}, \code{lad}, and \code{core}.
}

\references{
Cook, RD (2007): Fisher Lecture - Dimension Reduction in Regression (with discussion). Statistical Science, 22, 1--26.

Cook, R. D. and Forzani, L. (2009).  Likelihood-based sufficient dimension reduction.  
Journal of the American Statistical Association, Vol. 104, 485, pp 197--208.

Cook, R. D. and Forzani, L. (2009).  Covariance reducing models: An alternative to spectral modelling of 
covariance matrices. Biometrika 95, 799-812.

Cook, R. D. and Forzani, L. (2009).  Principal fitted components for dimension reduction in regression. 
Statistical Science 23, 485--501.
}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\seealso{\code{\link{pfc}}, \code{\link{lad}}, \code{\link{core}}
}
\examples{
data(bigmac)
fit1 <- ldr(X=bigmac[,-1], fy=bf(y=bigmac[,1], case="poly", degree=3), 
numdir=3, structure="iso", model="pfc", short=FALSE)
summary(fit1, y=bigmac[,1])

fit2 <- ldr(X=bigmac[,-1], fy=bf(y=bigmac[,1], case="pdisc", degree=0, 
nslices=5), numdir=3, structure="unstr", model="pfc")
summary(fit2, y=bigmac[,1])

Gm <- ldr(X=bigmac[,-1], y=bigmac[,1], ycat=FALSE, model="lad", nslices=3)$Gammahat 
V<- t(Gm)\%*\%t(bigmac[,-1])
par(mfrow=c(2,1))
plot(bigmac[,1], V[1,], xlab="y", ylab="First Direction")
plot(bigmac[,1], V[2,], xlab="y", ylab="Second Direction") 
}

