\name{lad}
\alias{lad}
\title{
Likelihood Acquired Directions
}
\description{
Method to estimate the central subspace, using inverse conditional mean and conditional variance functions.
}
\usage{
lad(X, y, ycat=TRUE, numdir=2, nslices=NULL, verbose=FALSE, 
    short=TRUE, ...)
}
\arguments{
  \item{X}{Data matrix with \code{n} rows of observations and \code{p} columns of predictors. 
  The predictors are assumed to have a continuous distribution. 
  }
  \item{y}{Vector of numerical response observations.}
  \item{ycat}{Boolean. Must be \code{TRUE} if \eqn{y} is categorical.}
  \item{numdir}{Integer between 1 and \eqn{p}. It is the number of directions to use.}
  \item{nslices}{Integer. Number of slices, must be less than \eqn{n}}
  \item{verbose}{Boolean. If \code{TRUE}, the screen output provides some additional logging.}
  \item{short}{Boolean. If \code{TRUE}, \code{core} computes the reduction for the specific number of directions \code{numdir}. 
If \code{FALSE}, it does the computation of the reduction for the \code{numdir} directions, from 0 to \code{numdir}.}
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
where \eqn{|\mathbf{A}|_0} indicates the product of the non-zero eigenvalues of a positive semi-definite symmetric 
matrix \eqn{\mathbf{A}} and \eqn{P_{\mathcal{S}}} indicates the projection onto the subspace \eqn{\mathcal{S}} in the usual inner product.  
The desired reduction is then \eqn{\hat{\Gamma}^{T}X}, where the columns of \eqn{\hat{\Gamma}} are a basis for the maximum 
likelihood estimate of \eqn{\mathcal{S}_{Y|X}}.
}
\value{
This command returns a list object of class \code{ldr}. The output depends on the argument \code{short}. If \code{short=FALSE}, 
a list of matrices is provided corresponding to the \code{numdir} values (1 through \code{numdir}) for each of the parameters 
\eqn{\Gamma}, \eqn{\Delta}, and \eqn{\Delta_y}; otherwise, a single list of matrices for a single value of \code{numdir}. 
The output of \code{loglik}, \code{aic}, \code{bic}, \code{numpar} are vectors of \code{numdir} elements if \code{short=FALSE}, 
and scalars otherwise. Following are the components returned:
 \item{Gammahat}{Estimate of \eqn{\Gamma}}
 \item{Deltahat}{Estimate of \eqn{\Delta}}
 \item{Deltahat_y}{Estimate of \eqn{\Delta_y}}
 \item{loglik}{Maximized value of the CORE log-likelihood.}
 \item{aic}{Akaike information criterion value.}
 \item{bic}{Bayesian information criterion value.}
 \item{numpar}{Number of parameters in the model.}
}
\references{
Cook RD, Forzani L (2009). Likelihood-based Sufficient Dimension Reduction, J. of the American Statistical Association, Vol. 104, No. 485, 197--208.}
\author{
Kofi Placid Adragni <kofi@umbc.edu>
}

\examples{
\dontrun{
data(flea) 
X <- flea[,-1]
y<-as.integer(factor(flea[,1], levels=unique(flea[,1])))
fit <- lad(X=X, y=y, numdir=2, verbose=TRUE)
summary(fit)
V<- t(fit$Gammahat)\%*\%t(X)
mycol <- y
mycol[mycol==1]<-"black"
mycol[mycol==2]<-"blue"
mycol[mycol==3]<-"green"
plot(V[1,], V[2,], xlab="Dir1", ylab="Dir2", pch=as.numeric(y), 
cex=1.50, col=mycol)
}
}

