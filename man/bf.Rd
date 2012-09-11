\name{bf}
\alias{bf}
\title{Function to generate a basis function.  
}
\description{
This function is to construct a data matrix of basis function using the \code{n} response observations. 
The response can be continuous or categorical. The function returns a matrix of \code{n} rows and \code{r} columns. 
The number of columns \code{r} depends on the choice of basis function. 
Polynomial, piecewise polynomial continuous and discontinuous, and Fourier bases are implemented. 
For a polynomial basis, \code{r} is the degree of the polynomial. 
}
\usage{
bf(y, case = c("poly", "categ", "fourier", "pcont", "pdisc"),
degree = 1, nslices = 5, scale=FALSE)
}

\arguments{
  \item{y}{
A vector of \code{n} observations of the response. 
}
  \item{case}{
Take values \code{poly} for polynomial, \code{categ} for categorical, \code{fourier}, 
\code{pcont} for piecewise continuous, and \code{pdisc} for piecewise discountinuous bases. 
}
  \item{degree}{
For polynomial and piecewise polynomial bases, \code{degree} is the degree of the polynomial. 
With \code{pdisc}, \code{degree=0} corresponds to piecewise constant.  
}
  \item{nslices}{
The number of slices for piecewise bases.
}
  \item{scale}{
If TRUE, the columns of the basis function are scaled to have unit variance.
}
}
\details{
The basis function \eqn{f_y} is a vector-valued function of the response \eqn{y \in R}. 
There is an infinite number of basis functions, including the polynomial, piecewise polynomial, and Fourier. We implemented the following:

1. Polynomial Basis: \eqn{f_y=(y, y^2, ..., y^r)^T}. It corresponds to the \code{poly} argument of \code{bf}. 
The argument \code{degree} is \eqn{r} of the polynomial is provided by the user. 
The subsequent \eqn{n \times r} data-matrix is column-wise centered.

2. Piecewise Constant Basis: It corresponds to \code{pdisc} with \code{degree=0}. It is obtained by first slicing the range of \eqn{y} 
into \eqn{h} slices \eqn{H_1,...,H_k}. The \eqn{k^{th}} component \eqn{f_y \in \mathrm{R}^{h-1}} is 
\eqn{f_{y_k}=J(y \in H_k)-n_k/n, k=1, ..., h-1}, where \eqn{n_y} is the number of observations in \eqn{H_k}, and 
\eqn{J} is the indicator function. 

3. Piecewise Discontinuous Linear Basis: It corresponds to \code{pdisc} with \code{degree=1}. 
It is more elaborate than the piecewise constant basis. A linear function of \eqn{y} is fit within each slice. 
Let \eqn{\tau_i} be the knots, or endpoints of the slices. The components of \eqn{f_y \in \mathrm{R}^{2h-1}} are obtained with
    \eqn{f_{y_{(2i-1)}} = J(y \in H_i)}; \eqn{f_{y_{2i}} =  J(y \in H_i)(y-\tau_{i-1})} for \eqn{i=1,2,...,h-1} and     
    \eqn{f_{y_{(2h-1)}}  =  J(y \in H_{h})(y-\tau_{h-1})}. The subsequent \eqn{n \times (2h-1)} data-matrix is column-wise centered.

4. Piecewise Continuous Linear Basis: The general form of the components \eqn{f_{y_i}} of 
\eqn{f_y \in \mathrm{R}^{h+1}} is given by \eqn{f_{y_1}  =  J(y \in H_1)} and 
    \eqn{f_{y_{i+1}} =  J(y \in H_{i})(y-\tau_{i-1})} for \eqn{i=1,...,h.}. The subsequent \eqn{n \times (h-1)} data-matrix is column-wise centered. 
    This case corresponds to \code{pcont} with {degree=1}.

5. Fourier Basis: They consist of a series of pairs of sines and cosines of increasing frequency. 
A Fourier basis is given by \eqn{f_y=(\cos(2\pi y), \sin(2\pi y),..., \cos(2\pi ky), \sin(2\pi ky))^T.}
The subsequent \eqn{n \times 2k} data-matrix is column-wise centered.

6. Categorical Basis: It is obtained using \code{categ} option when \eqn{y} takes \eqn{h} distinct values \eqn{1, 2,..., h}, corresponding to 
the number of sub-populations or sub-groups. The number of slices is naturally \eqn{h}. The expression for the basis is identical to piecewise constant basis.
}

\value{
\item{fy}{A matrix with \code{n} rows and \code{r} columns.
}
\item{scale}{Boolean. If TRUE, the columns of the output are standardized to have unit variance.
}
}
\references{
Adragni, KP (2009) PhD Dissertation, University of Minnesota.

Adragni, KP and Cook, RD (2009): Sufficient dimension reduction and prediction in regression. Phil. Trans. R. Soc. A 367, 4385-4405.

Cook, RD (2007): Fisher Lecture - Dimension Reduction in Regression (with discussion).  Statistical Science, Vol. 22, 1--26.
}

\author{
Kofi Placid Adragni <kofi@umbc.edu>
}
\examples{
data(bigmac)

# Piecewise constant basis with 5 slices
fy=bf(y=bigmac[,1], case="pdisc", degree=0, nslices=5)
fit1 <- pfc(X=bigmac[,-1], fy=fy, numdir=3, structure="aniso", short=FALSE)
summary(fit1, y=bigmac[,1])

# Cubic polynomial basis
fy=bf(y=bigmac[,1], case="poly", degree=3)
fit2 <- pfc(X=bigmac[,-1], fy=fy, numdir=3, structure="aniso", short=FALSE)
summary(fit2, y=bigmac[,1])

# Piecewise linear continuous with 3 slices
fy=bf(y=bigmac[,1], case="pcont", degree=1, nslices=3)
fit3 <- pfc(X=bigmac[,-1], fy=fy, numdir=3, structure="unstr", short=FALSE)
summary(fit3, y=bigmac[,1])
}

