\name{WishartMaxPar}
\alias{WishartMaxPar}
\title{White Wishart Maximum Eigenvalue Centering and Scaling}
\description{
  Centering and scaling for the maximum eigenvalue from a white Wishart
  matrix (sample covariance matrix) with with \code{ndf} degrees of freedom,
  \code{pdim} dimensions, population variance \code{var}, and order
  parameter \code{beta}.
}
\usage{
  WishartMaxPar(ndf, pdim, var=1, beta=1)
}
\arguments{
  \item{ndf}{the number of degrees of freedom for the Wishart matrix.}
  \item{pdim}{the number of dimensions (variables) for the Wishart matrix.}
  \item{var}{the population variance.}
  \item{beta}{the order parameter (1 or 2).}
}
\value{
  \item{centering}{gives the centering.}
  \item{scaling}{gives the scaling.}
}
\details{
  If \code{beta} is not specified, it assumes the default value of \code{1}.
  Likewise, \code{var} assumes a default of \code{1}.
  
  The returned values give appropriate centering and scaling for the largest
  eigenvalue from a white Wishart matrix so that the centered and scaled
  quantity converges in distribution to a Tracy-Widom random variable.  We
  use the second-order accurate versions of the centering and scaling given
  in the references below.
}
\references{
  El Karoui, N. (2006). A rate of convergence result for the largest 
    eigenvalue of complex white Wishart matrices. 
    \emph{Annals of Probability} \bold{34}, 2077--2117.

  Ma, Z. (2008). Accuracy of the Tracy-Widom limit for the largest eigenvalue 
    in white Wishart matrices.
    \emph{arXiv:0810.1329v1 [math.ST]}.
}
\author{Iain M. Johnstone, Zongming Ma, Patrick O. Perry and Morteza Shahram}
\seealso{\link{WishartMax}, \link{TracyWidom}}
\keyword{distribution}
