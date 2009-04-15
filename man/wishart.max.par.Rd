\name{wishart.max.par}
\alias{wishart.max.par}
\title{The White Wishart Maximum Eigenvalue Centering and Scaling}
\description{
  Centering and scaling for the maximum eigenvalue from a white Wishart
  matrix (sample covariance matrix) with with \code{n.df} degrees of freedom,
  \code{p.dim} dimensions, variance \code{var}, and order parameter
  \code{beta}.
}
\usage{
  wishart.max.par(n.df, p.dim, var=1, beta=1)
}
\arguments{
  \item{n.df}{the number of degrees of freedom in the Wishart matrix}
  \item{p.dim}{the number of dimensions (variables) in the Wishart matrix}
  \item{var}{the variance of the elements in the data matrix}
  \item{beta}{the order parameter (1 or 2)}
}
\value{
  \item{center}{gives the centering}
  \item{scale}{gives the scaling}
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
    \emph{Annals of Probability}. 
    \bold{34} 2077--2117.

  Ma, Z. (2008). Accuracy of the Tracy-Widom limit for the largest eigenvalue 
    in white Wishart matrices.
    \emph{arXiv:0810.1329v1 [math.ST]}.
}
\author{Patrick O. Perry}
\seealso{\link{WishartMax}, \link{TracyWidom}}
\keyword{distribution}
