\name{WishartMax}
\alias{WishartMax}
\alias{dwishart.max}
\alias{pwishart.max}
\alias{qwishart.max}
\alias{rwishart.max}
\title{The White Wishart Maximum Eigenvalue Distributions}
\description{
  Density, distribution function, quantile function, and random
  generation for the maximum eigenvalue from a white Wishart matrix
  (sample covariance matrix) with \code{n.df} degrees of freedom,
  \code{p.dim} dimensions, variance \code{var}, and order parameter
  \code{beta}.
}
\usage{
dwishart.max(x, n.df, p.dim, var=1, beta=1, log = FALSE)
pwishart.max(q, n.df, p.dim, var=1, beta=1, lower.tail = TRUE, log.p = FALSE)
qwishart.max(p, n.df, p.dim, var=1, beta=1, lower.tail = TRUE, log.p = FALSE)
rwishart.max(n, n.df, p.dim, var=1, beta=1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{n.df}{the number of degrees of freedom in the Wishart matrix}
  \item{p.dim}{the number of dimensions (variables) in the Wishart matrix}
  \item{var}{the variance of the elements in the data matrix}
  \item{beta}{the order parameter (1 or 2)}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dwishart.max} gives the density,
  \code{pwishart.max} gives the distribution function,
  \code{qwishart.max} gives the quantile function, and
  \code{rwishart.max} generates random deviates.  
}
\details{
  If \code{beta} is not specified, it assumes the default value of \code{1}.
  Likewise, \code{var} assumes a default of \code{1}.
  
  A white Wishart matrix is equal in distribution to \eqn{ (1/n) X' X },
  where \eqn{X} is an \eqn{n\times p} matrix with elements i.i.d. Normal
  with mean zero and variance \code{var}.  These functions give the limiting
  distribution of the largest eigenvalue from the such a matrix when 
  \code{n.df} and \code{p.dim} both tend to infinity.

  Supported values for \code{beta} are \code{1} for real data and
  and \code{2} for complex data.
}
\source{
  The functions are calculated by applying the appropriate centering and
  scaling (determined by \link{wishart.max.par}), and then calling the 
  corresponding functions for the \link{TracyWidom} distribution.
}
\references{
  Johansson, K. (2000). Shape fluctuations and random matrices.
    \emph{Communications in Mathematical Physics.} \bold{209}
    437--476.
    
  Johnstone, I.M. (2001). On the ditribution of the largest eigenvalue in
    principal component analysis. \emph{Annals of Statistics.} \bold{29}
    295--327.
}
\author{Patrick O. Perry}
\seealso{\link{wishart.max.par}, \link{WishartSpike}, \link{TracyWidom}}
\keyword{distribution}
