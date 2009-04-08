\name{LBEMax}
\alias{LBEMax}
\alias{dlbe.max}
\alias{plbe.max}
\alias{qlbe.max}
\alias{rlbe.max}
\title{The Laguerre Beta Ensemble Maximum Eigenvalue Distributions}
\description{
  Density, distribution function, quantile function, and random
  generation for the maximum eigenvalue from the Laguerre Beta Ensemble
  (sample covariance matrix) with \code{nrow} rows, \code{ncol} columns,
  variance \code{var}, and order parameter \code{beta}.
}
\usage{
dlbe.max(x, nrow, ncol, var=1, beta=1, log = FALSE)
plbe.max(q, nrow, ncol, var=1, beta=1, lower.tail = TRUE, log.p = FALSE)
qlbe.max(p, nrow, ncol, var=1, beta=1, lower.tail = TRUE, log.p = FALSE)
rlbe.max(n, nrow, ncol, var=1, beta=1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{nrow, ncol}{the number of rows and columns in the data matrix}
  \item{var}{the variance of the elements in the data matrix}
  \item{beta}{the order parameter (1 or 2)}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dlbe.max} gives the density,
  \code{plbe.max} gives the distribution function,
  \code{qlbe.max} gives the quantile function, and
  \code{rlbe.max} generates random deviates.  
}
\details{
  If \code{beta} is not specified, it assumes the default value of \code{1}.
  Likewise, \code{var} assumes a default of \code{1}.
  
  The Laguerre Beta Ensemble (LBE) is a random sample covariance matrix from
  white noise with variance \code{var}.  The data matrix is assumed to
  have \code{nrow} observations in \code{ncol} dimensions.  Each element of
  the matrix is an independent normal random variable with variance \code{var}.
  These functions give the limiting distribution of the largest eigenvalue
  from such a matrix when \code{nrow} and \code{ncol} both tend to infinity.

  Supported values for \code{beta} are \code{1} for the Laguerre Orthogonal
  Ensemble and \code{2} for the Laguerre Unitary Ensemble.
  These values correspond to real and complex data, respectively.
}
\source{
  The functions are calculated by applying the appropriate centering and
  scaling (determined by \link{lbe.max}), and then calling the corresponding
  functions for the \link{TracyWidom} distribution.
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
\seealso{\link{lbe.max}, \link{TracyWidom}}
\keyword{distribution}
