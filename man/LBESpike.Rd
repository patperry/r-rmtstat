\name{LBESpike}
\alias{LBESpike}
\alias{dlbe.spike.eig}
\alias{plbe.spike.eig}
\alias{qlbe.spike.eig}
\alias{rlbe.spike.eig}
\title{The Laguerre Beta Ensemble Maximum Eigenvalue Distributions}
\description{
  Density, distribution function, quantile function, and random
  generation for the maximum eigenvalue from the spiked Laguerre Beta Ensemble
  (sample covariance matrix) with \code{nrow} rows, \code{ncol} columns,
  with shape parameter \code{beta} and variance parameter \code{var}.
  The spike comes from the population covariance matrix being equal to
  \code{diag(spike+var,var,var,...,var)}.
}
\usage{
dlbe.spike.eig(x, spike, nrow=NA, ncol=NA, var=1, log = FALSE)
plbe.spike.eig(q, spike, nrow=NA, ncol=NA, var=1, lower.tail = TRUE, log.p = FALSE)
qlbe.spike.eig(p, spike, nrow=NA, ncol=NA, var=1, lower.tail = TRUE, log.p = FALSE)
rlbe.spike.eig(n, spike, nrow=NA, ncol=NA, var=1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{spike}{the population covariance eigenvalue corresponding to the spike.}
  \item{nrow, ncol}{the number of rows and columns in the data matrix}
  \item{var}{the variance parameter}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dlbe.spike.eig} gives the density,
  \code{plbe.spike.eig} gives the distribution function,
  \code{qlbe.spike.eig} gives the quantile function, and
  \code{rlbe.spike.eig} generates random deviates.  
}
\details{
  The Laguerre Beta Ensemble (LBE) is a random sample covariance matrix from
  white noise with \code{nrow} observations in \code{ncol} dimensions.  The
  spiked LBE has one population covariance eigenvalue equal to 
  \code{spike+var} and the rest equal to \code{var}.  These functions are
  related to the limiting distribution of the largest eigenvalue from such a 
  matrix when \code{nrow} and \code{ncol} both tending to infinity, with their
  ratio tending to a nonzero constant.

  For spiked distribution to exist, \code{spike} must be greater than
  \code{sqrt(ncol/nrow)*var}.
}
\references{
  Paul, D. (2007). Asymptotics of sample eigenstructure for a large 
    dimensional spiked covariance model.
    \emph{Statistica Sinica.} \bold{17}
    1617--1642.
}
\author{Patrick O. Perry}
\seealso{\link{lbe.spike.eig}, \link{WishartMax}}
\keyword{distribution}
