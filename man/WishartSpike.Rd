\name{WishartSpike}
\alias{WishartSpike}
\alias{dwishart.spike}
\alias{pwishart.spike}
\alias{qwishart.spike}
\alias{rwishart.spike}
\title{The Spiked Wishart Maximum Eigenvalue Distributions}
\description{
  Density, distribution function, quantile function, and random
  generation for the maximum eigenvalue from a spiked Wishart matrix
  (sample covariance matrix) with \code{n.df} degrees of freedom,
  \code{p.dim} dimensions, and population covariance matrix 
  \code{diag(spike+var,var,var,...,var)}.
}
\usage{
dwishart.spike(x, spike, n.df=NA, p.dim=NA, var=1, log = FALSE)
pwishart.spike(q, spike, n.df=NA, p.dim=NA, var=1, lower.tail = TRUE, log.p = FALSE)
qwishart.spike(p, spike, n.df=NA, p.dim=NA, var=1, lower.tail = TRUE, log.p = FALSE)
rwishart.spike(n, spike, n.df=NA, p.dim=NA, var=1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{spike}{the value of the spike.}
  \item{n.df}{the number of degrees of freedom for the Wishart matrix.}
  \item{p.dim}{the number of dimensions (variables) for the Wishart matrix.}
  \item{var}{the population (noise) variance.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dwishart.spike} gives the density,
  \code{pwishart.spike} gives the distribution function,
  \code{qwishart.spike} gives the quantile function, and
  \code{rwishart.spike} generates random deviates.  
}
\details{
  The spiked Wishart is a random sample covariance matrix from
  multivariate normal data with \code{n.df} observations in \code{p.dim} 
  dimensions.  The spiked Wishart has one population covariance eigenvalue
  equal to \code{spike+var} and the rest equal to \code{var}.  These
  functions are related to the limiting distribution of the largest eigenvalue
  from such a matrix when \code{n.df} and \code{p.dim} both tending to
  infinity, with their ratio tending to a nonzero constant.

  For the spiked distribution to exist, \code{spike} must be greater than
  \code{sqrt(p.dim/n.df)*var}.
}
\references{
  Baik, J., Ben Arous, G., and \enc{Péché}{Peche}, S. (2005).
    Phase transition of the largest eigenvalue for non-null 
    complex sample covariance matrices.
    \emph{Annals of Probability} \bold{33}, 1643--1697.
    
  Baik, J. and Silverstein, J. W. (2006). 
    Eigenvalues of large sample covariance matrices of spiked 
    population models.
    \emph{Journal of Multivariate Analysis} \bold{97}, 1382-1408.
    
  Paul, D. (2007). Asymptotics of sample eigenstructure for a large 
    dimensional spiked covariance model.
    \emph{Statistica Sinica.} \bold{17},
    1617--1642.
}
\author{Patrick O. Perry}
\seealso{\link{wishart.spike.par}, \link{WishartMax}}
\keyword{distribution}
