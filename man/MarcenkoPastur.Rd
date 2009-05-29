\name{MarcenkoPastur}
\alias{MarcenkoPastur}
\alias{dmp}
\alias{pmp}
\alias{qmp}
\alias{rmp}
\title{The Marcenko-Pastur Distribution}
\description{
  Density, distribution function, quantile function and random generation for 
  the \enc{Marčenko}{Marcenko}-Pastur distribution, the limiting distribution 
  of the empirical spectral measure for a large white Wishart matrix.
}
\usage{
dmp( x, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, log = FALSE )
pmp( q, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, lower.tail = TRUE, log.p = FALSE )
qmp( p, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, lower.tail = TRUE, log.p = FALSE )
rmp( n, ndf=NA, pdim=NA, var=1, svr=ndf/pdim )
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{ndf}{the number of degrees of freedom for the Wishart matrix.}
  \item{pdim}{the number of dimensions (variables) for the Wishart matrix.}
  \item{var}{the population variance.}
  \item{svr}{samples to variables ratio; the number of degrees of freedom per dimension.}   
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dmp} gives the density,
  \code{pmp} gives the distribution function,
  \code{qmp} gives the quantile function, and
  \code{rmp} generates random deviates.  
}
\details{
  The \code{concentration} can either be given explicitly, or else computed
  from the given \code{ndf} and \code{pdim}.  If \code{var} is not 
  specified, it assumes the default of \code{1}.
  
  The \enc{Marčenko}{Marcenko}-Pastur law is the limit of the random 
  probability measure which puts equal mass on all \code{pdim} eigenvalues
  of a normalized \code{pdim}-dimensional white Wishart matrix with 
  \code{ndf} degrees of freedom and scale parameter
  \code{diag(var, var, ..., var)}.  It is assumed that \code{ndf} goes to 
  infinity, and \code{ndf/pdim} goes to nonzero constant called the 
  "samples-to-variables ratio" (\code{svr}).
}
\source{
  Other than the density, these functions are relatively slow and imprecise.
  
  The distribution function is computed with \link{integrate}.
  The quantiles are computed via bisection using \link{uniroot}.
  Random variates are generated using the inverse CDF.
}

\references{
  \enc{Marčenko}{Marcenko}, V.A. and Pastur, L.A. (1967).  Distribution of 
    eigenvalues for some sets of random matrices.
    \emph{Sbornik: Mathematics}
    \bold{1},
    457--483.
}
\author{Iain M. Johnstone, Zongming Ma, Patrick O. Perry and Morteza Shahram}
\keyword{distribution}
