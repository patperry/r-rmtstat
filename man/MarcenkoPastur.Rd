\name{MarcenkoPastur}
\alias{MarcenkoPastur}
\alias{dmp}
\title{The Marcenko-Pastur Distribution}
\description{
  Density for the \enc{Marčenko}{Marcenko}-Pastur distribution, the
  limiting distribution of the empirical spectral measure for a large
  white Wishart matrix.
}
\usage{
dmp( x, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, log = FALSE )
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{ndf}{the number of degrees of freedom for the Wishart matrix.}
  \item{pdim}{the number of dimensions (variables) for the Wishart matrix.}
  \item{var}{the population variance.}
  \item{svr}{samples to variables ratio; the number of degrees of freedom per dimension.}   
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{dmp} gives the density.
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
  \code{concentration}.
}
\references{
  \enc{Marčenko}{Marcenko}, V.A. and Pastur, L.A. (1967).  Distribution of 
    eigenvalues for some sets of random matrices.
    \emph{Sbornik: Mathematics}
    \bold{1},
    457--483.
}
\author{Patrick O. Perry}
\keyword{distribution}
