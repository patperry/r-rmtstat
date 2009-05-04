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
dmp( x, n.df=NA, p.dim=NA, var=1, concentration=n.df/p.dim, log = FALSE )
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{n.df}{the number of degrees of freedom for the Wishart matrix.}
  \item{p.dim}{the number of dimensions (variables) for the Wishart matrix.}
  \item{var}{the population variance.}
  \item{concentration}{the number of degrees of freedom per dimension.}   
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{dmp} gives the density.
}
\details{
  The \code{concentration} can either be given explicitly, or else computed
  from the given \code{n.df} and \code{p.dim}.  If \code{var} is not 
  specified, it assumes the default of \code{1}.
  
  The \enc{Marčenko}{Marcenko}-Pastur law is the limit of the random 
  probability measure which puts equal mass on all \code{p.dim} eigenvalues
  of a normalized \code{p.dim}-dimensional white Wishart matrix with 
  \code{n.df} degrees of freedom and scale parameter
  \code{diag(var, var, ..., var)}.  It is assumed that \code{n.df} goes to 
  infinity, and \code{n.df/p.dim} goes to nonzero constant called the 
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
