\name{lbe.spike}
\alias{lbe.spike.eig}
\alias{lbe.spike.vec}
\title{Spiked Laguerre Beta Ensemble Eigenvalue and Eigenvector Limits}
\description{
  Centering and scaling for the maximum eigenvalue from the spiked
  Laguerre Beta Ensample (sample covariance matrix) with \code{nrow} rows,
  \code{ncol} columns, with shape parameter \code{beta} and population
  covariance matrix \code{diag(spike+var,var,var,...,var)}.

  Almost-sure limit of the square of the dot product between the sample
  and the population eigenvector for the spike.
}
\usage{
  lbe.spike.eig( spike, nrow=NA, ncol=NA, var=1 )
  lbe.spike.vec( spike, nrow=NA, ncol=NA, var=1 )
}
\arguments{
  \item{spike}{the population covariance eigenvalue corresponding to the spike.}
  \item{nrow, ncol}{the number of rows and columns in the data matrix}
  \item{var}{the variance parameter}
}
\value{
  \item{center} gives the centering
  \item{scale}  gives the scale
}
\details{
  The returned values give appropriate centering and scaling for the largest
  eigenvalue from the spiked Laguerre Beta Ensemble so that the centered and 
  scaled quantity converges in distribution to a normal random variable
  with mean 0 and variance 1.
  
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
\seealso{\link{LBESpike}}
\keyword{distribution}
