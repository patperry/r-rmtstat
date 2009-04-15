\name{wishart.spike.par}
\alias{wishart.spike.par}
\title{Spiked Wishart Eigenvalue Centering and Scaling}
\description{
  Centering and scaling for the a sample eigenvalue from a spiked
  Wishart matrix (sample covariance matrix) with \code{n.df} degrees of
  freedom, \code{p.dim} dimensions, and population covariance matrix 
  \code{diag(spike+var,var,var,...,var)}.
}
\usage{
  wishart.spike.par( spike, n.df=NA, p.dim=NA, var=1 )
}
\arguments{
  \item{spike}{the value of the spike.}
  \item{n.df}{the number of degrees of freedom in the Wishart matrix.}
  \item{p.dim}{the number of dimensions (variables) in the Wishart matrix.}
  \item{var}{the variance parameter.}
}
\value{
  \item{center} gives the centering
  \item{scale}  gives the scale
}
\details{
  The returned values give appropriate centering and scaling for the largest
  eigenvalue from the spiked Wishart matrix so that the centered and 
  scaled quantity converges in distribution to a normal random variable
  with mean 0 and variance 1.
  
  For spiked distribution to exist, \code{spike} must be greater than
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
\seealso{\link{WishartSpike}}
\keyword{distribution}
