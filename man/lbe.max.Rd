\name{lbe.max}
\alias{lbe.max}
\title{The Laguerre Beta Ensemble Maximum Eigenvalue Centering and Scaling}
\description{
  Centering and scaling for the maximum eigenvalue from the Laguerre
  Beta Ensample (sample covariance matrix) with \code{nrow} rows,
  \code{ncol} columns, variance \code{var}, and order parameter
  \code{beta}.
}
\usage{
  lbe.max(nrow, ncol, var=1, beta=1)
}
\arguments{
  \item{nrow, ncol}{the number of rows and columns in the data matrix}
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
  eigenvalue from the Laguerre Beta Ensemble so that the centered and scaled
  quantity converges in distribution to a Tracy-Widom random variable.  We
  use the second-order accurate versions of the centering and scaling defined
  in the references below.
}
\references{
  El Karoui, N. (2006). A rate of convergence result for the largest eigenvalue
    of complex white wishart matrices. \emph{Annals of Probability}. 
    \bold{34} 2077--2117.

  Johnstone, I.M. (2006). High dimensional statistical inference and random
    matrices. \emph{Proc. International Congress of Mathematicians}.
}
\author{Patrick O. Perry}
\seealso{\link{LBEMax}, \link{TracyWidom}}
\keyword{distribution}
