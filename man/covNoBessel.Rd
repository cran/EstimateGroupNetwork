\name{covNoBessel}
\alias{covNoBessel}
\title{
Covariance matrix without Bessel's correction
}
\description{
Computes the Covariance matrix without Bessel's correction, for consistency with package \pkg{JGL}
}

\usage{
covNoBessel(x,...)
}

\arguments{
  \item{x}{
A dataframe of numeric values.
}
  \item{...}{
Arguments to be passed to \link{cov}
}

}

\value{
A covariance matrix
}

\author{
Giulio Costantini
}

\examples{
library(psych)
data(bfi)
covNoBessel(bfi, use = "complete.obs")
}
\keyword{ multivariate }
