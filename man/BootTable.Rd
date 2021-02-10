\name{BootTable}
\alias{BootTable}
\title{
Summary table for a bootstrapped Joint Graphical Lasso model 
}
\description{
Create a table of bootstrapped means and confidence intervals for all edges of a bootstrapped Joint Graphical Lasso model obtained through \link{GroupNetworkBoot}.
}
\usage{
BootTable(BootOut)
}
\arguments{
  \item{BootOut}{
The output from \link{GroupNetworkBoot}
}
}

\details{
Summary table of the output of \link{GroupNetworkBoot}
}
\value{
\item{Var1}{Nodes included in each edge}
\item{Var2}{Nodes included in each edge}
\item{edges}{Edge identifier}
\item{sample}{sample value of each edge}
\item{boot.mean}{mean of boostrapped values of each edge}
\item{ci.lb}{lower bound of the .95 confidence interval}
\item{ci.ub}{upper bound of the .95 confidence interval}
\item{boot.zero}{
proportion of bootstraps, in which an edge was estimated as equal to zero (i.e., 0= edge not estimated as zero throughout bootstraps; 1= edge estimated as zero in all bootstraps)
}
\item{boot.pos}{Proportion of bootstraps in which an edge was estimated as >0 (i.e., positive)}
\item{boot.neg}{Proportion of bootstraps in which an edge was estimated as <0 (i.e., negative)}
\item{g}{group in which the edge was estimated}

}
\references{
Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks and their accuracy: A tutorial paper. Behavior Research Methods, 50(1), 195–212. https://doi.org/10.3758/s13428-017-0862-1
Danaher, P., Wang, P., & Witten, D. M. (2014). The joint graphical lasso for inverse covariance estimation across multiple classes. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(2), 373–397. https://doi.org/10.1111/rssb.12033
}
\author{
Nils Kappelmann <n.kappelmann@gmail.com>, Giulio Costantini
}

\seealso{
\pkg{JGL}, \pkg{qgraph}, \pkg{parcor}
}

\keyword{graphs}
\keyword{multivariate}